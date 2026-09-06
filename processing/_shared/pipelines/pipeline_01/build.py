"""Versioned sample-level reconstruction. Never executes an old result-writing tail."""
from pathlib import Path
import ast, json, sys, importlib.util, gzip, re, traceback, hashlib, subprocess, socket
socket.setdefaulttimeout(30)
import numpy as np
import pandas as pd
ROOT=Path('E:/Biology'); HERE=Path(__file__).parent
OUT=ROOT/'outputs/stage1_three_effects_20260903_v1'; OUT.mkdir(parents=True,exist_ok=True)
AUD=ROOT/'outputs/dataset_contrast_audit_20260903_v1'
def loadmod(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
legacy=loadmod('stage1_legacy',HERE/'dependencies/legacy_sample_loaders.py')
repair=loadmod('stage1_repair',HERE/'dependencies/legacy_gene_mapping.py')
INV=json.loads((AUD/'audit_data.json').read_text(encoding='utf8'))['inventory']
EVID={x['directory']:x for x in json.loads((AUD/'evidence_snapshot.json').read_text(encoding='utf8'))}
def read(p,**kw):return pd.read_csv(p,sep='\t',**kw)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def ensure_platform(platform,directory=None):
 if not platform or legacy.PROBE_MAPS.get(platform):return
 annout=OUT/'annotations';annout.mkdir(exist_ok=True)
 candidates=list((ROOT/'dataset').glob('*/raw/'+platform+'*'))+list((ROOT/'tmp/eu58_full_model_20260827/annotations').glob(platform+'*'))
 if directory:candidates+=list((ROOT/'dataset'/directory/'raw').glob('*family.soft.gz'))
 table=None
 for p in candidates:
  try:table=legacy.read_annotation(p);break
  except Exception:pass
 if table is None:
  import urllib.request
  p=annout/(platform+'_platform_only.txt')
  if not p.exists():
   url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+platform+'&targ=self&form=text&view=full'
   with urllib.request.urlopen(url,timeout=30) as response:
    data=response.read(100*1024*1024+1)
   if len(data)>100*1024*1024:raise ValueError('Platform-only response exceeds size limit')
   p.write_bytes(data)
  table=legacy.read_annotation(p)
 cols={c.lower().replace('_',' ').strip():c for c in table.columns};ic=cols.get('id',table.columns[0]);gc=next((cols[k] for k in ['entrez gene id','entrez gene','gene id','entrezgeneid','entrez id','locuslink id'] if k in cols),None);sc=next((cols[k] for k in ['gene symbol','symbol','gene name','genename','gene assignment'] if k in cols),None)
 out={};valid=set(legacy.ORG['valid_entrez'])
 accession_genes={}
 if 'gb list' in cols or 'gb acc' in cols:
  import sqlite3
  with sqlite3.connect(legacy.ORG_DB) as con:
   for accession,gene in con.execute('select a.accession,g.gene_id from accessions a join genes g on a._id=g._id'):
    accession_genes.setdefault(str(accession).split('.')[0],set()).add(str(gene))
  if platform=='GPL4302':
   with sqlite3.connect('C:/Users/peter/AppData/Local/R/win-library/4.6/hgu133plus2.db/extdata/hgu133plus2.sqlite') as con:
    for accession,gene in con.execute('select a.accession,p.gene_id from accessions a join probes p on a.probe_id=p.probe_id where p.gene_id is not null'):
     if gene in valid:accession_genes.setdefault(str(accession).split('.')[0],set()).add(str(gene))
 sequence_genes={}
 if platform=='GPL32193':
  # Custom Agilent layout has no gene annotation. Transfer only exact probe sequences
  # with a unique Entrez assignment from two deposited human Agilent platforms.
  for pl,ap in [('GPL6480',ROOT/'tmp/eu58_full_model_20260827/annotations/GPL6480.annot.gz'),('GPL21185',ROOT/'dataset/gse124272_1/raw/GPL21185_platform.txt')]:
   if pl not in legacy.PROBE_MAPS:
    mp=read(annout/(pl+'_mapping.tsv'),dtype=str);legacy.PROBE_MAPS[pl]=dict(zip(mp.probe_id,mp.gene_id))
   at=legacy.read_annotation(ap);sq=next(c for c in at.columns if c.upper() in ['SEQUENCE','PLATFORM_SEQUENCE'])
   for record in at.to_dict('records'):
    gene=legacy.PROBE_MAPS[pl].get(str(record['ID']));seq=str(record[sq]).strip().upper()
    if gene and len(seq)>=40 and re.fullmatch('[ACGT]+',seq):sequence_genes.setdefault(seq,set()).add(gene)
 for rr in table.to_dict('records'):
  gs=set(re.findall(r'\d+',str(rr[gc]))) if gc else set()
  gs &= valid
  if not gs and 'gene assignment' in cols:
   for part in str(rr[cols['gene assignment']]).split(' /// '):
    token=part.split(' // ')[-1].strip()
    if token in valid:gs.add(token)
  if not gs and sc:
   syms=re.split(r'\s*///\s*|\s*//\s*|[;,]',str(rr[sc]));gs={legacy.ORG['symbol'].get(s.strip().upper(),'') for s in syms}-{''}
  if not gs and accession_genes:
   for key in ['gb list','gb acc']:
    if key in cols:
     for acc in re.split(r'[;,\s/]+',str(rr[cols[key]])):gs.update(accession_genes.get(acc.split('.')[0],set()))
  if not gs and sequence_genes:gs=sequence_genes.get(str(rr.get(cols.get('sequence',''),'')).strip().upper(),set())
  if len(gs)==1:out[str(rr[ic])]=next(iter(gs))
 legacy.PROBE_MAPS[platform]=out
 pd.DataFrame(list(out.items()),columns=['probe_id','gene_id']).to_csv(annout/(platform+'_mapping.tsv'),sep='\t',index=False)
def capture(script):
 tree=ast.parse(script.read_text(encoding='utf8'));main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
 body=[];found=[]
 for node in main.body:
  txt=ast.unparse(node)
  if isinstance(node,ast.Expr) and any(t in txt for t in ['.to_csv(','.mkdir(','.write_text(','print(']):continue
  if isinstance(node,ast.Assign):
   # Preserve pre-normalization counts; missing entries are never filled with zero.
   if 'np.log2(' in txt and '.div(' in txt:
    call=next(n for n in ast.walk(node) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='div')
    body+=ast.parse('__raw_counts = '+ast.unparse(call.func.value)+'.copy()').body
   targets=[n.id for t in node.targets for n in ast.walk(t) if isinstance(n,ast.Name)]
   if '.to_numpy(' in txt and any(t in ['A','B','a','b','case','ctrl'] for t in targets):
    calls=[n for n in ast.walk(node.value) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='to_numpy']
    for call in calls:
     if len(found)<2:
      nm='__case_frame' if not found else '__control_frame';found.append(nm)
      body+=ast.parse(nm+' = '+ast.unparse(call.func.value)).body
    if len(found)==2:
     body+=ast.parse('return locals()').body;break
  body.append(node)
 if len(found)!=2:raise ValueError('native capture needs explicit case/control matrices')
 main.body=body
 # Remove source writes in helpers as well, and stop silent missing-to-zero coercions.
 class Safe(ast.NodeTransformer):
  def visit_Expr(self,n):
   if isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute) and n.value.func.attr in ['to_csv','write_text','write_bytes','mkdir']:return None
   return self.generic_visit(n)
  def visit_Call(self,n):
   n=self.generic_visit(n)
   if isinstance(n.func,ast.Attribute) and n.func.attr=='fillna' and n.args and isinstance(n.args[0],ast.Constant) and n.args[0].value==0:return n.func.value
   return n
 tree=Safe().visit(tree);tree.body=[n for n in tree.body if not isinstance(n,ast.If)];ast.fix_missing_locations(tree)
 ns={'__file__':str(script),'__name__':'stage1_capture'};exec(compile(tree,str(script),'exec'),ns)
 return ns['main'](),ast.unparse(tree)
def native(row):
 d=ROOT/'dataset'/row['directory'];local,code=capture(d/'process/run_analysis.py')
 a=local['__case_frame'];b=local['__control_frame']
 if not isinstance(a,pd.DataFrame) or not isinstance(b,pd.DataFrame):raise ValueError('captured selection is not a DataFrame')
 x=pd.concat([a,b],axis=1);kind='normalized';src='native log2/normalized expression'
 raw=local.get('__raw_counts')
 if isinstance(raw,pd.DataFrame) and set(x.columns)<=set(raw.columns):
  raw=raw[x.columns].copy()
  if isinstance(raw.index,pd.RangeIndex):
   gene=local.get('gene')
   if gene is None or len(gene)!=len(raw):raise ValueError('raw counts lost feature IDs')
   raw.index=np.asarray(gene)
  x=raw;kind='counts';src='native raw counts intercepted before normalization'
 if d.name=='gse34404_1':src='non_normalized Illumina signal; native sample selection'
 platform={'gse17048_1':'GPL6947','gse48060_1':'GPL570','gse137340_1':'GPL10558','gse51404_1':'GPL6947','gse51405_2':'GPL6947','gse54248_1':'GPL10558','gse94648_1':'GPL19109'}.get(row['directory'])
 if not platform:
  alltext=' '.join(EVID[row['directory']]['geo_headers'])+' '+code
  plats=re.findall(r'GPL\d+',alltext);platform=plats[0] if plats else None
 return dict(x=x,cases=list(a.columns),controls=list(b.columns),kind=kind,platform=platform,feature_type='auto',source=src,adapter_code=code)
def legacy_capture(eu,frame,cases,controls,feature_col,feature_type,kind,source,platform=None,symbol_col=None,feature_space=None):
 return dict(x=frame.set_index(feature_col)[cases+controls],cases=cases,controls=controls,kind=kind,source=source,platform=platform,feature_type=feature_type)
legacy.build_eu=legacy_capture
def special(row):
 d=row['directory'];base=ROOT/'dataset'/d;acc=row['accession']
 if d in ['gse44314_1','gse94648_1','gse42825_1']:
  from extra_loaders import extra
  return extra(row,sys.modules[__name__])
 if row['contrast_id']=='GSE69683_severe_asthma_vs_HC':
  old=EVID[d]['prior57'];p=ROOT/'dataset'/d/'raw/GSE69683_series_matrix.txt.gz';f,_=legacy.read_geo_series_matrix(p);m=read(Path(old['证据路径']).parent/'sample_manifest.tsv')
  return legacy_capture('',f,m.loc[m.group.eq('severe_asthma'),'gsm'].tolist(),m.loc[m.group.eq('control'),'gsm'].tolist(),'feature_id','probe','normalized',str(p),'GPL13158')
 if row['contrast_id']=='GSE68004_incomplete_KD_vs_HC':
  old=EVID[d]['prior57'];p=ROOT/'dataset'/d/'raw/GSE68004_family.soft.gz';f=legacy.read_family_soft_values(p);m=read(Path(old['证据路径']).parent/'sample_manifest.tsv');r=m[m.contrast.eq('incomplete_KD_vs_healthy')].iloc[0]
  return legacy_capture('',f,r.case_gsm_list.split(';'),r.control_gsm_list.split(';'),'feature_id','probe','normalized',str(p),'GPL10558')
 if d.startswith('gse112057_'):
  m=read(base/'process/sample_metadata.tsv');target={1:"Crohn's Disease",2:'Ulcerative Colitis',3:'Systemic JIA',4:'Oligoarticular JIA',5:'Polyarticular JIA'}[int(d[-1])]
  m['col']=m.title.str.split('_',n=1).str[0];p=base/'raw/GSE112057_RawCounts_dataset.txt.gz';f=read(p);f.iloc[:,0]=f.iloc[:,0].astype(str).str.replace(r'_\d+$','',regex=True)
  return legacy_capture('',f,m.loc[m.diagnosis.eq(target),'col'].tolist(),m.loc[m.diagnosis.eq('Control'),'col'].tolist(),f.columns[0],'symbol','counts',str(p))
 if d.startswith('gse186507_'):
  base=ROOT/'dataset/gse186507_1';m=read(base/'process/sample_manifest.tsv');p=base/'raw/GSE186507_MSCCR_Blood_counts.txt.gz'
  with gzip.open(p,'rt') as f:names=[t.strip('"') for t in f.readline().split()]
  cases=m.loc[m.group.eq('CD' if d.endswith('_1') else 'UC'),'subject'].tolist();controls=m.loc[m.group.eq('Healthy_Control'),'subject'].tolist()
  f=pd.read_csv(p,sep=r'\s+',header=None,skiprows=1,names=['gene_id']+names,usecols=['gene_id']+cases+controls)
  return legacy_capture('',f,cases,controls,'gene_id','ensembl','counts',str(p))
 if d=='gse181228_1':
  m=read(base/'process/sample_manifest.tsv');m=m[m.inclusion.eq('include')];cases=m.loc[m.group.eq('case') & m.characteristics.str.contains('visit: Baseline',regex=False),'gsm'].tolist();ctrl=m.loc[m.group.eq('healthy_control'),'gsm'].tolist();p=base/'raw/counts.tsv.gz';f=read(p)
  return legacy_capture('',f,cases,ctrl,f.columns[0],'entrez','counts',str(p))
 return None
def prepare(row):
 d=row['directory'];target=OUT/'contrasts'/row['contrast_id'];target.mkdir(parents=True,exist_ok=True)
 if (target/'input.json').exists():
  info=json.loads((target/'input.json').read_text(encoding='utf8'));info.update(row)
  (target/'input.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf8');return info
 r=special(row)
 if r is None:
  try:r=native(row)
  except Exception as ne:
   old=EVID[d]['prior57']
   if old:r=legacy.load_missing_eu(pd.Series(old))
   else:
    from extra_loaders import extra
    r=extra(row,sys.modules[__name__])
    if r is None:raise ValueError(str(ne)+'; no additional loader')
 x=r['x'].apply(pd.to_numeric,errors='coerce');cases=list(map(str,r['cases']));ctrl=list(map(str,r['controls']));x.columns=x.columns.astype(str)
 assert len(cases)>=2 and len(ctrl)>=2 and not set(cases)&set(ctrl) and len(set(cases+ctrl))==len(cases+ctrl),'invalid/overlapping samples'
 for k,n in [('n_case_audit',len(cases)),('n_control_audit',len(ctrl))]:
  if row.get(k) is not None:assert int(row[k])==n,f'{k}: {n} != {row[k]}'
 typ=r['feature_type'];ids=x.index.astype(str)
 if typ=='auto':
  mapped=repair.auto_map_features(ids,legacy)[0]
  if (mapped!='').mean()<0.5 and r.get('platform'):
   ensure_platform(r['platform'],d);mapped=np.array([legacy.map_feature(v,'probe',r['platform']) or '' for v in ids])
 else:
  if typ=='probe':ensure_platform(r.get('platform'),d)
  mapped=np.array([legacy.map_feature(v,typ,r.get('platform')) or '' for v in ids])
 mapt=pd.DataFrame({'original_id':ids,'gene_id':mapped});mapt.to_csv(target/'feature_mapping.tsv.gz',sep='\t',index=False)
 x.index=mapped;x=x[x.index!=''];pre=len(x)
 # All feature collapse precedes effects. Counts sum; normalized probes average.
 x=x.groupby(level=0).sum(min_count=1) if r['kind']=='counts' else x.groupby(level=0).mean()
 x.index.name='gene_id';x=x[cases+ctrl];missing=int(x.isna().sum().sum());x=x.replace([np.inf,-np.inf],np.nan).dropna()
 if r['kind']=='counts':assert (x.to_numpy()>=0).all(),'negative counts'
 assert len(x)>=1000,f'Only {len(x)} mapped complete genes'
 x.to_csv(target/'input_expression.tsv.gz',sep='\t',float_format='%.10g')
 sm=pd.DataFrame({'sample_id':cases+ctrl,'group':['case']*len(cases)+['control']*len(ctrl),'accession':row['accession'],'study_group':row['group']});sm.to_csv(target/'samples.tsv',sep='\t',index=False)
 if r.get('adapter_code'):(target/'sample_loader_snapshot.py').write_text(r['adapter_code'],encoding='utf8')
 info={**row,'n_case':len(cases),'n_control':len(ctrl),'input_kind':r['kind'],'source':r['source'],'platform':r.get('platform'),'input_genes':len(x),'mapped_rows_before_collapse':pre,'missing_values_removed_with_genes':missing,'input_sha256':sha(target/'input_expression.tsv.gz'),'scope':'core' if row['status'][0] in 'AB' else 'supplemental'}
 (target/'input.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf8');return info
def main():
 import sqlite3
 universe=OUT/'canonical_gene_universe.tsv'
 if not universe.exists():
  with sqlite3.connect(legacy.ORG_DB) as con:
   genes=pd.read_sql_query('select g.gene_id,i.symbol from genes g join gene_info i on g._id=i._id',con)
  genes.sort_values('gene_id',key=lambda s:s.astype(int)).to_csv(universe,sep='\t',index=False)
 selected=[dict(r) for r in INV if r['status'][0] in 'ABC']
 from extra_loaders import amend
 selected=amend(selected)
 (OUT/'planned_contrasts.json').write_text(json.dumps(selected,ensure_ascii=False,indent=2),encoding='utf8')
 filt=set(sys.argv[1:]);results=[]
 for i,r in enumerate(selected):
  if filt and r['directory'] not in filt and r['contrast_id'] not in filt:continue
  try:
   info=prepare(r);results.append({'contrast_id':r['contrast_id'],'status':'INPUT_READY','genes':info['input_genes']});print(f'{i+1}/{len(selected)} READY {r["contrast_id"]} {info["input_genes"]}',flush=True)
  except Exception as e:
   results.append({'contrast_id':r['contrast_id'],'directory':r['directory'],'status':'ERROR','error':str(e)});print(f'{i+1}/{len(selected)} ERROR {r["contrast_id"]} {e}',flush=True);traceback.print_exc()
  (OUT/'last_preparation_run.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf8')
if __name__=='__main__':main()
