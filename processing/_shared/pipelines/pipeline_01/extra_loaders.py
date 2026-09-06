import re,gzip,json
import numpy as np
import pandas as pd
def amend(rows):
 added=[]
 for r in rows:
  if r['directory']=='gse181228_1':r.update(contrast_id='GSE181228_SScILD_baseline_vs_HC',disease='系统性硬化症相关间质性肺病（基线）',n_case_audit=132,n_control_audit=43)
  if r['directory']=='gse94648_1':
   r.update(contrast_id='GSE94648_CD_unique_vs_HC',disease='克罗恩病',n_case_audit=48,n_control_audit=20)
   added.append(dict(r,contrast_id='GSE94648_UC_unique_vs_HC',disease='溃疡性结肠炎',n_case_audit=25))
  if r['directory']=='gse44314_1':
   r.update(contrast_id='GSE44314_classic1A_T1D_vs_HC',disease='经典1A型糖尿病',n_case_audit=5,n_control_audit=6)
   added.append(dict(r,contrast_id='GSE44314_fulminant_T1D_vs_HC',disease='暴发型1型糖尿病',status='C_supplemental_subtype',domain='',n_case_audit=5))
  for source,cid,n,domain in [('gse42825_1','GSE42825_TB_vs_HC',8,'感染性疾病'),('gse68004_1','GSE68004_incomplete_KD_vs_HC',13,'系统性自身免疫与风湿性疾病'),('gse69683_1','GSE69683_severe_asthma_vs_HC',334,'慢性呼吸系统疾病')]:
   if r['directory']==source:added.append(dict(r,contrast_id=cid,disease={'gse42825_1':'肺结核','gse68004_1':'不完全型川崎病','gse69683_1':'重度哮喘'}[source],n_case_audit=n,domain=domain))
 return rows+added
def extra(row,api):
 root=api.ROOT;d=row['directory'];acc=row['accession'];base=root/'dataset'/d;raw=base/'raw';cid=row['contrast_id'];read=api.read
 freeze=root/'census_step0/candidate_extensions_2026-08-22_v1'/(acc+'_lock_2026-08-22_v1')
 cand=root/'phase1b_EU_coordination_v1_2026-08-23/01_candidates'/acc
 def finish(f,ca,co,kind='normalized',typ='probe',platform=None,source=''):
  return api.legacy_capture(cid,f,list(ca),list(co),f.columns[0],typ,kind,str(source),platform)
 def series(p):
  f,meta=api.legacy.read_geo_series_matrix(p);plats=meta.get('!Sample_platform_id',[[]])[0];return f,meta,plats[0] if plats else None
 if d in ['gse201332_1','gse42331_1','gse5812_1','gse44314_1','gse69063_1','gse94648_1']:
  m=read(base/'process/sample_manifest.tsv');p=next(raw.glob('*series_matrix*'));f,meta,pl=series(p)
  if d=='gse44314_1':
   ca=m.loc[m.title.str.contains('Fulminant' if 'fulminant' in cid else 'Classical 1A',case=False),'sample_id'];co=m.loc[m.role.eq('control'),'sample_id']
  elif d=='gse94648_1':
   m=m[~m.title.str.contains('bis',case=False)];ca=m.loc[m.title.str.contains('CD' if '_CD_' in cid else 'UC'),'gsm'];co=m.loc[m.group.eq('healthy_control'),'gsm']
  elif d=='gse69063_1':
   ca=m.loc[m.role.eq('case') & m.included.eq(1),'sample_id'];co=m.loc[m.role.eq('control') & m.included.eq(1),'sample_id']
  else:ca=m.loc[m.group.eq('case'),'gsm'];co=m.loc[m.group.eq('control'),'gsm']
  return finish(f,ca,co,platform=pl,source=p)
 if d=='gse42825_1':
  m=read(base/'process/sample_manifest.tsv');sm=api.loadmod('gse42825_src',base/'process/run_analysis.py');ca=sm.CASE;co=sm.CONTROL
  if '_TB_' in cid:
   cm=read(cand/'contrast_manifest.tsv');rr=cm[cm.contrast.str.contains('TB',case=False)].iloc[0];ca=str(rr.case_gsm_list).split(';')
  mp=dict(zip(m.gsm,m.matrix_column));p=raw/'non_normalized_signal.tsv.gz';return finish(read(p),[mp[s] for s in sorted(ca)],[mp[s] for s in sorted(co)],platform='GPL10558',source=p)
 if d in ['gse136371_1','gse185855_1','gse270454_1','gse270454_2']:
  if d=='gse136371_1':
   m=read(freeze/(acc+'_sample_metadata.tsv'));ca=m.loc[m.diagnosis.ne('healthy control'),'gsm'];co=m.loc[m.diagnosis.eq('healthy control'),'gsm'];p=raw/(acc+'_counts_matrix.tsv.gz');typ='symbol'
  elif d=='gse185855_1':
   m=read(freeze/(acc+'_baseline_sample_lock.tsv'));ca=m.loc[m.group.eq('MDD'),'gsm'];co=m.loc[m.group.ne('MDD'),'gsm'];p=raw/(acc+'_baseline_counts_matrix.tsv.gz');typ='ensembl'
  else:
   m=read(freeze/(acc+'_sample_subject_data_lock_v1.tsv'));m=m[m.modality.eq('RNA-seq')].drop_duplicates('sample_label');ca=m.loc[m.group.eq('AD' if d.endswith('_1') else 'MCI'),'sample_label'];co=m.loc[m.group.eq('ASO'),'sample_label'];p=raw/(acc+'_RNAseq-combined-counts-matrix.csv.gz');typ='symbol'
  f=pd.read_csv(p,sep=',' if '270454' in d else '\t');f.columns=f.columns.str.strip()
  if d=='gse136371_1':f.iloc[:,0]=f.iloc[:,0].str.replace(r'_\d+$','',regex=True)
  return finish(f,ca,co,'counts',typ,source=p)
 if d.startswith('gse18781_') or d.startswith('gse61281_'):
  m=read(freeze/(acc+'_sample_subject_data_lock_v1.tsv'));p=next(raw.glob('*series_matrix*'));f,meta,pl=series(p)
  if acc=='GSE18781':m=m[m.contrast_id.eq(cid)];ca=m.loc[m.disease_group.eq('case'),'sample_id'];co=m.loc[m.disease_group.eq('healthy_control'),'sample_id']
  else:ca=m.loc[m.disease_group.eq('PsA' if d.endswith('_1') else 'PsC'),'sample_id'];co=m.loc[m.disease_group.eq('HC'),'sample_id']
  return finish(f,ca,co,platform=pl,source=p)
 if d.startswith('gse196822_'):
  m=read(root/'dataset/gse196822_1/process/sample_manifest.tsv');target={2:'Asymptomatic',3:'Covid-19 Bacterial Coinfection',4:'Mild',5:'Severe'}[int(d[-1])];p=raw/'GSE196822_Raw_counts_matrix.csv.gz';f=pd.read_csv(p)
  return finish(f,m.loc[m.condition.str.startswith(target),'sequencing_id'],m.loc[m.condition.eq('Healthy'),'sequencing_id'],'counts','ensembl',source=p)
 if d=='gse177044_1':
  m=read(cand/'sample_manifest.tsv');fs=[pd.read_csv(p).set_index('Geneid').drop(columns=['gene_name'],errors='ignore') for p in sorted(raw.glob('*.csv.gz'))];f=pd.concat(fs,axis=1);assert not f.columns.duplicated().any()
  return finish(f.reset_index(),m.loc[m.disease.eq('UC'),'title'],m.loc[m.disease.eq('Control'),'title'],'counts','ensembl',source=';'.join(str(p) for p in raw.glob('*.csv.gz')))
 if d=='gse244915_1':
  p=raw/'GSE244915_RNAseq_count_matrix.txt.gz'
  with gzip.open(p,'rt') as h:cols=h.readline().split()
  f=pd.read_csv(p,sep='\t',skiprows=1,names=cols);return finish(f,f.columns[4:],f.columns[1:4],'counts','symbol',source=p)
 if d=='gse190125_1':
  p=raw/'GSE190125_Counts_for_GEO.txt.gz';m=read(base/'process/sample_manifest.tsv');arrays={};order=None
  for chunk in pd.read_csv(p,sep='\t',chunksize=200000):
   for sid,part in chunk.groupby(chunk.columns[0],sort=False):
    ser=pd.Series(pd.to_numeric(part.iloc[:,-1]).values,index=part.iloc[:,1].astype(str));arrays.setdefault(sid,[]).append(ser)
  f=pd.DataFrame({sid:pd.concat(parts).groupby(level=0).sum(min_count=1) for sid,parts in arrays.items()}).reset_index(names='gene_id')
  return finish(f,m.loc[m.primary_role.eq('case'),'sample_id'],m.loc[m.primary_role.eq('control'),'sample_id'],'counts','ensembl',source=p)
 if d=='gse130953_1':
  m=read(base/'process/sample_manifest.tsv');p=raw/'GSE130953_Non-normalized_data.txt.gz';f=api.legacy.extract_signal_columns(read(p));cols={str(c).replace('.AVG_Signal',''):c for c in f.columns}
  return finish(f,[cols[s] for s in m.loc[m.group.eq('SSc'),'array_id']],[cols[s] for s in m.loc[m.group.eq('Healthy'),'array_id']],platform='GPL10558',source=p)
 if d in ['gse38267_1','gse124272_1']:
  import tarfile
  sm=api.loadmod('raw_'+d,base/'process/run_analysis.py');m=read(base/'process/sample_manifest.tsv')
  ca=m.loc[m.group.eq('case' if d=='gse38267_1' else 'IDD'),'gsm' if d=='gse38267_1' else 'geo_accession'].tolist();co=m.loc[m.group.ne('case' if d=='gse38267_1' else 'IDD'),'gsm' if d=='gse38267_1' else 'geo_accession'].tolist()
  with tarfile.open(raw/'raw.tar') as tar:
   names=tar.getnames();arrays={s:pd.Series(sm.read_agilent(tar,next(n for n in names if n.startswith(s+'_')))) for s in ca+co}
  f=pd.DataFrame(arrays).reset_index(names='probe_id');return finish(f,ca,co,platform='GPL13607' if d=='gse38267_1' else 'GPL21185',source='non_normalized '+str(raw/'raw.tar'))
 if d in ['gse28750_1','gse55201_1']:
  p=api.OUT/'cel_normalized'/d/'rma.tsv.gz';f=read(p);m=read(base/'process/sample_manifest.tsv')
  if d=='gse28750_1':m=m[m.include_locked.astype(str).str.upper().eq('TRUE')];key='sample_id'
  else:key='sample_key'
  return finish(f,m.loc[m.disease_group.eq('case'),key],m.loc[m.disease_group.eq('healthy_control'),key],platform='GPL570',source='Affymetrix RMA from '+str(raw))
 if d=='gse154851_1':
  p=next(raw.glob('*series_matrix*'));f,meta,pl=series(p);ids=meta['!Sample_geo_accession'][0];chars=meta['!Sample_characteristics_ch1'];text=[' | '.join(v[i] for v in chars) for i in range(len(ids))]
  titles=meta['!Sample_title'][0];ca=[s for s,t in zip(ids,titles) if t.startswith('Patient ')];co=[s for s,t in zip(ids,titles) if t.startswith('Control ')];return finish(f,ca,co,platform=pl,source=p)
 if d=='gse160207_1':
  p=raw/'counts.txt.gz';m=read(base/'process/sample_manifest.tsv');f=read(p).drop(columns=['Unnamed: 0','Chr']);return finish(f,m.loc[m.group.eq('case'),'sample_name'],m.loc[m.group.eq('control'),'sample_name'],'counts','symbol',source=p)
 if d=='gse184876_1':
  p=raw/'counts.txt.gz';m=read(base/'process/sample_manifest.tsv');f=read(p);return finish(f,m.loc[m.group.eq('case'),'count_column'],m.loc[m.group.eq('control'),'count_column'],'counts','ensembl',source=p)
 return None
