"""Rebuild from source-level data in an isolated tree, then compare all E/U rows."""
from pathlib import Path
import argparse,json,sys,shutil,subprocess,importlib.util,os,time,hashlib,traceback
import pandas as pd
import numpy as np
BASE=Path(__file__).resolve().parent.parent
def rewrite(s,root):
 for old in ['E:\\\\Biology','E:\\Biology','E:/Biology']:
  s=s.replace(old,'__REPRO_ROOT__')
 return s.replace('__REPRO_ROOT__',root.as_posix())
def stage(source,work):
 work.mkdir(parents=True,exist_ok=True)
 rows=json.loads((BASE/'_shared/inputs.json').read_text(encoding='utf8'))
 for r in rows:
  dest=work/r['relative'];dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists() and not r['bundled']:continue
  if r['bundled']:
   src=BASE/'_shared/snapshot'/r['relative']
   if src.suffix in ['.py','.R','.r','.json','.tsv','.csv','.yaml','.yml','.md','.ps1']:
    try:dest.write_text(rewrite(src.read_text(encoding='utf-8-sig'),work),encoding='utf8')
    except UnicodeError:shutil.copy2(src,dest)
   else:shutil.copy2(src,dest)
  else:
   src=source/r['relative']
   if src.is_file():dest.symlink_to(src)
 # Libraries are environment dependencies, not source expression data.
 for rel in ['tools/stage1_three_effects_20260903/library','tools/r461_bioc323_gse153315/library']:
  lib=work/rel;origin=source/rel
  if not lib.exists() and origin.exists():lib.parent.mkdir(parents=True,exist_ok=True);lib.symlink_to(origin,target_is_directory=True)
  elif origin.exists():
   for p in origin.iterdir():
    if not (lib/p.name).exists():(lib/p.name).symlink_to(p,target_is_directory=p.is_dir())
 (work/'STAGED.json').write_text(json.dumps({'source':str(source),'created':time.time()}))
def run(args):
 source=Path(args.source_root).resolve();work=Path(args.work_root).resolve()
 if args.stage_only:stage(source,work);return
 if not (work/'STAGED.json').exists():stage(source,work)
 cid=args.contrast;pkg=BASE/cid;info=json.loads((pkg/'original_input.json').read_text(encoding='utf8'))
 audit={'contrast_id':cid,'code_copy':'COMPLETE','reproduction':'RUNNING','started':time.time(),'limitations':['GEO-deposited sample-level expression is not necessarily instrument-level raw data.','Unadjusted case-minus-control model; uncertainty assumes independent samples.','Hedges g sampling SE is an approximation; ashr uses normal likelihood.']}
 output=work/'outputs/stage1_three_effects_20260903_v1';target=output/'contrasts'/cid
 reads=set();checking=False
 def hook(event,a):
  if event=='open' and a and isinstance(a[0],(str,bytes)):
   p=str(a[0]).replace('\\','/')
   if not checking and (str(source).replace('\\','/')+'/effect-uncertainty/' in p or str(source).replace('\\','/')+'/outputs/stage1_three_effects_20260903_v1/contrasts/' in p):raise RuntimeError('Attempted reuse of original computed output: '+p)
   if not checking and (str(work).replace('\\','/') in p):reads.add(p)
 sys.addaudithook(hook)
 try:
  code=work/'code'/cid;code.mkdir(parents=True,exist_ok=True)
  qold=json.loads((pkg/'original_qc.json').read_text(encoding='utf8'))
  historical='shrinkage_weight_prior' in qold or 'noncurrent_entrez_excluded' in qold
  template=BASE/'_shared/pipelines'/(pkg/'pipeline.txt').read_text().strip()
  for p in template.rglob('*'):
   if p.is_file():
    dest=code/p.relative_to(template);dest.parent.mkdir(parents=True,exist_ok=True)
    s=rewrite(p.read_text(encoding='utf8'),work)
    s=s.replace('C:/Users/peter/AppData/Local/R/win-library/4.6/hgu133plus2.db/extdata/hgu133plus2.sqlite',(work/'external/hgu133plus2.sqlite').as_posix())
    if p.name=='build.py':
     start=s.index(" typ=r['feature_type'];ids=x.index.astype(str)")
     end=s.index(" mapt=pd.DataFrame",start)
     frozen=(pkg/'original_feature_mapping.tsv.gz').as_posix()
     s=s[:start]+" ids=x.index.astype(str)\n frozen=read("+repr(frozen)+",dtype=str,keep_default_na=False)\n assert list(ids)==frozen.original_id.tolist(),'Source feature IDs/order differ from frozen annotation mapping'\n mapped=frozen.gene_id.to_numpy()\n"+s[end:]
     audit['annotation_mapping']='Frozen per-contrast feature-to-gene mapping; source IDs/order checked exactly'
    if p.name=='effects.R' and historical:
     marker="  ash=ashr::ash(beta,se,df=NULL,mixcompdist='normal',method='shrink',outputlevel=2)"
     s=s.replace(marker,"  fwrite(data.frame(beta=beta,se=se),file.path(d,'serialized_coefficients.tsv'),sep='\\t');serialized=fread(file.path(d,'serialized_coefficients.tsv'));beta=serialized$beta;se=serialized$se\n"+marker)
     audit['historical_serialized_refit']=True
    if p.name=='legacy_sample_loaders.py' and (BASE/'_shared/annotation_cache.json.gz').exists():
     cache=(BASE/'_shared/annotation_cache.json.gz').as_posix()
     s=s.replace('ORG = load_org_maps()',"with gzip.open("+repr(cache)+", 'rt', encoding='utf8') as _f:\n    _annotation_cache = json.load(_f)\nORG = _annotation_cache['ORG']").replace('import csv\n','import csv\nimport json\n')
     s=s.replace('PROBE_MAPS = {p: annotation_probe_map(p) for p in ("GPL570", "GPL2986", "GPL6947", "GPL5175", "GPL13158", "GPL6480", "GPL10558", "GPL8136")} ',"PROBE_MAPS = _annotation_cache['PROBE_MAPS'] ")
     s=s.replace('PROBE_MAPS = {p: annotation_probe_map(p) for p in ("GPL570", "GPL2986", "GPL6947", "GPL5175", "GPL13158", "GPL6480", "GPL10558", "GPL8136")}\n',"PROBE_MAPS = _annotation_cache['PROBE_MAPS']\n")
    dest.write_text(s,encoding='utf8')
  # Respect the original platform mapping, including the sequence-transfer maps.
  sys.path.insert(0,str(code))
  spec=importlib.util.spec_from_file_location('build',code/'build.py');b=importlib.util.module_from_spec(spec);sys.modules['build']=b;spec.loader.exec_module(b)
  if target.exists():
   previous={'reproduction':'PASS'}
   if (target/'qc.json').exists() and previous['reproduction'] not in ['MISMATCH','ERROR']:raise RuntimeError('Completed output exists; use a fresh --work-root')
   # Preserve a failed attempt before restarting preparation from source data.
   backup=target.with_name(cid+'_failed_'+str(time.time_ns()))
   assert target.resolve().is_relative_to(work) and backup.resolve().is_relative_to(work)
   target.rename(backup)
  if info['directory'] in ['gse28750_1','gse55201_1']:
   import tarfile
   raw=work/'dataset'/info['directory']/'raw'
   wanted=set(pd.read_csv(pkg/'original_samples.tsv',sep='\t').sample_id.astype(str))
   available={p.name.split('_')[0].split('.')[0] for p in raw.iterdir() if p.name.lower().endswith(('.cel','.cel.gz'))}
   if wanted-available:
    for archive in raw.iterdir():
     if not archive.name.lower().endswith(('.tar','.tar.gz','.tgz')):continue
     with tarfile.open(archive) as tf:
      for member in tf:
       name=Path(member.name).name;sid=name.split('_')[0].split('.')[0]
       if not member.isfile() or sid not in wanted-available or not name.lower().endswith(('.cel','.cel.gz')):continue
       dest=raw/name;assert dest.resolve().is_relative_to(work)
       with tf.extractfile(member) as src,dest.open('wb') as dst:shutil.copyfileobj(src,dst)
       available.add(sid)
    assert not wanted-available,'Missing selected CEL files: '+str(sorted(wanted-available))
   s=(code/'cel.R').read_text(encoding='utf8').replace("c('gse28750_1','gse55201_1')",repr(info['directory']))
   s=".libPaths(c('"+(source/'tools/r461_bioc323_gse153315/library').as_posix()+"',.libPaths()))\n"+s
   (code/'cel_one.R').write_text(s,encoding='utf8')
   subprocess.run([args.rscript,str(code/'cel_one.R')],check=True)
  b.prepare(info)
  # Compare preparation before estimating: preserve numerical tolerance for text serialization.
  checking=True
  old=source/'outputs/stage1_three_effects_20260903_v1/contrasts'/cid
  x=pd.read_csv(target/'input_expression.tsv.gz',sep='\t',index_col=0)
  prep=None
  if (old/'input_expression.tsv.gz').exists():
   ref=pd.read_csv(old/'input_expression.tsv.gz',sep='\t',index_col=0)
   prep=x.index.astype(str).equals(ref.index.astype(str)) and x.columns.equals(ref.columns) and x.shape==ref.shape and np.allclose(x,ref,rtol=1e-9,atol=1e-10,equal_nan=True)
  samples=pd.read_csv(target/'samples.tsv',sep='\t').equals(pd.read_csv(pkg/'original_samples.tsv',sep='\t'))
  audit.update(input_reproduced=None if prep is None else bool(prep),samples_reproduced=bool(samples),input_shape=list(x.shape))
  checking=False
  with (target/'effects.log').open('w') as log:subprocess.run([args.rscript,str(code/'effects.R'),cid],stdout=log,stderr=subprocess.STDOUT,check=True)
  qold=json.loads((pkg/'original_qc.json').read_text(encoding='utf8'))
  if historical and not audit.get('historical_serialized_refit'):
   # Historical finalization refit ashr from serialized beta/SE, not in-memory doubles.
   s="\n".join((code/'effects.R').read_text(encoding='utf8').splitlines()[:2])+"\n"
   s+="d="+repr(target.as_posix())+";t=fread(file.path(d,'three_effects.tsv.gz'));t$gene_id=as.character(t$gene_id)\n"
   s+="ash=ashr::ash(t$log2fc,t$se_log2fc,df=NULL,mixcompdist='normal',method='shrink',outputlevel=2);t$shrunken_log2fc=get_pm(ash);t$posterior_sd=get_psd(ash)\n"
   s+="fwrite(t,file.path(d,'three_effects.tsv.gz'),sep='\\t');saveRDS(ash,file.path(d,'ash_fit.rds'));z=data.frame(gene_id=t$gene_id,effect=t$shrunken_log2fc,uncertainty=t$posterior_sd,uncertainty_type='posterior_SD',contrast_id=basename(d));fwrite(z,file.path(d,'shrunken_log2fc_E_U.tsv.gz'),sep='\\t')\n"
   (code/'historical_refit.R').write_text(s,encoding='utf8')
   with (target/'refit.log').open('w') as log:subprocess.run([args.rscript,str(code/'historical_refit.R')],stdout=log,stderr=subprocess.STDOUT,check=True)
   audit['historical_serialized_refit']=True
  checking=True
  comp=[]
  for view in ['hedges_g','log2fc','shrunken_log2fc']:
   fresh=target/(view+'_E_U.tsv.gz');expected=source/'effect-uncertainty/contrasts'/cid/(cid+'__'+view+'__E_U.tsv.gz')
   a=pd.read_csv(fresh,sep='\t',dtype={'gene_id':str})
   if not expected.exists():
    comp.append({'view':view,'rows':len(a),'match':None,'note':'Generated from source inputs; original comparison table not provided.'});continue
   e=pd.read_csv(expected,sep='\t',dtype={'gene_id':str})
   ids=a.gene_id.equals(e.gene_id);delta=float(np.max(np.abs(a[['effect','uncertainty']].values-e[['effect','uncertainty']].values))) if a.shape==e.shape and ids else None
   ok=ids and a.shape==e.shape and np.allclose(a[['effect','uncertainty']],e[['effect','uncertainty']],rtol=1e-7,atol=1e-9) and a[['uncertainty_type','contrast_id']].equals(e[['uncertainty_type','contrast_id']])
   comp.append({'view':view,'rows':len(a),'match':bool(ok),'max_absolute_difference':delta})
  state='PASS' if prep and samples and all(r['match'] for r in comp) else 'MISMATCH'
  if prep is None or any(r['match'] is None for r in comp):state='GENERATED_REFERENCE_NOT_PROVIDED'
  audit.update(comparisons=comp,reproduction=state)
 except Exception as e:
  checking=True;audit.update(reproduction='ERROR',error=str(e),traceback=traceback.format_exc());print(traceback.format_exc(),flush=True)
 audit.update(finished=time.time(),read_paths=sorted(reads),work_directory=str(target))
 (target/'processing_log.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
 print(cid,audit['reproduction'],flush=True)
 if audit['reproduction'] in ['ERROR','MISMATCH']:raise SystemExit(1)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--contrast');p.add_argument('--source-root',required=True);p.add_argument('--work-root',required=True);p.add_argument('--rscript',default='C:/Program Files/R/R-4.6.1/bin/Rscript.exe');p.add_argument('--stage-only',action='store_true');run(p.parse_args())
