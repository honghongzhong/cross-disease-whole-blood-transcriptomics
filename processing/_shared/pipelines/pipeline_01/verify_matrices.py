from pathlib import Path
import json
import pandas as pd,numpy as np
OUT=Path('E:/Biology/outputs/stage1_three_effects_20260903_v1')
summary=json.loads((OUT/'verification_summary.json').read_text(encoding='utf8'));assert summary['status']=='COMPLETE'
views={'hedges_g':('hedges_g','se_hedges_g'),'log2fc':('log2fc','se_log2fc'),'shrunken_log2fc':('shrunken_log2fc','posterior_sd')}
checks=[]
for scope in ['core','supplemental','all']:
 d=OUT/'matrices'/scope;mask=pd.read_csv(d/'observed_mask.tsv.gz',sep='\t',index_col=0);mask.index=mask.index.astype(str)
 genes=pd.read_csv(d/'gene_index.tsv',sep='\t',dtype=str).gene_id.tolist();ids=pd.read_csv(d/'contrast_index.tsv',sep='\t').contrast_id.tolist()
 assert list(mask.index)==genes and list(mask.columns)==ids
 tables={cid:pd.read_csv(OUT/'contrasts'/cid/'three_effects.tsv.gz',sep='\t',dtype={'gene_id':str}).set_index('gene_id').reindex(genes) for cid in ids}
 for name,(e,u) in views.items():
  for kind,col in [('E',e),('U',u)]:
   p=d/(name+'_'+kind+'.tsv.gz');x=pd.read_csv(p,sep='\t',index_col=0);x.index=x.index.astype(str)
   assert x.index.equals(mask.index) and x.columns.equals(mask.columns),p
   assert np.array_equal(np.isfinite(x.to_numpy()),mask.to_numpy().astype(bool)),p
   if kind=='U':assert (x.to_numpy()[mask.to_numpy().astype(bool)]>0).all(),p
   for cid in ids:assert np.allclose(x[cid],tables[cid][col],equal_nan=True,rtol=1e-11,atol=1e-12),(p,cid)
   checks.append({'scope':scope,'view':name,'matrix':kind,'genes':len(genes),'contrasts':len(ids),'observed_entries':int(mask.to_numpy().sum()),'status':'PASS'})
 print('MATRICES_VERIFIED',scope,mask.shape,flush=True)
pd.DataFrame(checks).to_csv(OUT/'matrix_roundtrip_validation.tsv',sep='\t',index=False)
summary['matrix_roundtrip_validation']={'status':'PASS','matrices_checked':len(checks)}
(OUT/'verification_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
print('ALL_18_MATRICES_PASS',flush=True)
