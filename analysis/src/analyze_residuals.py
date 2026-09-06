"""Descriptive held-study residual analysis; no causal or individual inference."""
from pathlib import Path
import os
import numpy as np,pandas as pd
from scipy.stats import spearmanr
P=Path(os.environ['BLOOD_OUTPUT']); A=P/'analysis'; R=Path(os.environ['BLOOD_DATA_ROOT'])/'contrasts'
out=P/'analysis';out.mkdir(exist_ok=True)
te=pd.read_csv(A/'test_manifest.tsv',sep='\t');rows=[];summary=[]
sets={p[0]:set(p[2:]) for p in [line.split('\t') for line in (P/'resources/h.all.v2025.1.Hs.entrez.gmt').read_text().splitlines()] if len(p)>2}
for eff,k in [('hedges_g',5),('log2fc',4),('shrunken_log2fc',3)]:
 z=np.load(A/f'{eff}_reference.npz');genes=z['genes'];q=z['L'][:,:k]
 for _,t in te.iterrows():
  cid=t.contrast_id;f=pd.read_csv(R/cid/f'{cid}__{eff}__E_U.tsv.gz',sep='\t',dtype={'gene_id':str}).set_index('gene_id').reindex(genes)
  y=f.effect.to_numpy();u=f.uncertainty.to_numpy();ok=np.isfinite(y);v=y[ok];l=q[ok];p=l@np.linalg.lstsq(l,v,rcond=None)[0];res=v-p;gg=genes[ok]
  summary.append(dict(effect_type=eff,contrast_id=cid,n_case=t.n_case,n_control=t.n_control,effective_n=t.effective_n,R2=1-np.sum(res**2)/np.sum((v-v.mean())**2),residual_energy_fraction=np.sum(res**2)/np.sum(v**2),median_uncertainty=np.nanmedian(u),gene_residual_uncertainty_rho=spearmanr(np.abs(res),u[ok]).statistic))
  for name,gs in sets.items():
   ix=np.isin(gg,list(gs));energy=np.sum(v[ix]**2)
   rows.append(dict(effect_type=eff,contrast_id=cid,pathway=name,n_genes=int(ix.sum()),observed_rms=np.sqrt(np.mean(v[ix]**2)),residual_rms=np.sqrt(np.mean(res[ix]**2)),residual_energy_fraction=np.sum(res[ix]**2)/energy if energy>0 else np.nan,mean_observed=v[ix].mean(),mean_residual=res[ix].mean()))
pd.DataFrame(rows).to_csv(out/'pathway_residuals.tsv',sep='\t',index=False)
s=pd.DataFrame(summary);s.to_csv(out/'residual_summary.tsv',sep='\t',index=False)
for eff,d in s.groupby('effect_type'):
 print(eff,'effectiveN vs R2',spearmanr(d.effective_n,d.R2))
print(pd.DataFrame(rows).query("effect_type == 'hedges_g'").groupby('pathway').residual_energy_fraction.median().sort_values(ascending=False).head(10).to_string())
print(s[['effect_type','contrast_id','R2','effective_n']].to_string(index=False))
