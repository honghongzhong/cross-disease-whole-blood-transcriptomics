from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw';OUT=HERE.parent/'result'
def main():
 m=pd.read_csv(HERE/'sample_manifest.tsv',sep='\t');m=m[m.include_locked.astype(str).str.upper().eq('TRUE')];x=pd.read_csv(RAW/'GSE123658_read_counts.gene_level.txt.gz',sep='\t',compression='gzip').set_index('Samples');ids=m.sample_key.tolist();x=x[ids].apply(pd.to_numeric,errors='coerce').fillna(0);x=np.log2(x.div(x.sum(),axis=1)*1e6+1);ca=m.loc[m.disease_group.eq('case'),'sample_key'].tolist();co=m.loc[m.disease_group.eq('healthy_control'),'sample_key'].tolist();A=x[ca].to_numpy();B=x[co].to_numpy();e=A.mean(1)-B.mean(1);u=np.sqrt(A.var(1,ddof=1)/len(ca)+B.var(1,ddof=1)/len(co));p=ttest_ind(A,B,axis=1,equal_var=False,nan_policy='omit').pvalue;o=np.argsort(np.nan_to_num(p,nan=1));q=np.empty_like(p);q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=pd.DataFrame({'gene_id':x.index,'E':e,'U':u,'p_value':p,'q_value':np.clip(q,0,1)});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].sort_values('gene_id');OUT.mkdir(exist_ok=True);path=OUT/'GSE123658_T1D_vs_HC_E_U.tsv';out.to_csv(path,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} genes to {path}')
if __name__=='__main__':main()
