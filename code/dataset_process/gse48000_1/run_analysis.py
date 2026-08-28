from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind,rankdata
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw';OUT=HERE.parent/'result'
def qn(a):
 o=np.argsort(a,axis=0);s=np.take_along_axis(a,o,axis=0);v=s.mean(1);z=np.empty_like(a)
 for j in range(a.shape[1]):z[o[:,j],j]=v
 return z
def main():
 m=pd.read_csv(HERE/'sample_manifest.tsv',sep='\t');x=pd.read_csv(RAW/'GSE48000_non_normalized_set1.txt.gz',sep='\t',compression='gzip').set_index('ID_REF');ids=m.subject_id.tolist();x=x[ids].apply(pd.to_numeric,errors='coerce');x=pd.DataFrame(qn(np.log2(x.clip(lower=1).to_numpy())),index=x.index,columns=ids);a=pd.read_csv(RAW/'GPL10558.annot.gz',sep='\t',skiprows=28,dtype=str,encoding='latin1',low_memory=False).rename(columns={'ID':'probe_id','Gene symbol':'gene_id'})[['probe_id','gene_id']].dropna().drop_duplicates('probe_id');x=x.reset_index(names='probe_id').merge(a,on='probe_id').drop(columns='probe_id').groupby('gene_id')[ids].mean();ca=m.loc[m.group.eq('High-risk VTE'),'subject_id'].tolist();co=m.loc[m.group.eq('Healthy'),'subject_id'].tolist();A=x[ca].to_numpy();B=x[co].to_numpy();e=A.mean(1)-B.mean(1);u=np.sqrt(A.var(1,ddof=1)/len(ca)+B.var(1,ddof=1)/len(co));p=ttest_ind(A,B,axis=1,equal_var=False,nan_policy='omit').pvalue;o=np.argsort(np.nan_to_num(p,nan=1));q=np.empty_like(p);q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=pd.DataFrame({'gene_id':x.index,'E':e,'U':u,'p_value':p,'q_value':np.clip(q,0,1)});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].sort_values('gene_id');OUT.mkdir(exist_ok=True);path=OUT/'GSE48000_high_risk_VTE_vs_HC_E_U.tsv';out.to_csv(path,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} genes to {path}')
if __name__=='__main__':main()
