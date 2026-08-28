from pathlib import Path
import gzip,numpy as np,pandas as pd
from scipy.stats import ttest_ind
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw';OUT=HERE.parent/'result'
def main():
 m=pd.read_csv(HERE/'sample_manifest.tsv',sep='\t');m=m[m.selected.astype(str).str.upper().eq('TRUE')];
 with gzip.open(RAW/'GSE186507_MSCCR_Blood_counts.txt.gz','rt',encoding='latin1') as f: names=[z.strip('"') for z in f.readline().split()]
 tab=pd.read_csv(RAW/'GSE186507_MSCCR_Blood_counts.txt.gz',sep=r'\s+',header=None,skiprows=1,compression='gzip',engine='python');tab.columns=['gene_id']+names;ids=m.subject.tolist();x=tab[['gene_id']+ids].copy();x[ids]=x[ids].apply(pd.to_numeric,errors='coerce').fillna(0);x=x.set_index('gene_id');x=np.log2(x.div(x.sum(),axis=1)*1e6+1);ca=m.loc[m.group.eq('CD'),'subject'].tolist();co=m.loc[m.group.eq('Healthy_Control'),'subject'].tolist();A=x[ca].to_numpy();B=x[co].to_numpy();e=A.mean(1)-B.mean(1);u=np.sqrt(A.var(1,ddof=1)/len(ca)+B.var(1,ddof=1)/len(co));p=ttest_ind(A,B,axis=1,equal_var=False,nan_policy='omit').pvalue;o=np.argsort(np.nan_to_num(p,nan=1));q=np.empty_like(p);q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=pd.DataFrame({'gene_id':x.index,'E':e,'U':u,'p_value':p,'q_value':np.clip(q,0,1)});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].sort_values('gene_id');OUT.mkdir(exist_ok=True);path=OUT/'GSE186507_CD_vs_HC_E_U.tsv';out.to_csv(path,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} genes to {path}')
if __name__=='__main__':main()
