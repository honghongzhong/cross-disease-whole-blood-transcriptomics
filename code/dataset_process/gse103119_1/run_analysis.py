from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import ttest_ind

HERE=Path(__file__).resolve().parent; RAW=HERE.parent/'raw'; OUT=HERE.parent/'result'
def main():
    m=pd.read_csv(HERE/'sample_manifest.tsv',sep='\t'); m=m[m.sample_lock_status.str.startswith('pass')].copy()
    raw=pd.read_csv(RAW/'GSE103119_non-normalized.txt.gz',sep='\t',compression='gzip')
    ids=m.sample_key.tolist(); sig=ids
    if any(x not in raw.columns for x in sig): raise ValueError('sample column missing')
    pos=[raw.columns.get_loc(x) for x in ids]; dcols=[raw.columns[i+1] for i in pos]
    s=raw[sig].apply(pd.to_numeric,errors='coerce').to_numpy(); d=raw[dcols].apply(pd.to_numeric,errors='coerce').to_numpy()
    keep=(np.sum(d<0.01,axis=1)>=5)&np.all(np.isfinite(s)&(s>0),axis=1)
    x=pd.DataFrame(s[keep],columns=ids,index=raw.loc[keep,'ID_REF'].astype(str))
    a=pd.read_csv(RAW/'GPL10558.annot.gz',sep='\t',skiprows=28,dtype=str,encoding='latin1',low_memory=False).rename(columns={'ID':'probe_id','Gene symbol':'gene_id'})[['probe_id','gene_id']].dropna().drop_duplicates('probe_id')
    x=x.reset_index(names='probe_id').merge(a,on='probe_id').drop(columns='probe_id'); x=x[x.gene_id.str.match(r'^[A-Za-z0-9_.-]+$',na=False)].groupby('gene_id')[ids].mean()
    ca=m.loc[m.disease_group.eq('case'),'sample_key'].tolist(); co=m.loc[m.disease_group.eq('healthy_control'),'sample_key'].tolist(); A=x[ca].to_numpy();B=x[co].to_numpy()
    e=A.mean(1)-B.mean(1); u=np.sqrt(A.var(1,ddof=1)/len(ca)+B.var(1,ddof=1)/len(co)); p=ttest_ind(A,B,axis=1,equal_var=False,nan_policy='omit').pvalue
    o=np.argsort(np.nan_to_num(p,nan=1)); q=np.empty_like(p); q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=pd.DataFrame({'gene_id':x.index,'E':e,'U':u,'p_value':p,'q_value':np.clip(q,0,1)});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].sort_values('gene_id');OUT.mkdir(exist_ok=True);path=OUT/'GSE103119_Pneumonia_vs_HC_E_U.tsv';out.to_csv(path,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} genes to {path}')
if __name__=='__main__':main()
