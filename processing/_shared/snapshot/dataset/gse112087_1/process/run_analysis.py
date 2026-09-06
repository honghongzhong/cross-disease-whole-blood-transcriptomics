from pathlib import Path
import gzip, io, tarfile
import numpy as np, pandas as pd
from scipy.stats import ttest_ind

HERE=Path(__file__).resolve().parent; RAW=HERE.parent/'raw'; OUT=HERE.parent/'result'
def main():
    m=pd.read_csv(HERE/'sample_manifest.tsv',sep='\t');m=m[m.include_locked.astype(str).str.upper().eq('TRUE')].copy(); lane_to_donor={}
    for _,r in m.iterrows():
        for gsm in str(r.lane_gsm_ids).split(';'): lane_to_donor[gsm]=r.sample_key
    donor_counts={};
    with tarfile.open(RAW/'GSE112087_RAW.tar') as t:
        for name in t.getnames():
            gsm=name.split('_',1)[0]
            if gsm not in lane_to_donor: continue
            raw=gzip.decompress(t.extractfile(name).read()); d=pd.read_csv(io.BytesIO(raw),sep='\t',usecols=['EnsemblID','ReadCount'])
            d['ReadCount']=pd.to_numeric(d.ReadCount,errors='coerce').fillna(0); d=d.groupby('EnsemblID').ReadCount.sum()
            donor_counts[lane_to_donor[gsm]]=d if lane_to_donor[gsm] not in donor_counts else donor_counts[lane_to_donor[gsm]].add(d,fill_value=0)
    x=pd.DataFrame(donor_counts).fillna(0); lib=x.sum();x=np.log2(x.div(lib,axis=1)*1e6+1)
    ca=m.loc[m.disease_group.eq('case'),'sample_key'].tolist();co=m.loc[m.disease_group.eq('healthy_control'),'sample_key'].tolist();A=x[ca].to_numpy();B=x[co].to_numpy();e=A.mean(1)-B.mean(1);u=np.sqrt(A.var(1,ddof=1)/len(ca)+B.var(1,ddof=1)/len(co));p=ttest_ind(A,B,axis=1,equal_var=False,nan_policy='omit').pvalue;o=np.argsort(np.nan_to_num(p,nan=1));q=np.empty_like(p);q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=pd.DataFrame({'gene_id':x.index,'E':e,'U':u,'p_value':p,'q_value':np.clip(q,0,1)});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].sort_values('gene_id');OUT.mkdir(exist_ok=True);path=OUT/'GSE112087_SLE_vs_HC_E_U.tsv';out.to_csv(path,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} genes to {path}')
if __name__=='__main__':main()
