from pathlib import Path
from io import StringIO
import gzip,numpy as np,pandas as pd
from scipy.stats import ttest_ind
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw';OUT=HERE.parent/'result'
def main():
 m=pd.read_csv(HERE/'sample_manifest.tsv',sep='\t');m=m[m.selected.astype(str).str.upper().eq('TRUE')];lines=[];active=False
 with gzip.open(RAW/'GSE137340_series_matrix.txt.gz','rt',encoding='latin1',errors='replace') as f:
  for l in f:
   if l.startswith('!series_matrix_table_begin'):active=True;continue
   if l.startswith('!series_matrix_table_end'):break
   if active:lines.append(l)
 x=pd.read_csv(StringIO(''.join(lines)),sep='\t',quotechar='"').rename(columns={'ID_REF':'gene_id'});ids=m.gsm.tolist();x=x[['gene_id']+ids];x[ids]=x[ids].apply(pd.to_numeric,errors='coerce');x=x.dropna().groupby('gene_id')[ids].mean();ca=m.loc[m.analysis_group.eq('Sepsis_D1'),'gsm'].tolist();co=m.loc[m.analysis_group.eq('Healthy_Control'),'gsm'].tolist();A=x[ca].to_numpy();B=x[co].to_numpy();e=A.mean(1)-B.mean(1);u=np.sqrt(A.var(1,ddof=1)/len(ca)+B.var(1,ddof=1)/len(co));p=ttest_ind(A,B,axis=1,equal_var=False,nan_policy='omit').pvalue;o=np.argsort(np.nan_to_num(p,nan=1));q=np.empty_like(p);q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=pd.DataFrame({'gene_id':x.index,'E':e,'U':u,'p_value':p,'q_value':np.clip(q,0,1)});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].sort_values('gene_id');OUT.mkdir(exist_ok=True);path=OUT/'GSE137340_Sepsis_Day1_vs_HC_E_U.tsv';out.to_csv(path,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} genes to {path}')
if __name__=='__main__':main()
