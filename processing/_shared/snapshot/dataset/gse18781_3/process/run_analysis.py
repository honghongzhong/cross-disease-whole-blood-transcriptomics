from pathlib import Path
import gzip,numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent;OUT=HERE.parent/'result'
def read(p): return pd.read_csv(p,sep='\t',compression='infer')
def main():
 e=read(HERE/'input_E.tsv.gz');u=read(HERE/'input_U.tsv.gz');
 candidates=['gene_id','entrez_id','entrez','probe_id','ID_REF']; key=next(c for c in candidates if c in e.columns and c in u.columns)
 excluded=set(candidates)|{'canonical_symbol','symbol'}
 ek=next(c for c in e.columns if c not in excluded); uk=next(c for c in u.columns if c not in excluded)
 left=e[[key,ek]].copy().rename(columns={ek:'E'}); right=u[[key,uk]].copy().rename(columns={uk:'U'}); left[key]=left[key].astype(str); right[key]=right[key].astype(str)
 out=left.merge(right,on=key,how='inner'); out=pd.DataFrame({'gene_id':out[key], 'E':pd.to_numeric(out['E'],errors='coerce'),'U':pd.to_numeric(out['U'],errors='coerce')});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].drop_duplicates('gene_id').sort_values('gene_id');p=OUT/'processed_E_U.tsv';out.to_csv(p,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} rows to {p}')
if __name__=='__main__':main()
