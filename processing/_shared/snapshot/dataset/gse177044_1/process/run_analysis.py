from pathlib import Path
import gzip,numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent;OUT=HERE.parent/'result'
def read(p): return pd.read_csv(p,sep='\t',compression='infer')
def main():
 e=read(HERE/'input_E.tsv.gz');u=read(HERE/'input_U.tsv.gz');
 ek=next(c for c in e.columns if c not in {'gene_id','canonical_symbol','probe_id','entrez_id'});uk=next(c for c in u.columns if c not in {'gene_id','canonical_symbol','probe_id','entrez_id'});key='gene_id' if 'gene_id' in e.columns else ('entrez_id' if 'entrez_id' in e.columns else 'probe_id');
 out=pd.DataFrame({'gene_id':e[key].astype(str),'E':pd.to_numeric(e[ek],errors='coerce'),'U':pd.to_numeric(u[uk],errors='coerce')});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].drop_duplicates('gene_id').sort_values('gene_id');p=OUT/'processed_E_U.tsv';out.to_csv(p,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} rows to {p}')
if __name__=='__main__':main()
