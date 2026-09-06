from pathlib import Path
import pandas as pd,numpy as np
HERE=Path(__file__).resolve().parent;OUT=HERE.parent/'result'
def main():
 d=pd.read_csv(HERE/'input_rma_effects.tsv.gz',sep='\t',compression='gzip');d=d[d.mapping_eligible.astype(str).str.upper().eq('TRUE')].copy();out=pd.DataFrame({'gene_id':d.entrez_id.astype(str),'E':pd.to_numeric(d.effect,errors='coerce'),'U':pd.to_numeric(d.se,errors='coerce'),'p_value':pd.to_numeric(d.p_value,errors='coerce'),'q_value':pd.to_numeric(d.adjusted_p_value,errors='coerce')});out=out[np.isfinite(out.E)&np.isfinite(out.U)&(out.U>0)].drop_duplicates('gene_id').sort_values('gene_id');OUT.mkdir(exist_ok=True);p=OUT/'GSE28750_sepsis_vs_HC_E_U.tsv';out.to_csv(p,sep='\t',index=False,float_format='%.10g');print(f'Wrote {len(out)} genes to {p}')
if __name__=='__main__':main()
