"""Training-only parallel analysis; assesses rather than assumes dimension five."""
from run_additional_experiments import *
def main():
 torch.manual_seed(20260907)
 p={'purpose':'training-only dimension assessment requested after exploratory K5 analyses','criterion':'retain leading consecutive eigenvalues exceeding rank-matched 95th percentile of 499 gene-permuted training matrices','primary_scaling':'native effect scale, same as target matrix','sensitivity':['equal column L2 norm','one maximum-effective-n contrast per independence block'],'rank_search':'1 to 15 reported; all eigenvalues used for null','test_scores_used':False,'seed':20260907,'warning':'retrospectively specified analysis, not preregistered'}
 (A/'dimension_protocol.json').write_text(json.dumps(p,indent=2),encoding='utf8')
 genes=pd.read_csv(A/'common_genes.tsv',sep='\t',dtype=str).gene_id.tolist();tr=pd.read_csv(A/'training_manifest.tsv',sep='\t')
 tr['neff']=tr.n_case*tr.n_control/(tr.n_case+tr.n_control);representatives=tr.sort_values(['neff','contrast_id'],ascending=[False,True]).drop_duplicates('independence_block').index.to_numpy()
 rows=[];summary=[]
 for eff in EFFECTS:
  cols=[pd.read_csv(ROOT/f'contrasts/{cid}/{cid}__{eff}__E_U.tsv.gz',sep='\t',dtype={'gene_id':str}).set_index('gene_id').effect.reindex(genes) for cid in tr.contrast_id]
  x=torch.tensor(np.column_stack(cols),device=DEVICE,dtype=DT)
  for label,xx in [('native',x),('unit_column_norm',x/torch.linalg.norm(x,dim=0)),('one_per_block',x[:,representatives])]:
   _,v=fit(xx);samples=[]
   for start in range(0,499,16):
    _,vv=null_batch(xx,min(16,499-start));samples.extend(vv.cpu().numpy())
   samples=np.array(samples);cut=np.quantile(samples,.95,axis=0);val=v.cpu().numpy();retain=0
   for k in range(15):
    if k==retain and val[k]>cut[k]:retain+=1
    rows.append({'effect_type':eff,'scaling':label,'rank':k+1,'eigenvalue':val[k],'null_median':np.median(samples[:,k]),'null95':cut[k],'exceeds95':bool(val[k]>cut[k]),'excess_ratio':val[k]/cut[k]})
   summary.append({'effect_type':eff,'scaling':label,'selected_K':retain,'energy_at_selected_K':float(val[:retain].sum()/val.sum())})
   print(summary[-1],flush=True)
 save('dimension_parallel_analysis.tsv',rows);save('dimension_decision.tsv',summary)
 (A/'DIMENSION_COMPLETE.json').write_text(json.dumps({'state':'COMPLETE','permutations_per_configuration':499,'configurations':9}),encoding='utf8')
if __name__=='__main__':main()
