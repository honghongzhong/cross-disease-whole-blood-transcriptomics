"""Run the primary analyses at training-supported dimensions; fold-specific K for domains."""
from run_additional_experiments import *
def main():
 torch.manual_seed(20260908);rng=np.random.default_rng(20260908)
 decision=pd.read_csv(A/'dimension_decision.tsv',sep='\t');ks=decision[decision.scaling=='native'].set_index('effect_type').selected_K.to_dict()
 (A/'selected_dimension_protocol.json').write_text(json.dumps({'K':ks,'domain_K':'rank-matched 95th percentile using 99 permutations of retained training only','primary_null':199,'bootstraps':100,'main_effect':'hedges_g, dimensionless standardized disease effect','retrospective':True,'no_gene_masking':True},indent=2),encoding='utf8')
 tr=pd.read_csv(A/'training_manifest.tsv',sep='\t');te=pd.read_csv(A/'test_manifest.tsv',sep='\t');genes=pd.read_csv(A/'common_genes.tsv',sep='\t',dtype=str).gene_id.tolist();blocks=tr.independence_block.to_numpy();ub=np.unique(blocks)
 sets={p[0]:set(p[2:])&set(genes) for p in [s.split('\t') for s in (OUT/'resources/h.all.v2025.1.Hs.entrez.gmt').read_text().splitlines()]};heme=np.array([g in sets['HALLMARK_HEME_METABOLISM'] for g in genes]);keeptr=~tr.accession.isin(['GSE244401','GSE262659','GSE26049']).to_numpy()
 metrics=[];nulls=[];domains=[];boots=[];sens=[];module=[];spaces={}
 for eff,k in ks.items():
  frames=[pd.read_csv(ROOT/f'contrasts/{cid}/{cid}__{eff}__E_U.tsv.gz',sep='\t',dtype={'gene_id':str}).set_index('gene_id') for cid in list(tr.contrast_id)+list(te.contrast_id)]
  en=np.column_stack([f.effect.reindex(genes) for f in frames]);un=np.column_stack([f.uncertainty.reindex(genes) for f in frames]);x=torch.tensor(en[:,:36],device=DEVICE,dtype=DT);ys=[torch.tensor(en[:,36+c],device=DEVICE,dtype=DT) for c in range(12)];us=[torch.tensor(un[:,36+c],device=DEVICE,dtype=DT) for c in range(12)];q,val=fit(x,k);spaces[eff]=q.cpu().numpy()
  for c,m in enumerate(evaluate(q,ys,us,k)):metrics.append({'effect_type':eff,'K':k,'contrast_id':te.contrast_id[c],'category_id':te.category_id[c],'R2':m[0],'weighted_R2':m[1],'RMSE':m[2]})
  for start in range(0,199,8):
   nq,_=null_batch(x,min(8,199-start),k);ns=null_evaluate(nq,ys,k)
   for j in range(len(nq)):nulls.append({'effect_type':eff,'K':k,'replicate':start+j,'median_R2':float(np.median(ns[j]))})
  for c in range(12):
   cat=te.category_id[c];exclude=set(tr.loc[tr.category_id==cat,'independence_block']);keep=~tr.independence_block.isin(exclude).to_numpy();xx=x[:,keep];qd,vd=fit(xx);nqs=[];nvs=[]
   for start in range(0,99,9):
    nq,nv=null_batch(xx,min(9,99-start));nqs.append(nq);nvs.append(nv)
   nq=torch.cat(nqs);nv=torch.cat(nvs);cut=torch.quantile(nv,.95,dim=0);kd=0
   for j in range(15):
    if j==kd and vd[j]>cut[j]:kd+=1
   if kd:
    m=evaluate(qd,[ys[c]],[us[c]],kd)[0];ns=null_evaluate(nq,[ys[c]],kd)[:,0];p=(1+int((ns>=m[0]).sum()))/100
   else:
    valid=torch.isfinite(ys[c]);v=ys[c][valid];m=[r2(v,torch.zeros_like(v)).item(),float('nan'),torch.sqrt((v*v).mean()).item()];ns=np.array([m[0]]*99);p=1
   domains.append({'effect_type':eff,'category_id':cat,'contrast_id':te.contrast_id[c],'selected_K':kd,'training_contrasts':int(keep.sum()),'R2':m[0],'weighted_R2':m[1],'RMSE':m[2],'null95':float(np.quantile(ns,.95)),'p_permutation':p})
   del nq,nv,nqs,nvs
  for rep in range(100):
   chosen=rng.choice(ub,len(ub),replace=True);cols=np.concatenate([np.flatnonzero(blocks==b) for b in chosen]);qb,_=fit(x[:,cols],k);ov=torch.linalg.svdvals(q.T@qb).square().mean().item();boots.append({'effect_type':eff,'K':k,'replicate':rep,'overlap':ov})
   high=set(np.array(genes)[np.argsort(-qb.square().sum(1).cpu().numpy())[:500]])
   for name in ['HALLMARK_HEME_METABOLISM','HALLMARK_ALLOGRAFT_REJECTION']:
    gs=sets[name];hits=len(high&gs);module.append({'effect_type':eff,'K':k,'replicate':rep,'pathway':name,'enrichment_fold':hits/(500*len(gs)/len(genes))})
  weights=np.array([1/np.sum(blocks==b) for b in blocks]);qw,_=fit(x*torch.tensor(np.sqrt(weights),device=DEVICE),k)
  qr,_=fit(x[:,keeptr],k);qn,_=fit(x[~heme,:],k)
  for label,qq,yy,uu in [('equal_block_weight',qw,ys,us),('exclude_sickle_MPN',qr,ys,us),('reference_on_nonheme',q[~heme,:],[y[~heme] for y in ys],[u[~heme] for u in us]),('refit_without_heme',qn,[y[~heme] for y in ys],[u[~heme] for u in us])]:
   for c,m in enumerate(evaluate(qq,yy,uu,k)):sens.append({'effect_type':eff,'K':k,'analysis':label,'contrast_id':te.contrast_id[c],'R2':m[0]})
  print(eff,'selected-dimension analysis complete',flush=True)
 save('selected_metrics.tsv',metrics);save('selected_null.tsv',nulls);save('selected_bootstrap.tsv',boots);save('selected_sensitivity.tsv',sens);save('selected_module_stability.tsv',module)
 dd=pd.DataFrame(domains);dd['q_permutation']=bh(dd.p_permutation);dd.to_csv(A/'selected_domain_transfer.tsv',sep='\t',index=False)
 ed=pd.read_csv(A/'hallmark_enrichment.tsv',sep='\t');ed=ed[(ed.top_n==200)&(ed.factor<=ed.effect_type.map(ks))].copy();ed['q_selected']=bh(ed.p);ed.to_csv(A/'selected_enrichment.tsv',sep='\t',index=False)
 rows=[]
 for a,b in [('hedges_g','log2fc'),('hedges_g','shrunken_log2fc'),('log2fc','shrunken_log2fc')]:
  s=np.linalg.svd(spaces[a].T@spaces[b],compute_uv=False);rows.append({'effect_a':a,'effect_b':b,'mean_squared_cosine_smaller_space':float(np.mean(s*s))})
 save('selected_cross_effect.tsv',rows)
 (A/'SELECTED_COMPLETE.json').write_text(json.dumps({'state':'COMPLETE','dimensions':ks,'domain_folds':36,'fold_specific_rank_selection':True}),encoding='utf8')
if __name__=='__main__':main()
