"""New whole-vector representation experiments. No gene masking or old runs changed."""
import os
os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'
from pathlib import Path
import json,hashlib,gzip,time
import numpy as np,pandas as pd,torch
from scipy.stats import hypergeom
from scipy.optimize import linear_sum_assignment

ROOT=Path(os.environ['BLOOD_DATA_ROOT']); OUT=Path(os.environ['BLOOD_OUTPUT']); A=OUT/'analysis'; INPUTS=Path(__file__).resolve().parents[1]/'inputs'
EFFECTS=['hedges_g','log2fc','shrunken_log2fc']; RNG=np.random.default_rng(20260905)
torch.set_num_threads(4)
DEVICE='cuda' if torch.cuda.is_available() else 'cpu'
DT=torch.float64
def save(name,rows): pd.DataFrame(rows).to_csv(A/name,sep='\t',index=False)
def bh(p):
 p=np.asarray(p); order=np.argsort(p); q=np.empty(len(p)); q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
def fit(x,k=15):
 # Exact uncentred SVD through a small contrast Gram matrix.
 values,vectors=torch.linalg.eigh(x.T@x)
 values=values.flip(0).clamp_min(0); vectors=vectors.flip(1)
 q=x@vectors[:,:k]/torch.sqrt(values[:k]).clamp_min(1e-12)
 signs=torch.sign(q[q.abs().argmax(dim=0),torch.arange(k,device=DEVICE)]); q=q*signs
 return q,values
def r2(y,p):return 1-((y-p)**2).sum()/((y-y.mean())**2).sum()
def evaluate(q,ys,us,k):
 result=[]
 for y,u in zip(ys,us):
  keep=torch.isfinite(y); l=q[keep,:k]; v=y[keep]; uv=u[keep]
  score=torch.linalg.lstsq(l,v[:,None]).solution[:,0] if DEVICE=='cpu' else torch.linalg.solve(l.T@l,l.T@v)
  p=l@score; w=1/(uv**2+1e-8)
  result.append([r2(v,p).item(),(1-(w*(v-p)**2).sum()/(w*(v-(w*v).sum()/w.sum())**2).sum()).item(),torch.sqrt(((v-p)**2).mean()).item()])
 return np.asarray(result)
def null_batch(x,n,k=15):
 # Shuffle E gene identity independently within each complete training column.
 # U would travel with E; U is not used by this unweighted reference.
 b=torch.stack([torch.stack([x[torch.randperm(x.shape[0],device=DEVICE),c] for c in range(x.shape[1])],dim=1) for _ in range(n)])
 val,vec=torch.linalg.eigh(b.transpose(1,2)@b);val=val.flip(1).clamp_min(0);vec=vec.flip(2)
 return b@vec[:,:,:k]/torch.sqrt(val[:,:k]).clamp_min(1e-12)[:,None,:],val
def null_evaluate(q,ys,k):
 result=[]
 for y in ys:
  keep=torch.isfinite(y);l=q[:,keep,:k];v=y[keep]
  s=torch.linalg.solve(l.transpose(1,2)@l,(l.transpose(1,2)@v[:,None]).expand(q.shape[0],-1,-1))
  p=(l@s)[:,:,0]; result.append((1-((v[None,:]-p)**2).sum(1)/((v-v.mean())**2).sum()).cpu().numpy())
 return np.array(result).T
def main():
 if (A/'COMPLETION.json').exists():raise RuntimeError('Completed run exists')
 A.mkdir(exist_ok=True);torch.manual_seed(20260905)
 protocol={'scope':'additional retrospective experiments; whole-vector projection, not unseen-gene prediction','reference':'uncentred exact SVD on common training gene panel','rank_primary':5,'ranks_descriptive':[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],'null_replicates':199,'domain_null_replicates':99,'bootstrap_independence_blocks':100,'gene_enrichment':'Hallmark ORA, positive/negative top 200 loading genes; BH over all primary tests; 100/300 sensitivity','seed':20260905,'device':DEVICE,'dtype':'float64','old_models_modified':False,'panel_limitation':'common panel defined from original 36 training contrasts; domain tests conditional on this fixed feature panel','created_unix':time.time()}
 (A/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf8')
 tr=pd.read_csv(INPUTS/'training_manifest.tsv',sep='\t')
 te=pd.read_csv(INPUTS/'test_manifest.tsv',sep='\t')
 tr.to_csv(A/'training_manifest.tsv',sep='\t',index=False);te.to_csv(A/'test_manifest.tsv',sep='\t',index=False)
 tables={};source=[]
 for eff in EFFECTS:
  frames=[]
  for cid in list(tr.contrast_id)+list(te.contrast_id):
   p=ROOT/f'contrasts/{cid}/{cid}__{eff}__E_U.tsv.gz'
   frame=pd.read_csv(p,sep='\t',dtype={'gene_id':str}).set_index('gene_id');frames.append(frame)

  tables[eff]=frames
 genes=sorted(set.intersection(*[set(f.index) for f in tables['hedges_g'][:36]]))
 assert len(genes)==6144
 pd.Series(genes,name='gene_id').to_csv(A/'common_genes.tsv',sep='\t',index=False)
 symbols=pd.read_csv(INPUTS/'gene_symbols.tsv',sep='\t',dtype=str).set_index('gene_id').symbol.to_dict()
 sets={}
 resource=OUT/'resources/h.all.v2025.1.Hs.entrez.gmt'
 for line in resource.read_text().splitlines():
  parts=line.split('\t');sets[parts[0]]=set(parts[2:])&set(genes)
 assert len(sets)==50

 metrics=[];spectra=[];nulls=[];boots=[];enrich=[];domains=[];domain_null=[];sensitivity=[];pairs=[];module_scores=[];allQ={}
 blocks=tr.independence_block.to_numpy();unique=np.unique(blocks)
 for eff in EFFECTS:
  frames=tables[eff]; en=np.column_stack([f.effect.reindex(genes) for f in frames]);un=np.column_stack([f.uncertainty.reindex(genes) for f in frames])
  x=torch.tensor(en[:,:36],dtype=DT,device=DEVICE);ys=[torch.tensor(en[:,36+c],dtype=DT,device=DEVICE) for c in range(12)];us=[torch.tensor(un[:,36+c],dtype=DT,device=DEVICE) for c in range(12)]
  q,vals=fit(x);allQ[eff]=q[:,:5].cpu().numpy()
  # Numerical check against independent numpy float64 eigensolver.
  ev=np.linalg.eigvalsh(en[:,:36].T@en[:,:36])[::-1];assert np.allclose(vals.cpu(),ev,rtol=1e-9,atol=1e-8)
  np.savez(A/f'{eff}_reference.npz',genes=np.array(genes),L=q.cpu().numpy(),scores=(q.T@x).cpu().numpy(),eigenvalues=vals.cpu().numpy())
  scores=(q[:,:5].T@x).cpu().numpy().T
  save(f'{eff}_training_scores.tsv',[{'contrast_id':cid,**{f'factor{j+1}':scores[c,j] for j in range(5)}} for c,cid in enumerate(tr.contrast_id)])
  for k in range(1,16):
   spectra.append({'effect_type':eff,'rank':k,'energy_fraction':(vals[:k].sum()/vals.sum()).item()})
   for c,row in enumerate(evaluate(q,ys,us,k)): metrics.append({'effect_type':eff,'rank':k,'contrast_id':te.contrast_id[c],'category_id':te.category_id[c],'R2':row[0],'weighted_R2':row[1],'RMSE':row[2],'genes':int(torch.isfinite(ys[c]).sum())})
  for start in range(0,199,8):
   nq,nv=null_batch(x,min(8,199-start))
   scoresnull=null_evaluate(nq,ys,5)
   for j in range(len(nq)):
    nulls.append({'effect_type':eff,'replicate':start+j,'train_energy_K5':(nv[j,:5].sum()/nv[j].sum()).item(),'median_test_R2':float(np.median(scoresnull[j]))})
  print(eff,'null finished',flush=True)
  for rep in range(100):
   chosen=RNG.choice(unique,len(unique),replace=True);col=np.concatenate([np.flatnonzero(blocks==b) for b in chosen]);qb,_=fit(x[:,col],5)
   overlap=torch.linalg.svdvals(q[:,:5].T@qb).square().mean().item();sim=np.abs(q[:,:5].cpu().numpy().T@qb.cpu().numpy());ri,ci=linear_sum_assignment(-sim)
   boots.append({'effect_type':eff,'replicate':rep,'subspace_overlap':overlap,'matched_factor_cosine_mean':float(sim[ri,ci].mean())})
  # Fixed-panel whole-category exclusion also removes related training blocks.
  for c in range(12):
   cat=te.category_id[c];excluded=set(tr.loc[tr.category_id==cat,'independence_block']);keep=~tr.independence_block.isin(excluded).to_numpy();xd=x[:,keep];qd,_=fit(xd)
   for k in [3,5,10,15]:
    m=evaluate(qd,[ys[c]],[us[c]],k)[0];domains.append({'effect_type':eff,'category_id':cat,'contrast_id':te.contrast_id[c],'rank':k,'training_contrasts':int(keep.sum()),'R2':m[0],'weighted_R2':m[1],'RMSE':m[2]})
   for start in range(0,99,9):
    nq,_=null_batch(xd,min(9,99-start),5);ns=null_evaluate(nq,[ys[c]],5)[:,0]
    for j,v in enumerate(ns):domain_null.append({'effect_type':eff,'category_id':cat,'replicate':start+j,'R2':v})
  print(eff,'domain finished',flush=True)
  # Block equalisation and removing any single training block.
  weights=np.array([1/np.sum(blocks==b) for b in blocks]);qb,_=fit(x*torch.tensor(np.sqrt(weights),device=DEVICE))
  for c,m in enumerate(evaluate(qb,ys,us,5)):sensitivity.append({'effect_type':eff,'analysis':'equal_block_weight','contrast_id':te.contrast_id[c],'R2':m[0]})
  for b in unique:
   qb,_=fit(x[:,blocks!=b],5)
   for c,m in enumerate(evaluate(qb,ys,us,5)):sensitivity.append({'effect_type':eff,'analysis':f'drop_{b}','contrast_id':te.contrast_id[c],'R2':m[0]})
  # Stable reference orientation: maximal absolute loading is positive.
  load=q[:,:5].cpu().numpy();save(f'{eff}_loadings.tsv',[{'gene_id':g,'symbol':symbols.get(g,''),**{f'factor{j+1}':load[i,j] for j in range(5)}} for i,g in enumerate(genes)])
  for factor in range(5):
   for direction,sign in [('positive',1),('negative',-1)]:
    order=np.argsort(-sign*load[:,factor])
    for top in [100,200,300]:
     chosen=set(np.array(genes)[order[:top]])
     for name,gs in sets.items():
      hit=chosen&gs;p=hypergeom.sf(len(hit)-1,len(genes),len(gs),top)
      enrich.append({'effect_type':eff,'factor':factor+1,'direction':direction,'top_n':top,'pathway':name,'set_genes':len(gs),'hits':len(hit),'p':p,'genes':';'.join(sorted(hit)),'symbols':';'.join(str(symbols.get(g,g)) for g in sorted(hit))})
  # Rotation-invariant biological readout: pathway average observed/reconstructed effects.
  for c in range(12):
   y=ys[c];keep=torch.isfinite(y);l=q[keep,:5];v=y[keep];sc=torch.linalg.solve(l.T@l,l.T@v);pred=np.full(len(genes),np.nan);pred[keep.cpu().numpy()]=(l@sc).cpu().numpy();real=y.cpu().numpy()
   for name,gs in sets.items():
    ix=np.array([g in gs for g in genes])&np.isfinite(real)
    if ix.sum()>=10:module_scores.append({'effect_type':eff,'contrast_id':te.contrast_id[c],'pathway':name,'n_genes':int(ix.sum()),'mean_observed_effect':float(real[ix].mean()),'mean_reconstructed_effect':float(pred[ix].mean()),'mean_residual':float((real[ix]-pred[ix]).mean())})
  save('whole_vector_metrics.tsv',metrics);save('spectrum.tsv',spectra);save('null_distribution.tsv',nulls);save('bootstrap.tsv',boots);save('domain_transfer.tsv',domains);save('domain_null.tsv',domain_null);save('block_sensitivity.tsv',sensitivity);save('pathway_projection.tsv',module_scores)
  print(eff,'done',flush=True)
 for i,a in enumerate(EFFECTS):
  for b in EFFECTS[i+1:]:pairs.append({'effect_a':a,'effect_b':b,'subspace_overlap_K5':float(np.mean(np.linalg.svd(allQ[a].T@allQ[b],compute_uv=False)**2))})
 save('cross_effect_subspaces.tsv',pairs)
 ed=pd.DataFrame(enrich);ed['q']=np.nan
 for top in [100,200,300]:m=ed.top_n==top;ed.loc[m,'q']=bh(ed.loc[m,'p'])
 ed.to_csv(A/'hallmark_enrichment.tsv',sep='\t',index=False)
 (A/'COMPLETION.json').write_text(json.dumps({'state':'COMPLETE','effects':3,'null_replicates_per_effect':199,'domain_null_per_effect_domain':99,'bootstrap_per_effect':100,'reference_models_and_domain_refits':'complete','gpu':torch.cuda.get_device_name(0) if DEVICE=='cuda' else None,'dtype':'float64','new_biological_data':False},indent=2),encoding='utf8')
if __name__=='__main__':main()
