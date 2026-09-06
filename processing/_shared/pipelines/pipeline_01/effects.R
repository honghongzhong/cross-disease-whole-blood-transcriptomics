.libPaths(c('E:/Biology/tools/stage1_three_effects_20260903/library','E:/Biology/tools/r461_bioc323_gse153315/library',.libPaths()))
suppressPackageStartupMessages({library(limma);library(edgeR);library(ashr);library(jsonlite);library(data.table)})
args=commandArgs(TRUE); root='E:/Biology/outputs/stage1_three_effects_20260903_v1'
canonical=as.character(read.delim(file.path(root,'canonical_gene_universe.tsv'))$gene_id)
paths=list.files(file.path(root,'contrasts'),pattern='input.json',recursive=TRUE,full.names=TRUE)
for(p in paths){
 d=dirname(p);id=basename(d);if(length(args)>0 && !id %in% args)next
 if(file.exists(file.path(d,'qc.json')))next
 tryCatch({
  info=fromJSON(p);x=as.matrix(read.delim(gzfile(file.path(d,'input_expression.tsv.gz')),row.names=1,check.names=FALSE));storage.mode(x)='double'
  n1=info$n_case;n0=info$n_control;grp=c(rep(1,n1),rep(0,n0));design=cbind(intercept=1,case=grp)
  before=nrow(x);transform='deposited_log2_expression';nf=NULL
  if(info$input_kind=='counts'){
   y=DGEList(x);keep=filterByExpr(y,design=design);y=y[keep,,keep.lib.sizes=FALSE];y=calcNormFactors(y,method='TMM')
   v=voom(y,design,plot=FALSE);x=v$E;fit=lmFit(v,design);transform='TMM_filterByExpr_voom_log2CPM';nf=y$samples
  } else {
   # Explicit RNA abundance routes and raw array routes. Their source/scale is audited separately.
   src=tolower(paste(info$source,info$directory));rawarray=grepl('non.normalized|non_normalized',src) || info$directory=='gse51404_1'
   abundance=grepl('fpkm|tpm|gse145412|gse174056|gse83563|gse202518|gse226238|gse236442',src)
   if(rawarray && max(x)>100){x=normalizeBetweenArrays(log2(pmax(x,1)),method='quantile');transform='log2_floor1_quantile_raw_array'}
   else if(abundance && max(x)>50){stopifnot(min(x)>=0);x=log2(x+1);transform='log2_normalized_abundance_plus1'}
   else if(quantile(x,.95)>100){x=normalizeBetweenArrays(log2(pmax(x,1)),method='quantile');transform='log2_floor1_quantile_linear_array'}
   stopifnot(max(abs(x))<100)
   fit=lmFit(x,design)
  }
  # Sampling (unmoderated) coefficient SE; shrinkage is performed separately by ashr.
  beta=fit$coefficients[,2];se=fit$stdev.unscaled[,2]*fit$sigma;df=fit$df.residual
  a=x[,seq_len(n1),drop=FALSE];b=x[,n1+seq_len(n0),drop=FALSE]
  v1=apply(a,1,var);v0=apply(b,1,var);dfg=n1+n0-2;sdpool=sqrt(((n1-1)*v1+(n0-1)*v0)/dfg)
  J=exp(lgamma(dfg/2)-.5*log(dfg/2)-lgamma((dfg-1)/2));delta=(rowMeans(a)-rowMeans(b))/sdpool;g=J*delta
  seg=J*sqrt((n1+n0)/(n1*n0)+delta^2/(2*dfg))
  noncurrent=!rownames(x)%in%canonical
  if(any(noncurrent))fwrite(data.frame(gene_id=rownames(x)[noncurrent]),file.path(d,'excluded_noncurrent_entrez.tsv.gz'),sep='\t')
  keep=is.finite(beta)&is.finite(se)&se>0&is.finite(g)&is.finite(seg)&seg>0&!noncurrent
  beta=beta[keep];se=se[keep];df=df[keep];g=g[keep];seg=seg[keep];x=x[keep,,drop=FALSE]
  stopifnot(length(beta)>=1000)
  stopifnot(length(unique(df))==1)
  ash=ashr::ash(beta,se,df=NULL,mixcompdist='normal',method='shrink',outputlevel=2)
  pm=get_pm(ash);ps=get_psd(ash);stopifnot(all(is.finite(pm)),all(is.finite(ps)),all(ps>0))
  tab=data.frame(gene_id=rownames(x),hedges_g=g,se_hedges_g=seg,log2fc=beta,se_log2fc=se,shrunken_log2fc=pm,posterior_sd=ps,df_log2fc=df,n_case=n1,n_control=n0)
  fwrite(tab,file.path(d,'three_effects.tsv.gz'),sep='\t')
  fwrite(data.frame(gene_id=rownames(x),x,check.names=FALSE),file.path(d,'normalized_expression.tsv.gz'),sep='\t')
  for(view in c('hedges_g','log2fc','shrunken_log2fc')){
   u=switch(view,hedges_g='se_hedges_g',log2fc='se_log2fc',shrunken_log2fc='posterior_sd')
   t=data.frame(gene_id=tab$gene_id,effect=tab[[view]],uncertainty=tab[[u]],uncertainty_type=if(view=='shrunken_log2fc')'posterior_SD' else 'sampling_SE',contrast_id=id)
   fwrite(t,file.path(d,paste0(view,'_E_U.tsv.gz')),sep='\t')
  }
  saveRDS(ash,file.path(d,'ash_fit.rds'));if(!is.null(nf))write.table(nf,file.path(d,'normalization_factors.tsv'),sep='\t',quote=FALSE)
  qc=list(status='THREE_EFFECTS_COMPLETE',contrast_id=id,n_case=n1,n_control=n0,input_genes=before,common_valid_genes=nrow(tab),expression_transform=transform,log2fc_model=if(info$input_kind=='counts')'voom_weighted_least_squares' else 'OLS_on_log2_expression',shrinkage='ashr_continuous_normal_mixture_normal_likelihood',posterior_sd_zero=sum(ps==0),min_psd=min(ps),max_abs_input=max(abs(x)),g_variance='J^2*((n1+n0)/(n1*n0)+d^2/(2*df))',scope=info$scope,limitations='Unadjusted disease contrast. Shared participants/controls require grouped later validation. All-gene empirical prior must be re-fit within training folds for blind evaluation.')
  write_json(qc,file.path(d,'qc.json'),pretty=TRUE,auto_unbox=TRUE);cat('COMPLETE',id,nrow(tab),'genes\n');flush.console()
 },error=function(e){cat('ERROR',id,conditionMessage(e),'\n');flush.console();writeLines(conditionMessage(e),file.path(d,'effect_error.txt'))})
}
writeLines(capture.output(sessionInfo()),file.path(root,'R_sessionInfo.txt'))




