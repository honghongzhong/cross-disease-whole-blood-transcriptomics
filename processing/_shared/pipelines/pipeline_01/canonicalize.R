.libPaths(c('E:/Biology/tools/stage1_three_effects_20260903/library','E:/Biology/tools/r461_bioc323_gse153315/library',.libPaths()))
suppressPackageStartupMessages({library(ashr);library(jsonlite);library(data.table)})
root='E:/Biology/outputs/stage1_three_effects_20260903_v1';canonical=as.character(fread(file.path(root,'canonical_gene_universe.tsv'))$gene_id);args=commandArgs(TRUE)
for(p in list.files(file.path(root,'contrasts'),pattern='^qc.json$',recursive=TRUE,full.names=TRUE)){
 d=dirname(p);if(sum(utf8ToInt(basename(d))) %% 3 != as.integer(args[1]))next
 q=fromJSON(p);t=fread(file.path(d,'three_effects.tsv.gz'));t$gene_id=as.character(t$gene_id);bad=!t$gene_id%in%canonical;if(!any(bad))next
 fwrite(t[bad,],file.path(d,'excluded_noncurrent_entrez.tsv.gz'),sep='\t');t=t[!bad,]
 ash=ashr::ash(t$log2fc,t$se_log2fc,df=NULL,mixcompdist='normal',method='shrink',outputlevel=2);t$shrunken_log2fc=get_pm(ash);t$posterior_sd=get_psd(ash);stopifnot(all(t$posterior_sd>0))
 fwrite(t,file.path(d,'three_effects.tsv.gz'),sep='\t');saveRDS(ash,file.path(d,'ash_fit.rds'))
 x=fread(file.path(d,'normalized_expression.tsv.gz'));x$gene_id=as.character(x$gene_id);x=x[x$gene_id%in%canonical,];fwrite(x,file.path(d,'normalized_expression.tsv.gz'),sep='\t')
 for(view in c('hedges_g','log2fc','shrunken_log2fc')){
  u=switch(view,hedges_g='se_hedges_g',log2fc='se_log2fc',shrunken_log2fc='posterior_sd');z=data.frame(gene_id=t$gene_id,effect=t[[view]],uncertainty=t[[u]],uncertainty_type=if(view=='shrunken_log2fc')'posterior_SD' else 'sampling_SE',contrast_id=basename(d));fwrite(z,file.path(d,paste0(view,'_E_U.tsv.gz')),sep='\t')
 }
 q$common_valid_genes=nrow(t);q$noncurrent_entrez_excluded=sum(bad);q$min_psd=min(t$posterior_sd);write_json(q,p,pretty=TRUE,auto_unbox=TRUE);cat('CANONICAL_COMPLETE',basename(d),sum(bad),'excluded\n');flush.console()
}
