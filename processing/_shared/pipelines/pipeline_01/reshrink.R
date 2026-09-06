.libPaths(c('E:/Biology/tools/stage1_three_effects_20260903/library','E:/Biology/tools/r461_bioc323_gse153315/library',.libPaths()))
suppressPackageStartupMessages({library(ashr);library(jsonlite);library(data.table)})
root='E:/Biology/outputs/stage1_three_effects_20260903_v1/contrasts'
args=commandArgs(TRUE)
for(p in list.files(root,pattern='^qc.json$',recursive=TRUE,full.names=TRUE)){
 q=fromJSON(p);d=dirname(p);if(q$shrinkage=='ashr_continuous_normal_mixture_normal_likelihood')next
 if(length(args)>0 && sum(utf8ToInt(basename(d))) %% 3 != as.integer(args[1]))next
 t=fread(file.path(d,'three_effects.tsv.gz'));t$gene_id=as.character(t$gene_id)
 ash=ashr::ash(t$log2fc,t$se_log2fc,df=NULL,mixcompdist='normal',method='shrink',outputlevel=2)
 t$shrunken_log2fc=get_pm(ash);t$posterior_sd=get_psd(ash);stopifnot(all(is.finite(t$posterior_sd)),all(t$posterior_sd>0))
 file.copy(p,file.path(d,'qc_before_continuous_prior.json'),overwrite=TRUE)
 fwrite(t,file.path(d,'three_effects.tsv.gz'),sep='\t')
 z=data.frame(gene_id=t$gene_id,effect=t$shrunken_log2fc,uncertainty=t$posterior_sd,uncertainty_type='posterior_SD',contrast_id=basename(d))
 fwrite(z,file.path(d,'shrunken_log2fc_E_U.tsv.gz'),sep='\t');saveRDS(ash,file.path(d,'ash_fit.rds'))
 q$shrinkage='ashr_continuous_normal_mixture_normal_likelihood';q$posterior_sd_zero=0;q$min_psd=min(t$posterior_sd);q$shrinkage_pointmass=FALSE;q$shrinkage_weight_prior='uniform'
 write_json(q,p,pretty=TRUE,auto_unbox=TRUE);cat('CONTINUOUS_PRIOR_COMPLETE',basename(d),'\n');flush.console()
}

