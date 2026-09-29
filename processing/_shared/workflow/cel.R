# Normalize selected CEL files before expression matrix construction.
suppressPackageStartupMessages({
  library(affy)
  library(jsonlite)
})
root <- Sys.getenv("BIOLOGY_WORK_ROOT")
if (!nzchar(root)) stop("BIOLOGY_WORK_ROOT is required")
out <- file.path(root, "outputs/stage1_three_effects_20260903_v1/cel_normalized")
for (d in c("gse28750_1", "gse55201_1")) {
  target <- file.path(out, d)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(target, "rma.tsv.gz"))) next
  m <- read.delim(file.path(root, "dataset", d, "process/sample_manifest.tsv"), check.names = FALSE)
  if (d == "gse28750_1") {
    m <- m[toupper(m$include_locked) == "TRUE", ]
    ids <- m$sample_id
  } else {
    ids <- m$sample_key
  }
  files <- list.files(file.path(root, "dataset", d, "raw"), pattern = "[.]CEL([.]gz)?$", full.names = TRUE, ignore.case = TRUE)
  gsm <- sub("_.*", "", basename(files))
  files <- files[match(ids, gsm)]
  stopifnot(!anyNA(files))
  a <- ReadAffy(filenames = files)
  z <- exprs(rma(a))
  colnames(z) <- ids
  write.table(data.frame(probe_id = rownames(z), z, check.names = FALSE), gzfile(file.path(target, "rma.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(target, "sessionInfo.txt"))
  cat("RMA_COMPLETE", d, nrow(z), ncol(z), "\n")
}
