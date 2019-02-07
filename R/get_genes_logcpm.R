#' Extract logCPM and raw counts for selected genes and samples
#' @param geneCountsPath path to `.rds` files with output from `featureCounts`
#' @param sampleNames charatcter vector with sample names to be given as output. Sample names do not have to match exactly the column names in the `featureCounts` output but they have at least to contain the string.
#' @param genes character vector with gene `Symbol` to return as output
#' @param ncbi data frame containing mathing `Symbol` and `GeneID` from NCBI.
#' @param label character. Group label to assign as a column to data frame returned as output.


get_genes_logcpm <- function(geneCountsPath,
                             sample_names,
                            genes,
                            ncbi,
                            label){

  # 1. Read in counts
  counts <- readRDS(geneCountsPath)

  # Only keep counts for samples defined in sample_names
  grep_sample_columns <- grep(paste(sample_names,collapse = "|"),colnames(counts$counts))
  counts$counts <- counts$count[,grep_sample_columns]

  counts1 <- sapply(sample_names,function(name){
    reassign_names(sample_names = name,counts = counts$counts)
  } )
  counts1 <- unlist(counts1) # to avoid coercing to list if some sample names are not found in the column names

  colnames(counts$counts)[counts1] <- names(counts1)
  rownames(counts$counts) <- ncbi$Symbol[match(rownames(counts$counts),ncbi$GeneID)]

  cpmCounts <- cpm(counts$counts,log = TRUE)
  cpmCounts <- cpmCounts[rownames(cpmCounts) %in% genes,]
  cpmCounts <- data.frame(t(cpmCounts))
  cpmCounts$label <- label
  cpmCounts$SampleName <- rownames(cpmCounts)

  rawCounts <- counts$counts[rownames(counts$counts) %in% genes,]
  rawCounts <- data.frame(t(log2(rawCounts+0.5)))
  rawCounts$label <- label
  rawCounts$SampleName <- rownames(rawCounts)

  list(cpmCounts=cpmCounts,rawCounts=rawCounts)

}
