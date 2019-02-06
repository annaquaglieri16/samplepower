#' Rename columns of `featureCounts` counts output
#' @param sampleNames character with sample name
#' @param counts matrix of counts from `featureCounts`


# 2. Extract and rename counts
reassign_names <- function(sample_names,counts){
  grep_sampleName <- grep(sample_names,colnames(counts))
  return(grep_sampleName)
}
