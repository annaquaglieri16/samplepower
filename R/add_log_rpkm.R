#################
## Extract gene expression, compute logRPKM and attach it to the matrix
#################
# 1. get vector of symbols of genes for which a variant was called
# 2. read the list from featurecounts
# 3. compute logRPKM
# 4. Subset only those genes
# 5. associate samples_genes and merge it with teh variants

# 2. Extract and rename counts
reassign_names <- function(sample_names,counts){
  grep_sampleName <- grep(sample_names,colnames(counts))
  return(grep_sampleName)
}


add_log_rpkm <- function(variants, # this is the dataframe
                         gene_expression,
                         sample_names){

  # 1. Read in counts obtained with featurecounts and saved as RDS file
  counts <- readRDS(gene_expression)

  counts1 <- sapply(sample_names,function(name){
    reassign_names(name,counts$counts)
  } )

  colnames(counts$counts)[counts1] <- names(counts1)

  # 3. Compute CPM and convert to gene symbols
  # gene.length = counts$annotation$Length
  rpkmCounts <- edgeR::cpm(counts$counts,log = TRUE)
  rpkmCounts <- rpkmCounts[rownames(rpkmCounts) %in% variants$GeneID,]

  # 4. Convert from wide to long format
  rpkmCounts <- data.frame(rpkmCounts)
  rpkmCounts$GeneID <- rownames(rpkmCounts)
  rpkmCounts_long <- rpkmCounts %>% tidyr::gather(SampleName,LogRpkm,sample_names)
  rpkmCounts_long$SYMBOL <- variants$SYMBOL[match(rpkmCounts_long$GeneID,variants$GeneID)]

  # 5. Associate the correct logRPKM to each sample
  variant_expr <- merge(variants,rpkmCounts_long,all.x=TRUE)
  return(variant_expr)

}

