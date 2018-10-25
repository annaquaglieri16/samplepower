

parse_gatk_coverage <- function(path_to_gatk_coverage,
                                truth_set,TCGA=FALSE){

# truthset should contain  
  
cov <- read.delim(path_to_gatk_coverage)
cov$Locus <- as.character(cov$Locus)

# Reshape from wide to long the total depth column for every sample
cov_long_tot_depth <- cov %>% tidyr::gather(base_tot_counts,tot_depth, starts_with("Depth_for_")) %>%
  dplyr::select(Locus,base_tot_counts,tot_depth)

# Reshape from wide to long the per allel base depth for every sample
cov_long_alt_depth <- cov %>% gather(base_alt_counts,alleles, contains("base_counts")) %>%
  select(Locus,base_alt_counts,alleles)

# Unique names of the sample to extract information from column
if(TCGA){
  names_samples <- unique(gsub("-",".",truth_set$SampleNameShort))
}else{
  names_samples <- unique(gsub("-",".",truth_set$SampleName))
}

create_match_names_tot <- data.frame(base_tot_counts = unique(cov_long_tot_depth$base_tot_counts))

for(i in 1:nrow(create_match_names_tot)){
  mat <- stringr::str_extract(create_match_names_tot$base_tot_counts[i], names_samples) # extract sample name from every column name
  mat <- mat[!is.na(mat)]
  if(length(mat) < 1) mat <- NA
  create_match_names_tot$SampleName[i] <- mat
}

create_match_names_alt <- data.frame(base_alt_counts = unique(cov_long_alt_depth$base_alt_counts))
for(i in 1:nrow(create_match_names_alt)){
  mat <- stringr::str_extract(create_match_names_alt$base_alt_counts[i], names_samples)
  mat <- mat[!is.na(mat)]
  if(length(mat) < 1) mat <- NA
  create_match_names_alt$SampleName[i] <- mat
}

# Match Samplename to alt and tot depth 
cov_long_alt_depth1 <- cov_long_alt_depth %>% dplyr::left_join(create_match_names_alt) %>% select(- base_alt_counts)
cov_long_tot_depth1 <- cov_long_tot_depth %>% dplyr::left_join(create_match_names_tot) %>% select(- base_tot_counts)

cov_long <- cov_long_alt_depth1 %>% dplyr::full_join(cov_long_tot_depth1)

# Trick for TCGA names
if(TCGA){
  cov_long$SampleName <- gsub("[.]","-",cov_long$SampleName)
  cov_long <- cov_long %>% dplyr::rename(SampleNameShort = SampleName)
}

cov_long_infos <- truth_set %>% left_join(cov_long)

cov_long_infos1 <- cov_long_infos %>% separate(alleles,into = c("A","C","G","T","N"),sep=" ",remove=FALSE)
#head(cov_long_infos1[is.na(cov_long_infos1$A),c("SampleNameShort","Locus","A","C")])
#cov_long_infos1[,c("SampleNameShort","Locus","A","C")]

cov_long_infos1 <- cov_long_infos1 %>% dplyr::rename(tot_depth_GATK = tot_depth)
cov_long_infos1$alt_depth_GATK <- NA
for(i in 1:nrow(cov_long_infos1)){
  allele_alt <- cov_long_infos1$alt_initial[i]
  allele_alt_dep <- cov_long_infos1[i,allele_alt]
  cov_long_infos1$alt_depth_GATK[i] <- as.numeric(as.character(gsub(paste0(allele,":"),"",allele_alt_dep)))
}

return(cov_long_infos1)

}


# Test with TCGA

# # TCGA
# path_to_gatk_coverage <- "/stornext/Genomics/data/AML_RNA/quaglieri.a/GEO_Leucegene_data/scripts/06-Streamline_downsampling/05-DepthOfCoverage/TCGA_initial/TCGA_initial"
# truth_set <- read_csv(file.path("/stornext/Genomics/data/AML_RNA/quaglieri.a/GEO_Leucegene_data/scripts/08-TCGA/02-prepare_for_sens_analysis_data/mut_infos_standardised.csv"))
# truth_set$LongSampleName <- as.character(names40M$LongSampleName)[match(truth_set$SampleName,names40M$SampleName)]
# length(unique(truth_set$LongSampleName)) # 139
# length(unique(truth_set$SampleName)) # 139
# truth_set <- truth_set %>% dplyr::rename(SampleNameShort = SampleName,
#                                          SampleName = LongSampleName,
#                                          Feature=ensembl_transcript_id) %>%
#   select(SampleName,SampleNameShort,Locus,chrom,pos,end,ref,alt_initial,NCBI,SYMBOL,variant_type,Feature,PresentInRNA,
#          NormalRefReads,NormalVarReads,NormalVAF,TumorRefReads,TumorVarReads,TumorVAF,RNARefReads,RNAVarReads,RNAVAF,RNAVAF_myEst)
# 
# cov_gatk <- parse_gatk_coverage(path_to_gatk_coverage = path_to_gatk_coverage,truth_set = subset(truth_set,variant_type %in% "SNV"),TCGA = TRUE)
# 
# # leucegene
# truth_set <- read_csv(file.path("/stornext/Genomics/data/AML_RNA/quaglieri.a/GEO_Leucegene_data/scripts/02-variant_loss_downsampling/02-variant_loss_downsampling_data","mut_for_paper_table.csv"))
# truth_set <- truth_set %>% dplyr::rename(Feature = ensembl_transcript_id,
#                                                 alt_initial = alt_initial_vardict) %>%
#   separate(Locus, into=c("chrom","pos"),remove=FALSE,sep=":")
# 
# path_to_gatk_coverage1 <- "/stornext/Genomics/data/AML_RNA/quaglieri.a/GEO_Leucegene_data/scripts/06-Streamline_downsampling/05-DepthOfCoverage/20M_100/20M_100"
# cov_gatk <- parse_gatk_coverage(path_to_gatk_coverage = path_to_gatk_coverage1,truth_set = subset(truth_set,variant_type %in% "SNV"),TCGA = FALSE)



