# Parse coverage file
# coverage file - using the mut_infos I should be able to annotate the coverage
# head(mut_infos)

# Output in the coverage file from GATK:
# one row per variant that was given to GATK
# Two columns for every sample (bamfile): 1) Depth_for_Sample 2) Sample_base_counts
# I need to

# Function to extract the depth of coverage of the alt allele for the sample of interest from every row
# To extract the exact depth. I need to know for each variant (row) what's the sample and the alternative allele for which I have to extract the depth
extract_depth_allele <- function(row_gatk_coverage){

  sample <- as.character(row_gatk_coverage$SampleName)
  column_name <- paste0(sample,"_base_counts")
  alt <- row_gatk_coverage$alt_initial

  if(is.na(alt)) { # if there is not available alt allele since I wasn't able to detect it in the initial set
    alt_dep <- NA } else {
      sample_base_counts <- row_gatk_coverage[,column_name]
      alleles <- strsplit(as.character(sample_base_counts),split=" ")[[1]]
      alleles_alt <- alleles[grep(alt,alleles)]
      alt_dep <- gsub(paste0(alt,":"),"",alleles_alt)
    }

  return(alt_dep)

}

# Function to extract the total depth of coverage for the sample of interest from every row
extract_tot_depth <- function(row_gatk_coverage){

  sample <- as.character(row_gatk_coverage$SampleName)
  column_name <- paste0("Depth_for_",sample)
  alleles_tot_depth <- row_gatk_coverage[,column_name]

  return(alleles_tot_depth)

}


# The GATK coverage file comes in a very specific annoying format. I know from the variants called in the initial deep sequenced RNA-Seq samples
# what are the Alt and Ref allele for each one of the variants for which I computed depth and VAF. The alt and ref allele were not reported in the truth set

parse_gatk_coverage <- function(path_to_gatk_coverage,
                           truth_set){

  # Read in coverage data
  cov <- read.delim(path_to_gatk_coverage)

  # Create key column to combine cov with mutation information and get the alt and ref alleles from the initial set of variants
  # since the allele were not reported in the paper (= truth set)

  cov <- merge(cov,truth_set,all.x=TRUE)
  #cov$Location <- gsub(":","_",cov$Locus)
  cov_initial_snv <- subset(cov,VARIANT_CLASS %in% "SNV")

  # Extract tot depth and alt depth
  alt_allele_counts <- NULL
  tot_depth_counts <- NULL
  for(i in 1:nrow(cov_initial_snv)){
    alt_allele_counts <- c(alt_allele_counts,extract_depth_allele(cov_initial_snv[i,]))
    tot_depth_counts <- c(tot_depth_counts,extract_tot_depth(cov_initial_snv[i,]))
  }

  cov_initial_snv$alt_depth_GATK <- as.numeric(as.character(alt_allele_counts))
  cov_initial_snv$tot_depth_GATK <- as.numeric(as.character(tot_depth_counts))

  base_counts_depth <- c(colnames(cov_initial_snv)[grep("Depth_for_",colnames(cov_initial_snv))],
                         colnames(cov_initial_snv)[grep("_base_counts",colnames(cov_initial_snv))])

  # Remove a lot of the extra columns
  cov_initial_snv <- cov_initial_snv[,!(colnames(cov_initial_snv)  %in% c(base_counts_depth,"Total_Depth","Average_Depth_sample"))]
  cov_initial_snv$VAF_GATK <- cov_initial_snv$alt_depth_GATK/cov_initial_snv$tot_depth_GATK

  return(cov_initial_snv)

}


