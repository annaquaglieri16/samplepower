# Parse coverage file
# coverage file - using the mut_infos I should be able to annotate the coverage
# head(mut_infos)

# Output in the coverage file from GATK:
# one row per variant that was given to GATK
# Two columns for every sample (bamfile): 1) Depth_for_Sample 2) Sample_base_counts
# I need to

# Function to extract the depth of coverage of the alt allele for the sample of interest from every row
# To extract the exact depth. I need to know for each variant (row) what's the sample and the alternative allele for which I have to extract the depth
extract_depth_allele <- function(row_gatk_coverage,TCGA=FALSE){

  if(TCGA){ # quite of a bad fix just for TCGA names
    sample <- gsub("-",".",as.character(row_gatk_coverage$SampleNameShort))
    base_counts_columns <- colnames(row_gatk_coverage)[grep("_base_counts",colnames(row_gatk_coverage))]
    sample_pmatch <- grep(sample,base_counts_columns)
    column_name <- base_counts_columns[sample_pmatch]
  }else{
    sample <- as.character(row_gatk_coverage$SampleName)
    column_name <- paste0(sample,"_base_counts")
  }

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
extract_tot_depth <- function(row_gatk_coverage,TCGA=FALSE){

  if(TCGA){ # quite of a bad fix just for TCGA names
    sample <- gsub("-",".",as.character(row_gatk_coverage$SampleNameShort))
    depth_columns <- colnames(row_gatk_coverage)[grep("Depth_for_",colnames(row_gatk_coverage))]
    sample_pmatch <- grep(sample,depth_columns)
    column_name <- depth_columns[sample_pmatch]
  }else{
    sample <- as.character(row_gatk_coverage$SampleName)
    column_name <- paste0(sample,"_base_counts")
  }

  alleles_tot_depth <- row_gatk_coverage[,column_name]

  return(alleles_tot_depth)

}


# The GATK coverage file comes in a very specific annoying format. I know from the variants called in the initial deep sequenced RNA-Seq samples
# what are the Alt and Ref allele for each one of the variants for which I computed depth and VAF. The alt and ref allele were not reported in the truth set

parse_gatk_coverage <- function(path_to_gatk_coverage,
                           truth_set,TCGA=FALSE){

  # Read in coverage data
  cov <- read.delim(path_to_gatk_coverage)
  cov$Locus <- as.character(cov$Locus)
  # Create key column to combine cov with mutation information and get the alt and ref alleles from the initial set of variants
  # since the allele were not reported in the paper (= truth set)

  cov <- merge(cov,truth_set,all.x=TRUE)
  #cov$Location <- gsub(":","_",cov$Locus)fVARI

  # Extract tot depth and alt depth
  alt_allele_counts <- NULL
  tot_depth_counts <- NULL
  for(i in 1:nrow(cov)){
    alt_allele_counts <- c(alt_allele_counts,extract_depth_allele(cov[i,],TCGA=TCGA))
    tot_depth_counts <- c(tot_depth_counts,extract_tot_depth(cov_initial_snv[i,],TCGA=TCGA))
  }

  cov$alt_depth_GATK <- as.numeric(as.character(alt_allele_counts))
  cov$tot_depth_GATK <- as.numeric(as.character(tot_depth_counts))

  base_counts_depth <- c(colnames(cov)[grep("Depth_for_",colnames(cov))],
                         colnames(cov)[grep("_base_counts",colnames(cov))])

  # Remove a lot of the extra columns
  cov <- cov[,!(colnames(cov)  %in% c(base_counts_depth,"Total_Depth","Average_Depth_sample"))]
  cov$VAF_GATK <- cov$alt_depth_GATK/cov$tot_depth_GATK

  return(cov)

}


