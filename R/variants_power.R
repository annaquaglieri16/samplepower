# The Variants called from the different callers that are used here would have been annotated with VEP and parsed with my functions so that I have control over the columns used for filtering
# What are the columns from the truth set that I use?
# What are the columns from the the variant called that I use?

# Columns from the call set that I use
# (Inside fieldsExtract)
# Existing_variation <- ifelse(Existing_variation == "",NA,Existing_variation)
# ExAC_AF <- as.numeric(as.character(ExAC_AF))
# EUR_AF <- as.numeric(as.character(EUR_AF))
# SYMBOL <- ifelse(SYMBOL == "",NA,SYMBOL)
# SampleName <- gsub("_","",SampleName)
# key_SampleName
# chrom
# pos
# Feature
# ADJVAF_ADJ_indels created by my parsing function

# Inside flag var
# 1. Location (chr_position)
# alt (alternative allele)
# 'nsam' is a column from normal_variants specifying how many nomral samples share that variant
# 'minVAF' is a column from normal_variants specifying the minimum VAF reported in one of the variants found in normals
# 'filter' column from my parsing
# 'qual', 'alt_depth', 'tot_depth', 'VAF' column from my parsing
# the file flag_patterns has specific columns used in the function

# Fields shared between truth and call sets
# 1.gene name always SYMBOL in both truth and called set
# 2. SampleName always SampleName in both truth and called set
# 3. variant_type column for SNV/INDEL in both truth and called set
# 4. chrom always chrom
# 5. pos always pos
# 6. alt and ref columns to parse the coverage file


# Fields from truth set needed to be that column name
# Locus: chr:position (same as ouput from GATK coverage)
# alt_initial: alt allele of initial dataset. Needs to be distinguished form the alt in the call set
# ref : reference allele
# SampleName : SRX name
# Feature : ensemble_transcript_id
# variant_type: SNV/INDEL
# chrom
# pos
# SYMBOL: gene name in SYMBOL

variants_power <- function(vcfdir = vcf_germline_var,
                          down_label,
                          germline = TRUE,
                          truth_set = mut_infos_ra,
                          vcfdir_initial= variant_initial,
                          caller="mutect",
                          gene_expression = gene_expr,
                          normal_variants = mut_pon,
                          exon_ranges=hg19_exons_combined,
                          homop_ranges=GRanges_homop,
                          RNAedit_ranges = RNAeditGR,
                          repeats_ranges = repeatsGR,
                          ncbi = ncbi,
                          flag_patterns = file.path(dir,"scripts/02-variant_loss_downsampling/flag_patterns_added2.csv"),
                          path_to_gatk_coverage = file.path(dir, "scripts/06-Streamline_downsampling/05-DepthOfCoverage/50M_100/50M_100")){

  #################################
  # 1. Read in Downsampled variants
  #################################

  # Germline - read variants in for the run
  names_samples <- gsub("_germline_final.txt","",basename(vcfdir))
  list_variants <- lapply(vcfdir,function(x){
    read.delim(x,fill=TRUE, stringsAsFactors = FALSE)}
    )

  # Check variant number and samples
  len_variants <- sapply(list_variants,nrow)
  len_col <- sapply(list_variants,ncol)

  if(length(len_variants) < length(names_samples)){
    message("WARNING: Some samples' variants might be corrupted. No variants returned.")
  }

  if(length(unique(len_col)) > 1){
    stop("Some samples have different number of columns reported. Check the variant annotation and parsing steps.")
  }

  print("Downsampled variants read in.")

  # Extract interesting fields from annotation and combine variants
  # In the step below variants not on genes are removed to reduce the size of the call set
  down_variants_sub <- lapply(1:length(list_variants),
                              function(index){
                                extract_fields(list_variants[[index]],label=down_label)}
                              )
  down_variants_sub <- do.call(rbind,down_variants_sub)

  #############################
  # 1. Read in Initial variants
  #############################

  names_samples <- gsub("_germline_final.txt","",basename(vcfdir_initial))
  list_variants <- lapply(vcfdir_initial,
                          function(x){
                            read.delim(x,fill=TRUE, stringsAsFactors = FALSE)}
                          )

  # Check variant number and samples
  len_variants <- sapply(list_variants,nrow)
  len_col <- sapply(list_variants,ncol)

  if(length(len_variants) < length(names_samples)){
    message("WARNING: Some samples' variants in the initial variant calls might be corrupted. No variants returned.")
  }

  if(length(unique(len_col)) > 1){
    stop("Some samples in the initial variant calls have different number of columns reported. Check the variant annotation and parsing steps.")
  }

  print("Initial variants read in.")

  # Extract interesting fields from annotation and combine variants
  # In the step below variants not on genes are removed to reduce the size of the call set
  variants_init <- lapply(1:length(list_variants),
                          function(index){
                            extract_fields(list_variants[[index]],label=down_label)}
                          )
  variants_init <- do.call(rbind,variants_init)

  ##############
  # 2. Filtering
  ##############

  # This step will add flagged column and Keep columns to decide whether a variant is passed by a caller
  # Downsampled
  variants_down_filtered <- flag_variants(variants = down_variants_sub,
                                    normal_variants = normal_variants,
                                    flag_patterns = flag_patterns,
                                    exon_ranges = exon_ranges,
                                    homop_ranges = homop_ranges,
                                    RNAedit_ranges = RNAeditGR,
                                    repeats_ranges = repeatsGR)

  # Only select variants falling on genes that we have to study
  variants_down_filtered <- subset(variants_down_filtered, SYMBOL %in% truth_set$SYMBOL)

  # Not really a good strategy but for the moment I will leave it. This step is beacuse later on I want to rbind all the calls
  #if(caller == "mutect" | caller == "varscan"){
  #  variants_down_filtered$ADJVAF_ADJ_indels <- NA
  #} # use bind_rows

  # Initial
  variants_init_filtered <- flag_variants(variants = variants_init,
                                    normal_variants = normal_variants,
                                    flag_patterns = flag_patterns,
                                    exon_ranges = exon_ranges,
                                    homop_ranges = homop_ranges,
                                    RNAedit_ranges = RNAeditGR,
                                    repeats_ranges = repeatsGR)
  variants_init_filtered <- subset(variants_init_filtered, SYMBOL %in% truth_set$SYMBOL)

  #if(caller == "mutect" | caller == "varscan"){
  #  variants_init_filtered$ADJVAF_ADJ_indels <- NA
  #}

  # Truth set: Variants in paper
  # 3. Read in coverage computed for every location - only for the 58 SNVs in the paper
  # Downsampled
  coverage_downsampled <- parse_gatk_coverage(truth_set = truth_set,
                                         path_to_gatk_coverage = as.character(path_to_gatk_coverage))

  # VAF_GATK will be from GATK depth of Cov and VAF the estimate from each caller
  coverage_downsampled <- unique(merge(coverage_downsampled,
                                       variants_down_filtered,all.x=TRUE))

  # Add label for downsampling run and caller
  coverage_downsampled$down_label <- down_label
  coverage_downsampled$caller <- unique(variants_down_filtered$caller)


  ###################################################
  # Define if a SNV was called or not in this run:
  # Using defualt and annot filters as well as the match with the alt allele
  ###################################################

  coverage_downsampled <- coverage_downsampled %>% mutate(Called_annot = ifelse(!is.na(Keep_annot) & Keep_annot,1,0),
                                                          Called_defaults = ifelse(!is.na(Keep_defaults) & Keep_defaults,1,0),
                                                          match_alt = ifelse(alt == alt_initial,1,0))

  # Alt is obtained from GATK both is the variant is called or not called so there will always an ALT
  coverage_downsampled <- coverage_downsampled %>% mutate(Called_annot =  ifelse(Called_annot == 1 & match_alt == 1,1,0),
                                                          Called_defaults = ifelse(Called_defaults == 1 & match_alt == 1,1,0))
  # add GeneID
  coverage_downsampled$GeneID <- ncbi$GeneID[match(coverage_downsampled$SYMBOL,ncbi$SYMBOL)]

  # 4. Add gene expression counts
  sampleNames <- as.character(unique(truth_set$SampleName))
  coverage_downsampled_expr <- add_log_rpkm(variants = coverage_downsampled,
                                            gene_expression = gene_expression,
                                          sample_names = sampleNames)

  # 5. Determine power of recovery
  power <- variant_recovery_power(variants = coverage_downsampled,
                                truth_set = truth_set)

  #############################################
  ## Sens and FP using initial set as truth set
  #############################################


  variants_init_filtered$key1_SampleName <- paste(variants_init_filtered$chrom,
                                                  variants_init_filtered$pos,
                                                  variants_init_filtered$SampleName,
                                                  variants_init_filtered$alt,
                                                  variants_init_filtered$SYMBOL,sep=";")

  variants_down_filtered$key1_SampleName <- paste(variants_down_filtered$chrom,
                                                  variants_down_filtered$pos,
                                                  variants_down_filtered$SampleName,
                                                  variants_down_filtered$alt,
                                                  variants_down_filtered$SYMBOL,sep=";")

  # I use the dbSNP, Exac databases to annotate the variants and set the flags.
  # The same variant can be annotated differently, in some cases they have or they don't have tdbSNP for instance.
  # This means that if I use the flag term I will have duplicated variants
  # What I am interested here is to see how much different library sizes affect the false positive rate at different levels of downsampling.
  variants_init_filtered_unique <- unique(variants_init_filtered[,c("Location","caller","chrom","pos","ref","alt",
                                                                    "key1_SampleName","down_label","Keep_defaults","Keep_annot")])



  variants_down_filtered_unique <- unique(variants_down_filtered[,c("Location","caller","chrom","pos","ref","alt",
                                                                    "key1_SampleName","Keep_defaults","Keep_annot")])

  # Match if a variant present in the initial run was found in the downsampled dataset.
  # Is there a match with the variants called in the downsampled dataset?
  # DownMatch_defaults is 1 if a variant in the initial file is called in the downsampled set

  variants_init_filtered_unique <- variants_init_filtered_unique %>%
    mutate(DownMatch_defaults = ifelse(key1_SampleName %in% variants_down_filtered_unique$key1_SampleName[variants_down_filtered_unique$Keep_defaults],TRUE,FALSE)) %>%
    mutate(DownCalled_defaults = ifelse(DownMatch_defaults & Keep_defaults,1,0))


  # Same for annot strategy
  variants_init_filtered_unique <- variants_init_filtered_unique %>%
    mutate(DownMatch_annot = ifelse(key1_SampleName %in% variants_down_filtered_unique$key1_SampleName[variants_down_filtered_unique$Keep_annot],TRUE,FALSE)) %>%
    mutate(DownCalled_annot = ifelse(DownMatch_annot & Keep_annot,1,0))

  ######################
  # Compute Sensitivity
  #####################

  sens_defaults <- sum(variants_init_filtered_unique$DownCalled_defaults)/sum(variants_init_filtered_unique$Keep_defaults)
  sens_annot <- sum(variants_init_filtered_unique$DownCalled_annot)/sum(variants_init_filtered_unique$Keep_annot)

  sens = c(sens_defaults,sens_annot)
  names(sens) <- c("defaults","annotations")

  #############
  # Compute FP
  #############
  # Sum all the variants called by a caller in the downsampled set that are not in the initial set
  # 1. only take those variants that matched a variant called in the initial set with default parameters
  # Sum the keep_defaults from the downsampled dataset

  variants_passed_in_initial_defaults <- variants_init_filtered_unique$key1_SampleName[variants_init_filtered_unique$Keep_defaults]

  # Number of variantscalled in the downsmapled call set nott present in the initial call set
  fpn_defaults <- sum(variants_down_filtered_unique$Keep_defaults[!(variants_down_filtered_unique$key1_SampleName %in% var_passed_in_initial_defaults)])
  fp_defaults <- fpn_defaults/sum(variants_down_filtered_unique$Keep_defaults)

  var_passed_in_initial_annot <- variants_init_filtered_unique$key1_SampleName[variants_init_filtered_unique$Keep_annot]
  fpn_annot <- sum(variants_down_filtered_unique$Keep_annot[!(variants_down_filtered_unique$key1_SampleName %in% var_passed_in_initial_annot)])
  fp_annot <- fpn_annot/sum(variants_down_filtered_unique$Keep_annot)

  fp = c(fpn_defaults,fp_defaults,fpn_annot,fp_annot)
  names(fp) <- c("Ndefaults","Pdefaults","Nannot","Pannot")

  power_init <- list(sens = sens,
                     fp = fp,
                     variants_init_filtered_unique = variants_init_filtered_unique)

  ######################################################
  # Add flag in the variant set: 1 if FP and 0 otherwise
  ######################################################

  # variants_down_filtered differently from variants_down_filtered_unique containes multiple rows for the same variant if annotated on different transcripts.
  variants_down_filtered$FP_defaults <- ifelse(!(variants_down_filtered$key1_SampleName %in% variants_init_filtered_unique$key1_SampleName) &
                                                 variants_down_filtered$Keep_defaults, 1, 0)
  variants_down_filtered$FP_annot <- ifelse(!(variants_down_filtered$key1_SampleName %in% variants_init_filtered_unique$key1_SampleName) &
                                              variants_down_filtered$Keep_annot, 1, 0)

  #############################################################
  # Add gene expression for all variants called by this caller
  ##############################################################

  # Add gene ID and expression
  variants_down_filtered$GeneID <- ncbi$GeneID[match(variants_down_filtered$SYMBOL,ncbi$SYMBOL)]
  variants_down_filtered_expr <- add_log_rpkm(variants = variants_down_filtered,
                                              gene_expression = gene_expression,
                                              sample_names = sampleNames)

  #####################
  # Recovery of INDELS
  #####################

  # Select everything thet it is not a SNV: deletions, insetions. indels, substitutions
  indels_paper <- subset(truth_set,VARIANT_CLASS %in% "INDEL")
  indels_caller <- subset(variants_down_filtered_expr,!(VARIANT_CLASS %in% "SNV"))

  # Check if the variant passes the annotation filter or the default filters: Keep it only if it passes
  indels_caller$Called_annot <- ifelse(!is.na(indels_caller$Keep_annot) & indels_caller$Keep_annot,1,0)
  indels_caller$Called_defaults <- ifelse(!is.na(indels_caller$Keep_defaults) & indels_caller$Keep_defaults,1,0)

  # See how many variants from the truth set are recovered
  indels_paper_recovery <- match_indels(truth_set = indels_paper,
                                        variants = indels_caller)

  # The sensitivity will be donw with Manual curation to see if the INDELs were called and outside of the function I will compute the sensitivity.

  # Reorder columns
  indels_paper_recovery <- indels_paper_recovery %>%
    dplyr::select(chrom,pos,ref,alt,alt_initial,Mutation,Locus,indel_features,VAF,VAF_paper,ref_depth,alt_depth,SYMBOL,
                  type_indel,type_indel2,caller,down_label,key_SampleName,SampleName,everything())


  #############################################################
  ## Sens and FP using initial set as truth for INDELs for which I don't need manual curation
  #############################################################

  indels_init_filtered <- subset(variants_init_filtered,!(VARIANT_CLASS %in% "SNV"))

  # called in the initial deep sequenced data
  indels_init_filtered$key1_SampleName <- paste(indels_init_filtered$chrom,
                                                indels_init_filtered$pos,
                                                indels_init_filtered$SampleName,
                                                indels_init_filtered$alt,
                                                indels_init_filtered$SYMBOL,sep=";")

  # called from a specific downsampling run
  indels_caller$key1_SampleName <- paste(indels_caller$chrom,
                                         indels_caller$pos,
                                         indels_caller$SampleName,
                                         indels_caller$alt,
                                         indels_caller$SYMBOL,sep=";")

  # I use the dbSNP, Exac databases to annotate the variants, set the flags. The same variant can be annotated differently, in some cases they have or they don't have tdbSNP for instance. This means that if I use the flag term I will have duplicated variants
  indels_init_filtered_unique <- unique(indels_init_filtered[,c("Location","caller","chrom","pos","ref","alt",
                                                                "key1_SampleName","down_label","Keep_defaults","Keep_annot")])

  indels_caller_unique <- unique(indels_caller[,c("Location","caller","chrom","pos","ref","alt",
                                                  "key1_SampleName","Keep_defaults","Keep_annot")])

  # Create vector
  # is there a match with the variants called in the downsampled dataset?
  # DownMatch_defaults is 1 if a variant in the initial file is called in the downsampled set

  indels_init_filtered_unique <- indels_init_filtered_unique %>%
    mutate(DownMatch_defaults = ifelse(key1_SampleName %in% indels_caller_unique$key1_SampleName[indels_caller_unique$Keep_defaults],TRUE,FALSE)) %>%
    mutate(DownCalled_defaults = ifelse(DownMatch_defaults & Keep_defaults,1,0))

  # same for annot strategy
  indels_init_filtered_unique <- indels_init_filtered_unique %>%
    mutate(DownMatch_annot = ifelse(key1_SampleName %in% indels_caller_unique$key1_SampleName[indels_caller_unique$Keep_annot],TRUE,FALSE)) %>%
    mutate(DownCalled_annot = ifelse(DownMatch_annot & Keep_annot,1,0))

  # Sensitivity
  sens_defaults_indels <- sum(indels_init_filtered_unique$DownCalled_defaults)/sum(indels_init_filtered_unique$Keep_defaults)
  sens_annot_indels <- sum(indels_init_filtered_unique$DownCalled_annot)/sum(indels_init_filtered_unique$Keep_annot)

  sens_indels = c(sens_defaults_indels,sens_annot_indels)
  names(sens_indels) <- c("defaults","annotations")

  # False positives
  # summ all the variants called by a caller in the downsampled set that are not in the initial set
  # 1. only take those variants that matched a variant called in the initial set with default parameters
  # sum the keep_defaults from  in the downsampled dataset
  var_passed_in_initial_defaults <- indels_init_filtered_unique$key1_SampleName[indels_init_filtered_unique$Keep_defaults]
  fpn_defaults_indels <- sum(indels_caller_unique$Keep_defaults[!(indels_caller_unique$key1_SampleName %in% var_passed_in_initial_defaults)])
  fp_defaults_indels <- fpn_defaults_indels/sum(indels_caller_unique$Keep_defaults)

  var_passed_in_initial_annot <- indels_init_filtered_unique$key1_SampleName[indels_init_filtered_unique$Keep_annot]
  fpn_annot_indels <- sum(indels_caller_unique$Keep_annot[!(indels_caller_unique$key1_SampleName %in% var_passed_in_initial_annot)])
  fp_annot_indels <- fpn_annot/sum(indels_caller_unique$Keep_annot)

  fp_indels = c(fpn_defaults_indels,fp_defaults_indels,fpn_annot_indels,fp_annot_indels)
  names(fp_indels) <- c("Ndefaults","Pdefaults","Nannot","Pannot")

  power_init_indels <- list(sens = sens_indels,
                            fp = fp_indels,
                            indels_init_filtered_unique = indels_init_filtered_unique)



  # return results
  list(power=power, # sensitivyt overall with paper variants
       power_init = power_init, # sensitivity using initial variants also to look at false positives
       power_init_indels = power_init_indels,
       variants=coverage_downsampled_expr,
       indels_paper_recovery=indels_paper_recovery, # to be saved and manually curated
       all_variants=variants_down_filtered_expr)
}
