#' Compute sensitivity and specificity for the a set of variants with respect to some truth sets
#' @param variant_files character vector containing the paths to all the tab delimited files of variants from the downsampled run. See details.
#' @param variant_files_initial character vector containing the paths to all the tab delimited files of variants from the initial run with deeply sequenced RNA-Seq samples. See details.
#' @param down_label character. Any label to be assigned to the current run
#' @param truth_set external truth set
#' @param caller one of `varscan`, `mutect` or `vardict`
#' @param gene_expression if available, `rds` objects with raw gene expression data where rows are genes and columns are samples
#' @param normal_variants link to panel of normal file parsed with `samplepower::parse_pon`
#' @param exon_ranges GRanges object with exon boundaries relative to the reference genome used
#' @param homop_ranges GRanges object with homopolymers regions of more than 5bp relative to the reference genome used
#' @param RNAedit_ranges GRanges object with RNA editing sites relative to the reference genome used
#' @param repeats_ranges GRanges object with repetitive regions of the reference genome used
#' @param ncbi NCBI annotation dataframe including at least the columns `GeneID` and `Symbol`
#' @param flag_patterns dataframe containing the flags to assign to variants. THIS NEEDS TO BE SET AS INTERNAL TO THE PACKAGE
#' @param path_to_gatk_coverage path to file where the GATK `DepthOfCoverage` output is saved
#' @param TCGA logical. If TRUE variants will be macthed with truth sets excluding transcript information
#'
#' @details
#' The files whose links are recorded in `variant_files_initial` and `variant_files` are the ouput files obtained with the `call_variants.R` function from the functions in https://github.com/annaquaglieri16/RNA-seq-variant-calling#panel-of-normals-pon. The VCF output from every callers are standardised to make the results easily comparable across callers.


variants_power <- function(variant_files,
                           variant_files_initial,
                           down_label = NA,
                           truth_set = NA,
                           caller,
                           gene_expression = NA,
                           normal_variants,
                           exon_ranges,
                           homop_ranges,
                           RNAedit_ranges,
                           repeats_ranges,
                           ncbi,
                           flag_patterns,
                           path_to_gatk_coverage,
                           TCGA=FALSE){

  #########################################
  # 0. Check column names and validity of files
  #########################################

  print(paste0("Run: ",down_label, " Caller: ",caller))
  print(Sys.time())

  if(length(variant_files) < 1){
    stop("No downsampled variant files available in input")
  }

  if(length(variant_files_initial) < 1){
    stop("No initial variant files available in input")
  }


  search_env <- pryr::where("down_label")
  if(!exists_with_class( "down_label","character",where = search_env)){
    stop("down_label (label for the downsampled run) was not specified")
  }

  # Truth set - existence and column names
  if(!exists("truth_set",where = search_env)){
    stop("No set of true variants provided.")
  } else {

    need_colmuns <- c("chrom","pos","Locus","alt_initial","ref_initial","variant_type","SYMBOL","Feature","SampleName")

    check_columns <- sum(!(need_colmuns %in% colnames(truth_set)))

    if(check_columns > 0){
      missing <- need_colmuns[!(need_colmuns %in% colnames(truth_set))]
      stop(paste0("Check requirements for column names of truth_set. The following columns are missing: ",missing))
    }
  }

  if(is.na(caller)){
    stop("Specify which 'caller' was used to produce the variants provided.")
  }

  #if(!exists("gene_expression",where = search_env)){
  #  add_gene_counts <- FALSE
  #  warning("The path to the gene_expression file is missing. The logRPKM of genes won't be added next to the variants reported.")
  #} else {
  #    if(!file.exists(gene_expression)){
  #      stop("Path to gene_expression provided but file does not exist.")
  #    }else{
  #      add_gene_counts = TRUE
  #    }
  #}

  if(is.na(gene_expression)){
    add_gene_counts <- FALSE
    warning("The path to the gene_expression file is missing. The logRPKM of genes won't be added next to the variants reported.")
  } else {
    if(!file.exists(gene_expression)){
      stop("Path to gene_expression provided but file does not exist.")
    }else{
      add_gene_counts = TRUE
    }
  }


  # existence and column names
  if(class(try(nrow(normal_variants))) %in% "try-error"){
    stop("Panel of normal variants not provided in normal_variants argument")
  }else{
    check_columns <-
      sum(!(c("nsam","minVAF") %in% colnames(normal_variants)))

    if(check_columns > 0){
      missing <- colnames(normal_variants)[!(c("nsam","minVAF") %in% colnames(normal_variants))]
      stop(paste0("Check requirements for column names of normal_variants The following columns are missing: ",missing))
    }
  }

  try_open1 <- exists_with_class(obj = "exon_ranges", check_class =  "GRanges",silent = TRUE,where = search_env)
  try_open2 <- exists_with_class(obj = "homop_ranges", check_class =  "GRanges",silent = TRUE,where = search_env)
  try_open3 <- exists_with_class(obj = "RNAedit_ranges", check_class =  "GRanges",silent = TRUE,where = search_env)
  try_open4 <- exists_with_class(obj = "repeats_ranges", check_class =  "GRanges",silent = TRUE,where = search_env)
  try_open <- c(try_open1,try_open2,try_open3,try_open4)
  names(try_open) <- c("exon_ranges","homop_ranges","RNAedit_ranges","repeats_ranges")

  if(!try_open1 | !try_open2 | !try_open3 | !try_open4){
    stop(paste0("The GRanges annotation objects are missing for: ",
                paste0(names(try_open)[c(!try_open1,!try_open2,!try_open3,!try_open4)],collapse=",")," are missing."))
  }

  #flag_patterns

  if(!exists("path_to_gatk_coverage")){
    add_gatk_dop <- FALSE
    warning("The path to the path_to_gatk_coverage file is missing. The alt and tot depth from GATK depth of coverage won't be reported.")
  } else {
    if(!file.exists(path_to_gatk_coverage)){
      stop("Path to path_to_gatk_coverage provided but file does not exist.")
    }else{
      add_gatk_dop = TRUE
    }
  }

  if(!exists("ncbi")){
    stop("The NCBI dataset providing GeneIDs and SYMBOLs is missing")
  }else{
    check_columns <-
      sum(!(c("GeneID","SYMBOL") %in% colnames(ncbi)))

    if(check_columns > 0){
      needs_column <- c("GeneID","SYMBOL")
      missing <- needs_column[!(needs_column %in% c("GeneID","SYMBOL"))]
      stop(paste0("Check requirements for column names of NCBI. The following columns are missing: ",missing))
    }
  }

  #################################
  # 1. Read in Downsampled variants
  #################################

  # Germline - read variants in for the run
  names_samples <- gsub("_germline_final.txt","",basename(variant_files))
  list_variants <- lapply(variant_files,function(x){
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
                                extract_fields(list_variants[[index]],label=down_label,use_transcript = !TCGA)}
  )
  down_variants_sub <- do.call(rbind,down_variants_sub)

  # Only consider canonical chromosomes
  chroms <- c(paste0("chr",1:22),"chrX","chrY","chrM")
  down_variants_sub <- down_variants_sub %>% dplyr::filter(chrom %in% chroms)

  #############################
  # 1. Read in Initial variants
  #############################

  names_samples <- gsub("_germline_final.txt","",basename(variant_files_initial))
  list_variants <- lapply(variant_files_initial,
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
                            extract_fields(list_variants[[index]],label=down_label,use_transcript = !TCGA)}
  )
  variants_init <- do.call(rbind,variants_init)

  # Only consider canonical chromosomes
  chroms <- c(paste0("chr",1:22),"chrX","chrY","chrM")
  variants_init <- variants_init %>% dplyr::filter(chrom %in% chroms)

  print("Important fields extracted.")

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
                                          RNAedit_ranges = RNAedit_ranges,
                                          repeats_ranges = repeats_ranges)

  print("Downsampled variant flagged.")

  # Only select variants falling on genes that we have to study
  variants_down_filtered <- subset(variants_down_filtered, SYMBOL %in% truth_set$SYMBOL)

  # Initial
  variants_init_filtered <- flag_variants(variants = variants_init,
                                          normal_variants = normal_variants,
                                          flag_patterns = flag_patterns,
                                          exon_ranges = exon_ranges,
                                          homop_ranges = homop_ranges,
                                          RNAedit_ranges = RNAedit_ranges,
                                          repeats_ranges = repeats_ranges)

  variants_init_filtered <- subset(variants_init_filtered, SYMBOL %in% truth_set$SYMBOL)

  print("Initial variant flagged.")

  # Truth set: Variants in paper
  # 3. Read in coverage computed for every Location - only for the 58 SNVs in the paper
  # Downsampled

  coverage_downsampled <- parse_gatk_coverage(truth_set = truth_set[truth_set$variant_type %in% "SNV",],
                                              path_to_gatk_coverage = as.character(path_to_gatk_coverage),
                                              TCGA=TCGA)

  # VAF_GATK will be from GATK depth of Cov and VAF the estimate from each caller
  if(TCGA){

    #tab <- variants_down_filtered %>% group_by(key_SampleName) %>% count()

    variants_down_filtered <- unique(variants_down_filtered %>%
                                       dplyr::select(SampleName,chrom,pos,ref,alt,qual,filter,genotype,
                                              tot_depth ,VAF,ref_depth,key_SampleName,
                                              alt_depth,ref_forw,ref_rev,alt_forw,alt_rev,
                                              SYMBOL,down_label,variant_type,Exon_edge,RepeatMasker,Homopolymers,
                                              Quality_defaults,Quality_annot,Flag,Keep_annot,Keep_defaults,
                                              germline_somatic,Flag_defaults,Flag_annot,
                                              PON,COSMIC,EXAC_rare,EXAC_common,dbSNP,RADAR,Location,caller))
    #tab <- variants_down_filtered1 %>% group_by(key_SampleName) %>% count()

    variants_init_filtered <- unique(variants_init_filtered %>%
                                       dplyr::select(SampleName,chrom,pos,ref,alt,qual,filter,
                                              genotype,tot_depth ,VAF,ref_depth,key_SampleName,
                                      alt_depth,ref_forw,ref_rev,alt_forw,alt_rev,SYMBOL,down_label,
                                      variant_type,Exon_edge,RepeatMasker,Homopolymers,
                                      Quality_defaults,Quality_annot,Flag,Keep_annot,Keep_defaults,
                                      germline_somatic,Flag_defaults,Flag_annot,
                                     PON,COSMIC,EXAC_rare,EXAC_common,dbSNP,RADAR,Location,caller))

    coverage_downsampled <- unique(merge(coverage_downsampled,
                                         variants_down_filtered,all.x=TRUE))

  }else{

    coverage_downsampled <- unique(merge(coverage_downsampled,
                                         variants_down_filtered,all.x=TRUE))
  }

  print("Added GATK total and alt depth.")

  # Add label for downsampling run and caller
  coverage_downsampled$down_label <- down_label
  coverage_downsampled$caller <- unique(variants_down_filtered$caller)

  ###################################################
  # Define if a SNV was called or not in this run:
  # Using defualt and annot filters as well as the match with the alt allele
  ###################################################

  coverage_downsampled <- coverage_downsampled %>%
    dplyr::mutate(Called_annot = ifelse(!is.na(Keep_annot) & Keep_annot,1,0),
                  Called_defaults = ifelse(!is.na(Keep_defaults) & Keep_defaults,1,0),
                  match_alt = ifelse(alt == alt_initial,1,0)) %>%
    dplyr::mutate(Flag_annot = case_when(match_alt == 0 ~ paste0(Flag_annot,";mismatch alt"),
                                         TRUE ~ Flag_annot))

  # Alt is obtained from GATK both if the variant is called or not called so there will always an ALT
  coverage_downsampled <- coverage_downsampled %>%
    dplyr::mutate(Called_annot =  ifelse(Called_annot == 1 & match_alt == 1,1,0),
                  Called_defaults = ifelse(Called_defaults == 1 & match_alt == 1,1,0))

  # 4. Add gene expression counts
  coverage_downsampled$GeneID <- ncbi$GeneID[match(coverage_downsampled$SYMBOL,ncbi$SYMBOL)]
  sampleNames <- as.character(unique(truth_set$SampleName))

  if(add_gene_counts){
    coverage_downsampled_expr <- add_log_rpkm(variants = coverage_downsampled,
                                              gene_expression = gene_expression,
                                              sample_names = sampleNames)
  }else{
    coverage_downsampled_expr <- coverage_downsampled
  }


  # 5. Determine power of recovery
  power <- sensitivity_by_thresholds(variants = coverage_downsampled,
                                     truth_set = truth_set[truth_set$variant_type %in% "SNV",])

  print("Sensitivity computed for SNVs.")

  #############################################
  ## Sens and FP using initial set as truth set
  #############################################

  variants_init_filtered$key1_SampleName <- paste(variants_init_filtered$chrom,
                                                  variants_init_filtered$pos,
                                                  variants_init_filtered$SampleName,
                                                  variants_init_filtered$alt,
                                                  variants_init_filtered$SYMBOL,sep=":")

  variants_down_filtered$key1_SampleName <- paste(variants_down_filtered$chrom,
                                                  variants_down_filtered$pos,
                                                  variants_down_filtered$SampleName,
                                                  variants_down_filtered$alt,
                                                  variants_down_filtered$SYMBOL,sep=":")

  # I use the dbSNP, Exac databases to annotate the variants and set the flags.
  # The same variant can be annotated differently, in some cases they have or they don't have tdbSNP for instance.
  # This means that if I use the flag term I will have duplicated variants
  # What I am interested here is to see how much different library sizes affect the false positive rate at different levels of downsampling.
  #variants_init_filtered_unique <- unique(variants_init_filtered[,c("Location","caller","chrom","pos","ref","alt", "key1_SampleName","down_label"#,"Keep_defaults","Keep_annot")])

  #  variants_down_filtered_unique <- unique(variants_down_filtered[,c("Location","caller","chrom","pos","ref","alt",
  #                                                                    "key1_SampleName","Keep_defaults","Keep_annot")])


  # Create unique rows for each variant by pasting the flags
  variants_init_filtered_unique <- variants_init_filtered %>%
    dplyr::group_by(key1_SampleName,chrom,pos,SampleName,alt,SYMBOL,caller,down_label,variant_type) %>%
    dplyr::filter(variant_type %in% "SNV") %>%
    dplyr::summarise(Flag_annot = paste(Flag_annot,collapse = ";"),
                     Flag_defaults = paste(Flag_defaults,collapse = ";"),
                     Keep_annot = max(Keep_annot),
                     Keep_defaults = max(Keep_defaults),
                     qual = round(mean(qual,na.rm=TRUE),3),
                     ref_depth = round(mean(ref_depth,na.rm=TRUE),3),
                     alt_depth = round(mean(alt_depth,na.rm=TRUE),3),
                     tot_depth = round(mean(tot_depth,na.rm=TRUE),3),
                     VAF = round(mean(VAF,na.rm=TRUE),3)) %>%
    dplyr::mutate(Keep_annot = ifelse(Keep_annot == 1, TRUE,FALSE),
                  Keep_defaults = ifelse(Keep_defaults == 1, TRUE,FALSE),
                  Flag_defaults = sapply(strsplit(Flag_defaults,split=";"),function(x) paste(sort(unique(x)),collapse=";")),
                  Flag_annot = sapply(strsplit(Flag_annot,split=";"),function(x) paste(sort(unique(x)),collapse=";")))


  # dowsampled
  variants_down_filtered_unique <- variants_down_filtered %>%
    dplyr::filter(variant_type %in% "SNV") %>%
    dplyr::group_by(key1_SampleName,chrom,pos,SampleName,alt,SYMBOL,caller,down_label,variant_type) %>%
    dplyr::summarise(Flag_annot = paste(Flag_annot,collapse = ";"),
                     Flag_defaults = paste(Flag_defaults,collapse = ";"),
                     Keep_annot = max(Keep_annot),
                     Keep_defaults = max(Keep_defaults),
                     qual = round(mean(qual,na.rm=TRUE),3),
                     ref_depth = round(mean(ref_depth,na.rm=TRUE),3),
                     alt_depth = round(mean(alt_depth,na.rm=TRUE),3),
                     tot_depth = round(mean(tot_depth,na.rm=TRUE),3),
                     VAF = round(mean(VAF,na.rm=TRUE),3)) %>%
    dplyr::rename(alt_depth_down = alt_depth,
                  VAF_down = VAF,
                  qual_down = qual,
                  Flag_annot_down = Flag_annot,
                  Flag_defaults_down = Flag_defaults,
                  Keep_annot_down = Keep_annot,
                  Keep_defaults_down = Keep_defaults) %>%
    dplyr::mutate(Keep_annot_down = ifelse(Keep_annot_down == 1, TRUE,FALSE),
                  Keep_defaults_down = ifelse(Keep_defaults_down == 1, TRUE,FALSE),
                  Flag_defaults_down = sapply(strsplit(Flag_defaults_down,split=";"),function(x) paste(sort(unique(x)),collapse=";")),
                  Flag_annot_down = sapply(strsplit(Flag_annot_down,split=";"),function(x) paste(sort(unique(x)),collapse=";")))

  # Match if a variant present in the initial run was found in the downsampled dataset.
  # Is there a match with the variants called in the downsampled dataset?
  # DownMatch_defaults is 1 if a variant in the initial file is called in the downsampled set
  #
    variant_called_down <- variants_down_filtered_unique$key1_SampleName[variants_down_filtered_unique$Keep_defaults_down]

  variants_init_filtered_unique <- variants_init_filtered_unique %>%
    dplyr::mutate(DownMatch_defaults = ifelse(key1_SampleName %in% variant_called_down,TRUE,FALSE)) %>%
    dplyr::mutate(DownCalled_defaults = ifelse(DownMatch_defaults & Keep_defaults,1,0))

  # Same for annot strategy
  variant_called_down <- variants_down_filtered_unique$key1_SampleName[variants_down_filtered_unique$Keep_annot_down]

  variants_init_filtered_unique <- variants_init_filtered_unique %>%
    dplyr::mutate(DownMatch_annot = ifelse(key1_SampleName %in% variant_called_down,TRUE,FALSE)) %>%
    dplyr::mutate(DownCalled_annot = ifelse(DownMatch_annot & Keep_annot,1,0))

  ######################
  # Compute Sensitivity
  #####################

  sens_defaults <- sum(variants_init_filtered_unique$DownCalled_defaults)/sum(variants_init_filtered_unique$Keep_defaults)
  sens_annot <- sum(variants_init_filtered_unique$DownCalled_annot)/sum(variants_init_filtered_unique$Keep_annot)

  sens = c(sens_defaults,sens_annot)
  names(sens) <- c("defaults","annotations")

  print("Sensitivity using the initial set as truth set computed for SNVs.")

  #############
  # Compute FP
  #############
  # Sum all the variants called by a caller in the downsampled set that are not in the initial set
  # 1. only take those variants that matched a variant called in the initial set with default parameters
  # Sum the keep_defaults from the downsampled dataset
  variant_called_init <- variants_init_filtered_unique$key1_SampleName[variants_init_filtered_unique$Keep_defaults]

  # Number of variants called in the downsmapled call set not present in the initial call set
  fpn_defaults <- sum(variants_down_filtered_unique$Keep_defaults_down[!(variants_down_filtered_unique$key1_SampleName %in% variant_called_init)])
  fp_defaults <- fpn_defaults/sum(variants_down_filtered_unique$Keep_defaults_down)

  # annot strategy
  variant_called_init <- variants_init_filtered_unique$key1_SampleName[variants_init_filtered_unique$Keep_annot]

  fpn_annot <- sum(variants_down_filtered_unique$Keep_annot_down[!(variants_down_filtered_unique$key1_SampleName %in% variant_called_init)])
  fp_annot <- fpn_annot/sum(variants_down_filtered_unique$Keep_annot_down)

  fp = c(fpn_defaults,fp_defaults,fpn_annot,fp_annot)
  names(fp) <- c("Ndefaults","Pdefaults","Nannot","Pannot")

  # Merge variants found in the initial dataset with the features (VAF/alt and tot depth) from the ones in the downsampled set
  # I will use this dataset to study the features of the variants lost in downsampled runs

  variants_init_filtered_unique <- variants_init_filtered_unique %>%
    dplyr::left_join(variants_down_filtered_unique[,c("key1_SampleName","chrom","pos","SampleName",
                                               "VAF_down","alt_depth_down",
                                               "SYMBOL","caller","down_label",
                                               "Flag_annot_down","Flag_defaults_down",
                                               "Keep_annot_down","Keep_defaults_down","variant_type")]) %>%
    dplyr::mutate(Flag_annot_down = ifelse(is.na(Flag_annot_down),"Not called",Flag_annot_down),
                  Flag_defaults_down = ifelse(is.na(Flag_defaults_down),"Not called",Flag_defaults_down))

  power_init <- list(sens = sens,
                     fp = fp,
                     variants_init_filtered_unique = variants_init_filtered_unique)

  print("False Positive rate using the initial set as truth set computed for SNVs.")


  #############################################################
  # Add gene expression for all variants called by this caller
  ##############################################################

  # Add gene ID and expression

  if ( add_gene_counts ) {
    variants_down_filtered$GeneID <- ncbi$GeneID[match(variants_down_filtered$SYMBOL,ncbi$SYMBOL)]
    variants_down_filtered_expr <- add_log_rpkm(variants = variants_down_filtered,
                                                gene_expression = gene_expression,
                                                sample_names = sampleNames)
  } else {
    variants_down_filtered_expr <- variants_down_filtered
  }

  print("Gene expression added.")

  #####################
  # Recovery of INDELS
  #####################

  # Select everything thet it is not a SNV: deletions, insetions. indels, substitutions
  indels_paper <- subset(truth_set,variant_type %in% "INDEL")
  indels_caller <- subset(variants_down_filtered_expr,!(variant_type %in% "SNV"))

  # Check if the variant passes the annotation filter or the default filters: Keep it only if it passes
  indels_caller$Called_annot <- ifelse(!is.na(indels_caller$Keep_annot) & indels_caller$Keep_annot,1,0)
  indels_caller$Called_defaults <- ifelse(!is.na(indels_caller$Keep_defaults) & indels_caller$Keep_defaults,1,0)

  # See how many variants from the truth set are recovered
  indels_paper_recovery <- match_indels(truth_set = indels_paper,
                                        variants = indels_caller,
                                        use_transcript = !TCGA)

  print("Indels matched with external truth set.")

  # The sensitivity will be donw with Manual curation to see if the INDELs were called and outside of the function I will compute the sensitivity.

  # Reorder columns
  if(!TCGA){
    indels_paper_recovery <- indels_paper_recovery %>%
      dplyr::select(chrom,pos,ref,alt,ref_initial,alt_initial,Mutation,Location,indel_features,VAF,VAF_paper,ref_depth,alt_depth,SYMBOL,
                    type_indel,type_indel2,caller,down_label,key_SampleName,SampleName,everything())
  } else {

    indels_paper_recovery <- indels_paper_recovery %>%
      dplyr::select(chrom,pos,ref,alt,ref_initial,alt_initial,Location,VAF,ref_depth,alt_depth,SYMBOL,
                    caller,down_label,key_SampleName,SampleName,everything())

  }


  #############################################################
  ## Sens and FP using initial set as truth for INDELs for which I don't need manual curation
  #############################################################

  indels_init_filtered <- subset(variants_init_filtered,!(variant_type %in% "SNV"))

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
#  indels_init_filtered_unique <- unique(indels_init_filtered[,c("Location","caller","chrom","pos","ref","alt",                                                                "key1_SampleName","down_label"#,"Keep_defaults","Keep_annot")])

#  indels_caller_unique <- unique(indels_caller[,c("Location","caller","chrom","pos","ref","alt",
 #                                                 "key1_SampleName","Keep_defaults","Keep_annot")])



  # Create unique rows for each variant by pasting the flags
  indels_init_filtered_unique <- indels_init_filtered %>%
    dplyr::group_by(key1_SampleName,chrom,pos,SampleName,alt,SYMBOL,caller,down_label) %>%
    dplyr::summarise(Flag_annot = paste(Flag_annot,collapse = ";"),
              Flag_defaults = paste(Flag_defaults,collapse = ";"),
              Keep_annot = max(Keep_annot),
              Keep_defaults = max(Keep_defaults),
              qual = round(mean(qual,na.rm=TRUE),3),
              ref_depth = round(mean(ref_depth,na.rm=TRUE),3),
              alt_depth = round(mean(alt_depth,na.rm=TRUE),3),
              tot_depth = round(mean(tot_depth,na.rm=TRUE),3),
              VAF = round(mean(VAF,na.rm=TRUE),3)) %>%
    dplyr::mutate(Keep_annot = ifelse(Keep_annot == 1, TRUE,FALSE),
           Keep_defaults = ifelse(Keep_defaults == 1, TRUE,FALSE),
           Flag_defaults = sapply(strsplit(Flag_defaults,split=";"),function(x) paste(sort(unique(x)),collapse=";")),
           Flag_annot = sapply(strsplit(Flag_annot,split=";"),function(x) paste(sort(unique(x)),collapse=";")))

  indels_caller_unique <- indels_caller %>%
    dplyr::group_by(key1_SampleName,chrom,pos,SampleName,alt,SYMBOL,caller,down_label) %>%
    dplyr::summarise(Flag_annot = paste(Flag_annot,collapse = ";"),
              Flag_defaults = paste(Flag_defaults,collapse = ";"),
              Keep_annot = max(Keep_annot),
              Keep_defaults = max(Keep_defaults),
              qual = round(mean(qual,na.rm=TRUE),3),
              ref_depth = round(mean(ref_depth,na.rm=TRUE),3),
              alt_depth = round(mean(alt_depth,na.rm=TRUE),3),
              tot_depth = round(mean(tot_depth,na.rm=TRUE),3),
              VAF = round(mean(VAF,na.rm=TRUE),3)) %>%
    dplyr::rename(alt_depth_down = alt_depth,
                  VAF_down = VAF,
                  qual_down = qual,
                  Flag_annot_down = Flag_annot,
                  Flag_defaults_down = Flag_defaults,
                  Keep_annot_down = Keep_annot,
                  Keep_defaults_down = Keep_defaults) %>%
    dplyr::mutate(Keep_annot_down = ifelse(Keep_annot_down == 1, TRUE,FALSE),
           Keep_defaults_down = ifelse(Keep_defaults_down == 1, TRUE,FALSE),
           Flag_defaults_down = sapply(strsplit(Flag_defaults_down,split=";"),function(x) paste(sort(unique(x)),collapse=";")),
           Flag_annot_down = sapply(strsplit(Flag_annot_down,split=";"),function(x) paste(sort(unique(x)),collapse=";")))

  # Create vector
  # is there a match with the variants called in the downsampled dataset?
  # DownMatch_defaults is 1 if a variant in the initial file is called in the downsampled set

  indels_init_filtered_unique <- indels_init_filtered_unique %>%
    dplyr::mutate(DownMatch_defaults = ifelse(key1_SampleName %in% indels_caller_unique$key1_SampleName[indels_caller_unique$Keep_defaults_down],TRUE,FALSE)) %>%
    dplyr::mutate(DownCalled_defaults = ifelse(DownMatch_defaults & Keep_defaults,1,0))

  # same for annot strategy
  indels_init_filtered_unique <- indels_init_filtered_unique %>%
    dplyr::mutate(DownMatch_annot = ifelse(key1_SampleName %in% indels_caller_unique$key1_SampleName[indels_caller_unique$Keep_annot_down],TRUE,FALSE)) %>%
    dplyr::mutate(DownCalled_annot = ifelse(DownMatch_annot & Keep_annot,1,0))

  # Sensitivity
  sens_defaults_indels <- sum(indels_init_filtered_unique$DownCalled_defaults)/sum(indels_init_filtered_unique$Keep_defaults)
  sens_annot_indels <- sum(indels_init_filtered_unique$DownCalled_annot)/sum(indels_init_filtered_unique$Keep_annot)

  sens_indels = c(sens_defaults_indels,sens_annot_indels)
  names(sens_indels) <- c("defaults","annotations")

  print("Sensitivity using the initial set as truth set computed for INDELs.")

  # False positives
  # summ all the variants called by a caller in the downsampled set that are not in the initial set
  # 1. only take those variants that matched a variant called in the initial set with default parameters
  # sum the keep_defaults from  in the downsampled dataset

  # Using default filters
  var_passed_in_initial_defaults <- indels_init_filtered_unique$key1_SampleName[indels_init_filtered_unique$Keep_defaults]
  fpn_defaults_indels <- sum(indels_caller_unique$Keep_defaults_down[!(indels_caller_unique$key1_SampleName %in% var_passed_in_initial_defaults)])
  fp_defaults_indels <- fpn_defaults_indels/sum(indels_caller_unique$Keep_defaults_down)

  # using filters which include external databases
  var_passed_in_initial_annot <- indels_init_filtered_unique$key1_SampleName[indels_init_filtered_unique$Keep_annot]
  fpn_annot_indels <- sum(indels_caller_unique$Keep_annot_down[!(indels_caller_unique$key1_SampleName %in% var_passed_in_initial_annot)])
  fp_annot_indels <- fpn_annot_indels/sum(indels_caller_unique$Keep_annot_down)

  fp_indels = c(fpn_defaults_indels,fp_defaults_indels,fpn_annot_indels,fp_annot_indels)
  names(fp_indels) <- c("Ndefaults","Pdefaults","Nannot","Pannot")

 # Merge initial and downsamples runs
  indels_init_filtered_unique <- indels_init_filtered_unique %>%
    dplyr::left_join(indels_caller_unique[,c("key1_SampleName","chrom","pos","SampleName",
                                               "VAF_down","alt_depth_down",
                                               "SYMBOL","caller","down_label",
                                               "Flag_annot_down","Flag_defaults_down",
                                               "Keep_annot_down","Keep_defaults_down")]) %>%
    dplyr::mutate(Flag_annot_down = ifelse(is.na(Flag_annot_down),"Not called",Flag_annot_down),
                  Flag_defaults_down = ifelse(is.na(Flag_defaults_down),"Not called",Flag_defaults_down))

  power_init_indels <- list(sens = sens_indels,
                            fp = fp_indels,
                            indels_init_filtered_unique = indels_init_filtered_unique)

  print("FP rate using the initial set as truth set computed for INDELs.")

  ######################################################
  # Add flag in the variant set: 1 if FP and 0 otherwise
  ######################################################

  # I will use the dataset below to study the features of the false positives in the downsampled runs
  variants_down_filtered_expr$FP_defaults <- ifelse(!(variants_down_filtered_expr$key1_SampleName %in% variants_init_filtered$key1_SampleName) &
                                                      variants_down_filtered_expr$Keep_defaults, 1, 0)
  variants_down_filtered_expr$FP_annot <- ifelse(!(variants_down_filtered_expr$key1_SampleName %in% variants_init_filtered$key1_SampleName) &
                                                   variants_down_filtered_expr$Keep_annot, 1, 0)


  print("Analysis completed")

  # return results
  list(power=power, # sensitivyt overall with paper variants
       power_init = power_init, # sensitivity using initial variants also to look at false positives
       power_init_indels = power_init_indels,
       variants=coverage_downsampled_expr,
       indels_paper_recovery=indels_paper_recovery, # to be saved and manually curated
       all_variants=variants_down_filtered_expr)
}

