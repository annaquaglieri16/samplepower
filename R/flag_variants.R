############
## Parse PON
############

create_vaf <- function(vaf){
  as.numeric(substr(vaf,1,nchar(vaf)-1))/100
}

# parsePON_tsv <- function(linkPON,caller="vardict"){
#   pon <- fread(linkPON,header = TRUE)
#   pon <- data.frame(pon)
#
#   Nsam <- function(set){
#     length(gregexpr("SRX", set)[[1]])
#   }
#
#   pon$nsam <- sapply(pon$set,function(entry_set) Nsam(entry_set))
#   pon$nsam_prop <- pon$nsam/17
#   return(pon)
# }

parsePON <- function(linkPON,caller="vardict"){
  pon <- fread(linkPON,header = TRUE)
  pon <- data.frame(pon)

  # Extract VAF from VarScan output
  VAFvarscan <- function(freq){
    freq <- gsub("%","",freq)
    freq <- as.numeric(as.character(freq))/100
    return(freq)
  }

  VAFother <- function(freq){
    freq <- as.numeric(as.character(freq))
    return(freq)
  }

  # Extract VAF for every variant called by the caller
  # dtermine number of normal samples where eah variant was called
  # determine the minimum VAF that each variant was called with

  if(caller %in% "varscan"){
    FREQcol <- grep("FREQ",colnames(pon))
    pon_FREQ <- pon[,FREQcol]
    pon_VAF <- apply(pon_FREQ,2,function(x) VAFvarscan(x))
    minVAF <- apply(pon_VAF,1,min,na.rm=TRUE)
    combine_vaf <- cbind(pon[,c("CHROM","POS","REF","ALT","set")],minVAF)
  }

  if(caller %in% c("mutect","vardict")){
    FREQcol <- grep(".AF",colnames(pon))
    pon_FREQ <- pon[,FREQcol]
    pon_VAF <- apply(pon_FREQ,2,function(x) VAFother(x))
    minVAF <- apply(pon_VAF,1,min,na.rm=TRUE)
    combine_vaf <- cbind(pon[,c("CHROM","POS","REF","ALT","set")],minVAF)
  }

  # Nset
  Nsam <- function(set){
    length(gregexpr("SRX", set)[[1]])
  }
  combine_vaf$nsam <- sapply(combine_vaf$set,function(entry_set) Nsam(entry_set))
  combine_vaf$nsam_prop <- combine_vaf$nsam/17

  return(combine_vaf)
}


###########################################################
# Extract only fields that I am interested in at the moment
# Byotype is not unique so I would have multiple hits for the same variant]
# Variants should have been VEP annotated
##########################################################

extract_fields <- function(variants,label){

  variants <- variants %>%
    mutate(Existing_variation = ifelse(Existing_variation == "",NA,Existing_variation), # set to NA if no Existing_variation provided
           ExAC_AF = as.numeric(as.character(ExAC_AF)),  # convert to numeric
           EUR_AF = as.numeric(as.character(EUR_AF)), # convert to numeric
           SYMBOL = ifelse(SYMBOL == "",NA,SYMBOL),  # set to NA if no Existing_variation provided
           down_label = label,
           variant_type = case_when(VARIANT_CLASS %in% "SNV" ~ "SNV", # standardise the way in which we define variant_type
                                    !(VARIANT_CLASS %in% "SNV") ~ "INDEL",
                                    TRUE ~ VARIANT_CLASS)) %>%
           filter(!is.na(SYMBOL)) %>%                                      # remove variants not on genes
            unite("key_SampleName", c("chrom", "pos", "SampleName", "SYMBOL", "Feature"), remove =FALSE,sep=":") # create key

  return(unique(variants))

}

#######################################
## Add flag and Keep for every variant
#######################################

flag_variants <- function(variants,
                      normal_variants,
                      flag_patterns,
                      exon_ranges,
                      homop_ranges,
                      RNAedit_ranges,
                      repeats_ranges){


  # Add flags based on a filtering pattern that I created

  # VARIANTS FOUND IN NORMALS
  # If Variants > 2 Normals (= common in normals) OR Variants < 2 normals but with VAF > 3% then define it as found in normals
  # The match in the normals should be based on chrom pos ref and alt. This coulde be biased for INDELs.

  variants <- variants %>% unite(Location, alt , col = "Location_alt", sep="_",remove = FALSE)

  normal_variants <- normal_variants %>% unite(Location, alt , col = "Location_alt", sep="_",remove = FALSE)
  normal_variants_morethan2 <- normal_variants %>% filter(nsam > 2)
  normal_variants_lessthan2 <- normal_variants %>% filter(nsam <= 2 & normal_variants$minVAF > 0.03)

  variants <- variants %>% mutate(PON = case_when( (Location_alt %in% normal_variants_morethan2$Location_alt) |
                                                   (Location_alt %in% normal_variants_lessthan2$Location_alt) ~ 1,
                                                    TRUE ~ 0))

  # 2. COSMIC
  variants$COSMIC <- 0
  variants$COSMIC[grep("COSM",variants$Existing_variation)] <- 1

  # 3. EXAC rare
  variants$EXAC_rare <- ifelse(!is.na(variants$ExAC_AF) & variants$ExAC_AF <= 0.01, 1,0)
  variants$EXAC_common <- ifelse(!is.na(variants$ExAC_AF) & variants$ExAC_AF > 0.01, 1,0)

  # 4. dbSNP
  variants$dbSNP <- 0
  variants$dbSNP[grep("rs",variants$Existing_variation)] <- 1

  # Exon boundaries
  tot_range <- GRanges(seqnames = variants$chrom,
                       IRanges(start=variants$pos,end = variants$pos))

  # Over Exon boundaries
  over_end <- findOverlaps(tot_range,exon_ranges)
  variants$Exon_edge <- 0
  if(length(queryHits(over_end)) > 0){
    variants$Exon_edge[queryHits(over_end)] <- 1
  }

  # Over RNA editing sites
  over_RNAedit <- findOverlaps(tot_range,RNAedit_ranges)
  variants$RADAR <- 0
  if(length(queryHits(over_RNAedit)) > 0){
    variants$RADAR[queryHits(over_RNAedit)] <- 1
  }

  # Repetitive regions
  over_end <- findOverlaps(tot_range,repeats_ranges)
  variants$RepeatMasker <- 0
  if(length(queryHits(over_end)) > 0){
    variants$RepeatMasker[queryHits(over_end)] <- 1
  }

  # Homoplymers > 5
  over_end <- findOverlaps(tot_range,homop_ranges)
  variants$Homopolymers <- 0
  if(length(queryHits(over_end)) > 0){
    variants$Homopolymers[queryHits(over_end)] <- 1
  }

  ########################################################
  # Add quality flag based on default filtering strategies
  ########################################################

  if(as.character(variants$caller)[1] == "freebayes"){
    # QUAL field: mix of both > 20 from github page
    variants$Quality_defaults <- ifelse(variants$qual >= 20, 1,0)
  } else {
    variants$Quality_defaults <- ifelse(variants$filter %in% "PASS", 1,0)
  }

  #########################
  # Add quality flag based
  #########################

  # Vardict only has one MQ quality average mapping quality at that position
  variants <- variants %>%
    mutate(Quality_annot = case_when((qual >= 18) & (alt_depth > 5) & (tot_depth > 15) & (VAF > 0.03) ~ 1,
                                     TRUE ~ 0))

  #########################################################
  ## Add annotation flag and Keep column
  #########################################################

  ## NB: the field Keep is initialised by the PON and annotation databases
  patterns <- read.csv(flag_patterns,stringsAsFactors = FALSE)
  variants_flagged <- merge(variants,patterns,all.x=TRUE)

  # Update Keep filed based on default filtering strategies by the callers
  # 1. Use default filters for the caller + present in normal (PON already considered in the flag patterns csv): Keep if the quality is good and not present in normal
  # 2. Use more stringent heuristic filters + flag with annotation databased

  # Update keep_defaults using PON with default PASS filters from callers
  variants_flagged$Keep_defaults <- ifelse((variants_flagged$Quality_defaults == 0),FALSE, variants_flagged$Keep_defaults)
  # Update keep_annot using PON + annotations with quality filters
  variants_flagged$Keep_annot <- ifelse((variants_flagged$Quality_annot == 0),FALSE,variants_flagged$Keep_annot)


  # Flag based on default filtering by callers and based on our thresholds
  variants_flagged$Flag_defaults <- ifelse(variants_flagged$Quality_defaults == 0,"filtered by caller",as.character(variants_flagged$Flag))
  variants_flagged$Flag_annot <- ifelse(variants_flagged$Quality_annot == 0,"filtered by quality thresholds",as.character(variants_flagged$Flag))

  # Exon edges
  variants_flagged$Keep_annot <- ifelse(variants_flagged$Exon_edge == 1,FALSE,variants_flagged$Keep_annot)
  variants_flagged$Flag_annot <- ifelse(variants_flagged$Exon_edge == 1,"exon_edge",as.character(variants_flagged$Flag_annot))

  # Repetitive regions
  variants_flagged$Keep_annot <- ifelse(variants_flagged$RepeatMasker == 1,FALSE,variants_flagged$Keep_annot)
  variants_flagged$Flag_annot <- ifelse(variants_flagged$RepeatMasker == 1,"repetitive_region",as.character(variants_flagged$Flag_annot))

  # Homopolymers
  variants_flagged$Keep_annot <- ifelse(variants_flagged$Homopolymers == 1,FALSE,variants_flagged$Keep_annot)
  variants_flagged$Flag_annot <- ifelse(variants_flagged$Homopolymers == 1,"homopolymers",as.character(variants_flagged$Flag_annot))

  return(variants_flagged)

}


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
  rpkmCounts <- cpm(counts$counts,log = TRUE)
  rpkmCounts <- rpkmCounts[rownames(rpkmCounts) %in% variants$GeneID,]

  # 4. Convert from wide to long format
  rpkmCounts <- data.frame(rpkmCounts)
  rpkmCounts$GeneID <- rownames(rpkmCounts)
  rpkmCounts_long <- rpkmCounts %>% gather(SampleName,LogRpkm,sample_names)
  rpkmCounts_long$SYMBOL <- variants$SYMBOL[match(rpkmCounts_long$GeneID,variants$GeneID)]

  # 5. Associate the correct logRPKM to each sample
  variant_expr <- merge(variants,rpkmCounts_long,all.x=TRUE)
  return(variant_expr)

}


####################################################################
## True/false positives using the initial bamfiles as gold standard
####################################################################


