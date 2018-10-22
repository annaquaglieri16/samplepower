

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

  variants <- variants %>% tidyr::unite(Location, alt , col = "Location_alt", sep="_",remove = FALSE)

  normal_variants <- normal_variants %>% tidyr::unite(Location, alt , col = "Location_alt", sep="_",remove = FALSE)
  normal_variants_morethan2 <- normal_variants %>% dplyr::filter(nsam > 2)
  normal_variants_lessthan2 <- normal_variants %>% dplyr::filter(nsam <= 2 & normal_variants$minVAF > 0.03)

  variants <- variants %>%
    dplyr::mutate(PON = dplyr::case_when( (Location_alt %in% normal_variants_morethan2$Location_alt) |
                                                   (Location_alt %in% normal_variants_lessthan2$Location_alt) ~ 1,
                                                    TRUE ~ 0))

  print("Variants present in normals flagged.")

  # 2. COSMIC
  variants$COSMIC <- 0
  variants$COSMIC[grep("COSM",variants$Existing_variation)] <- 1
  print("Variants present in COSMIC flagged.")


  # 3. EXAC rare
  variants$EXAC_rare <- ifelse(!is.na(variants$ExAC_AF) & variants$ExAC_AF <= 0.01, 1,0)
  variants$EXAC_common <- ifelse(!is.na(variants$ExAC_AF) & variants$ExAC_AF > 0.01, 1,0)
  print("Rare and common ExAC variants flagged.")

  # 4. dbSNP
  variants$dbSNP <- 0
  variants$dbSNP[grep("rs",variants$Existing_variation)] <- 1
  print("dbSNP variants flagged.")

  # Exon boundaries
  tot_range <- GenomicRanges::GRanges(seqnames = variants$chrom,
                                      IRanges::IRanges(start=variants$pos,end = variants$pos))

  # Over Exon boundaries
  over_end <- GenomicRanges::findOverlaps(tot_range,exon_ranges)
  variants$Exon_edge <- 0
  if(length(S4Vectors::queryHits(over_end)) > 0){
    variants$Exon_edge[S4Vectors::queryHits(over_end)] <- 1
  }
  print("Variants over exon boundaries flagged.")


  # Over RNA editing sites
  over_RNAedit <- GenomicRanges::findOverlaps(tot_range,RNAedit_ranges)
  variants$RADAR <- 0
  if(length(S4Vectors::queryHits(over_RNAedit)) > 0){
    variants$RADAR[S4Vectors::queryHits(over_RNAedit)] <- 1
  }
  print("Variants over RNA editing sites flagged.")

  # Repetitive regions
  over_end <- GenomicRanges::findOverlaps(tot_range,repeats_ranges)
  variants$RepeatMasker <- 0
  if(length(S4Vectors::queryHits(over_end)) > 0){
    variants$RepeatMasker[S4Vectors::queryHits(over_end)] <- 1
  }
  print("Variants over repetitive regions flagged.")

  # Homoplymers > 5
  over_end <- GenomicRanges::findOverlaps(tot_range,homop_ranges)
  variants$Homopolymers <- 0
  if(length(S4Vectors::queryHits(over_end)) > 0){
    variants$Homopolymers[S4Vectors::queryHits(over_end)] <- 1
  }
  print("Variants over homopolymer regions > 5bp flagged.")

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
    dplyr::mutate(Quality_annot = dplyr::case_when((qual >= 18) & (alt_depth > 5) & (tot_depth > 15) & (VAF > 0.03) ~ 1,
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


