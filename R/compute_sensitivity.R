#################
## Outcome Loss
#################


sens_compute <- function(variants,
                         truth_set){

  # Sensitivity for SNVs using default filters + PON
  sens_defaults = sum(variants$Called_defaults)/nrow(truth_set[truth_set$variant_type %in% "SNV",])

  # Sensitivity for SNVs using annotations + quality measures
  sens_annot = sum(variants$Called_annot)/nrow(truth_set[truth_set$variant_type %in% "SNV",])
  sens = c(sens_defaults,sens_annot)
  names(sens) <- c("defaults","annotations")

  return(sens)

}

sensitivity_by_thresholds <- function(variants,
                         truth_set){

  # by different levels of depth
  depth_seq=seq(0,200,5)
  sens_by_depth <- sapply(depth_seq,function(d) {
    sens_compute(variants=subset(variants,tot_depth >= d),
                  truth_set=truth_set)}
  )
  colnames(sens_by_depth) <- depth_seq

  # by different levels of VAF
  vaf_seq=seq(0,1,by = 0.01)
  sens_by_vaf <- sapply(vaf_seq,function(vaf) {
    sens_compute(variants=subset(variants,VAF >= vaf),
                 truth_set=truth_set)}
  )
  colnames(sens_by_vaf) <- vaf_seq

  # overall by the caller
  sens_overall <- sens_compute(variants=variants,
                               truth_set=truth_set)

  # Output variants missed
  variants_missed_defaults <- variants[variants$Called_defaults != 1,]
  variants_missed_annot <- variants[variants$Called_annot != 1,]

  list(sens_by_depth=sens_by_depth,
       sens_by_vaf=sens_by_vaf,
       sens_overall=sens_overall,
       missed_defaults=variants_missed_defaults,
       missed_annot=variants_missed_annot)

}
