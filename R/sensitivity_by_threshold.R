#' This function is a wrapper for `sens_compute()` and computes sensitivity for varying VAF and total depth.
#' @param variants dataframe created within the `variant_power.R` function. The columns required in this dataframe are `Called_defaults` and `Called_defaults` which are numeric 0/1 defining whther a variant is called. If `indel = TRUE` also the column `key` is required. `key` is a unique key for a variant: `chr:POS:SampleName:Gene:Transcript`.
#' @param truth_set external truth set
#' @param indel logical. TRUE is INDELs are analysed.

sensitivity_by_thresholds <- function(variants,
                                      truth_set,indel = FALSE){

  # by different levels of depth
  depth_seq=seq(0,200,5)
  sens_by_depth <- sapply(depth_seq,function(d) {
    sens_compute(variants =subset(variants,tot_depth >= d),
                 truth_set = truth_set,indel = indel)}
  )
  colnames(sens_by_depth) <- depth_seq

  # by different levels of VAF
  vaf_seq=seq(0,1,by = 0.01)
  sens_by_vaf <- sapply(vaf_seq,function(vaf) {
    sens_compute(variants=subset(variants,VAF >= vaf),
                 truth_set=truth_set,indel = indel)}
  )
  colnames(sens_by_vaf) <- vaf_seq

  # overall by the caller
  sens_overall <- sens_compute(variants=variants,
                               truth_set=truth_set,indel = indel)

  # Output variants missed
  variants_missed_defaults <- variants[variants$Called_defaults != 1,]
  variants_missed_annot <- variants[variants$Called_annot != 1,]

  list(sens_by_depth=sens_by_depth,
       sens_by_vaf=sens_by_vaf,
       sens_overall=sens_overall,
       missed_defaults=variants_missed_defaults,
       missed_annot=variants_missed_annot)

}
