#################
## Outcome Loss
#################



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
