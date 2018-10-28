sens_compute <- function(variants,
                         truth_set){

  # Sensitivity for SNVs using default filters + PON
  sens_defaults = sum(variants$Called_defaults)/nrow(truth_set)

  # Sensitivity for SNVs using annotations + quality measures
  sens_annot = sum(variants$Called_annot)/nrow(truth_set)
  sens = c(sens_defaults,sens_annot)
  names(sens) <- c("defaults","annotations")

  return(sens)

}
