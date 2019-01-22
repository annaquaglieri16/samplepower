sens_compute <- function(variants,
                         truth_set,
                         indel = FALSE){

  if(!indel){

    # Sensitivity for SNVs using default filters + PON
    sens_defaults = sum(variants$Called_defaults)/nrow(truth_set)

    # Sensitivity for SNVs using annotations + quality measures
    sens_annot = sum(variants$Called_annot)/nrow(truth_set)
    sens = c(sens_defaults,sens_annot)
    names(sens) <- c("defaults","annotations")

  } else {

    # Sensitivity for SNVs using default filters + PON
    ncalled_annot <- length(unique(variants$key[variants$Called_annot == 1]))
    ncalled_default <- length(unique(variants$key[variants$Called_defaults == 1]))
    ncalled_annot05 <- length(unique(variants$key[variants$Called_annot05 == 1]))
    ncalled_default05 <- length(unique(variants$key[variants$Called_defaults05 == 1]))

    # Sensitivity for SNVs using annotations + quality measures
    sens_annot = ncalled_annot/nrow(truth_set)
    sens_annot05 = ncalled_annot05/nrow(truth_set)
    sens_defaults = ncalled_default/nrow(truth_set)
    sens_defaults05 = ncalled_default05/nrow(truth_set)

    sens = c(sens_defaults,sens_annot,sens_defaults05,sens_annot05)
    names(sens) <- c("defaults","annotations","defaults05","annotations05")

  }

  return(sens)

}
