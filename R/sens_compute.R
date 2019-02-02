#' Computes sensitivity with respect to a truth set
#' @param variants dataframe created within the `variant_power.R` function. The columns required in this dataframe are `Called_defaults` and `Called_defaults` which are numeric 0/1 defining whther a variant is called. If `indel = TRUE` also the column `key` is required. `key` is a unique key for a variant: `chr:POS:SampleName:Gene:Transcript`.
#' @param truth_set external truth set. The truth set os only needed to know the total number of variants in it.
#' @param indel logical. TRUE is INDELs are analysed.

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
