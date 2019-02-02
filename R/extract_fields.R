#' Extract specific fields from `variant_files` that will be needed for the sensitivity analysis
#' @param variants dataframe containing all the variants in the cohort. This is created within the `variants_power()` function.
#' @param label character. Any label to be assigned to the current run. This is passed from the argument `down_label` when calling the `variants_power()` function.
#' @param use_transcript logical. Whether transcript information is kept in the sensitivity analysis. This will be passed from `variants_power()` and it will be equivalent to `!TCGA`.

#' @details It is expected that variants were annotated with the Variant Effect Predictor (VEP) https://asia.ensembl.org/info/docs/tools/vep/index.html. This function works with VEP 0.89

extract_fields <- function(variants,label,use_transcript){

  variants <- variants %>%
    mutate(Existing_variation = ifelse(Existing_variation == "",NA,Existing_variation), # set to NA if no Existing_variation provided
           ExAC_AF = as.numeric(as.character(ExAC_AF)),  # convert to numeric
           EUR_AF = as.numeric(as.character(EUR_AF)), # convert to numeric
           SYMBOL = ifelse(SYMBOL == "",NA,SYMBOL),  # set to NA if no Existing_variation provided
           down_label = label,
           variant_type = dplyr::case_when(VARIANT_CLASS %in% "SNV" ~ "SNV", # standardise the way in which we define variant_type
                                    !(VARIANT_CLASS %in% "SNV") ~ "INDEL",
                                    TRUE ~ VARIANT_CLASS)) %>%
    dplyr::filter(!is.na(SYMBOL))                                    # remove variants not on genes

  if(use_transcript){
    variants <- variants %>%
      tidyr::unite("key_SampleName", c("chrom", "pos", "SampleName", "SYMBOL", "Feature"), remove =FALSE,sep=":") # create key
  }else{
    variants <- variants %>%
      tidyr::unite("key_SampleName", c("chrom", "pos", "SampleName", "SYMBOL"), remove =FALSE,sep=":") # create key
  }



  return(unique(variants))

}
