###########################################################
# Extract only fields that I am interested in at the moment
# Byotype is not unique so I would have multiple hits for the same variant]
# Variants should have been VEP annotated
##########################################################

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

  if(!use_transcript){
    variants <- variants %>%
      tidyr::unite("key_SampleName", c("chrom", "pos", "SampleName", "SYMBOL", "Feature"), remove =FALSE,sep=":") # create key
  }else{
    variants <- variants %>%
      tidyr::unite("key_SampleName", c("chrom", "pos", "SampleName", "SYMBOL"), remove =FALSE,sep=":") # create key
  }



  return(unique(variants))

}
