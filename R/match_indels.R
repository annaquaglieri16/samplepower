
match_indels <- function(truth_set,variants){
  
  if ("pos" %in% colnames(truth_set)){
    start_var <- as.numeric(as.character(truth_set$pos))
    end_var <- as.numeric(as.character(truth_set$pos))
  } else {
    if ( "start" %in% colnames(truth_set) & "end" %in% colnames(truth_set) ) {
      start_var <- as.numeric(as.character(truth_set$start))
      end_var <- as.numeric(as.character(truth_set$end))
    } else {
      stop("truth_set needs to include either a 'pos' column or a 'start' and 'end' column.")
    }
  }
  
  if( !("chrom" %in% colnames(truth_set))){
    stop("truth_set requires a 'chr' column")
  }
  
  if( !("SYMBOL" %in% colnames(truth_set))){
    stop("truth_set requires a 'SYMBOL' column")
  }
  
  if( !("SampleName" %in% colnames(truth_set))){
    stop("truth_set requires a 'SampleName' column")
  }
  
  if( !("Feature" %in% colnames(truth_set))){
    stop("truth_set requires a 'Feature' which indicates the 'ensembl_transcript_id' column.")
  }
  
  if( !("SampleName" %in% colnames(variants_called)) | 
      !("Feature" %in% colnames(variants_called)) | 
      !("SYMBOL" %in% colnames(variants_called))  ){
    stop("variants_called requires: 'SampleName', 'Feature' for ensembl_transcript_id transcript name and 'SYMBOL' columns")
  }
  
  # Create GRanges of truth set
  df <- truth_set %>% dplyr::select(-chrom,-pos)
  truth_set_GR <- GRanges(seqnames = truth_set$chrom,
                          IRanges(start = start_var,
                                  end = end_var),mcols=df)
  
  # as a final result I want to obtain a data frame with ncol >= than ncol(truth_set) that for every variant in the truth set it gives me
  # the match with a variant in the downsampled set and all the information next to it.
  
  # for each variant in the truth set: maybe I should use the foreach function
  match_truth_set <- sapply(1:length(truth_set_GR), function(index_indel){
    
    truth_indel <- truth_set_GR[index_indel]
    
    # 1. Reduce search space
    samplename <- as.character(mcols(truth_indel)$mcols.SampleName)
    gene <- as.character(mcols(truth_indel)$mcols.SYMBOL)
    transcript <- as.character(mcols(truth_indel)$mcols.Feature)
    
    callset_indel <- variants_called %>% filter(SampleName %in% samplename,
                                                SYMBOL %in% gene,
                                                Feature %in% transcript) # Feature :/
    
    df <- callset_indel %>% dplyr::select(-chrom,-pos)
    callset_GR <- GRanges(seqnames = callset_indel$chrom,
                          IRanges(start = as.numeric(as.character(callset_indel$pos)),
                                  end = as.numeric(as.character(callset_indel$pos))),mcols=df)
    
    over <- findOverlaps(query = truth_indel,subject = callset_GR,maxgap = 30)
    
    if(length(over) < 1 ){
      
      df_indel <- data.frame(truth_indel)
      colnames(df_indel) <- str_remove(colnames(df_indel),"mcols.")
      
    } else{
      
      df_indel <- cbind(data.frame(truth_indel),data.frame(mcols(callset_GR[subjectHits(over)])))
      colnames(df_indel) <- str_remove(colnames(df_indel),"mcols.")
      
    }
    
    return(df_indel)
    
  })
  
  df_indel_merge <- do.call(bind_rows,match_truth_set)  %>%
    rename(chrom = seqnames,
           pos = start) %>%
    dplyr::select(-width,-strand)
  return(df_indel_merge)
  
}

