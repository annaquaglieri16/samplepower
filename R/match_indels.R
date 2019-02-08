
match_indels <- function(truth_set,variants,use_transcript,maxgap=50){

  # If an INDELS is reported in the truth set with a start != end then I will use this range to find an overlap?
  if ( "pos" %in% colnames(truth_set) & "end" %in% colnames(truth_set) ) {
    start_var <- as.numeric(as.character(truth_set$pos))
    end_var <- as.numeric(as.character(truth_set$end))
  } else {
    if ("pos" %in% colnames(truth_set)){
      start_var <- as.numeric(as.character(truth_set$pos))
      end_var <- as.numeric(as.character(truth_set$pos))
    } else {
      stop("truth_set needs to include either a 'pos' column or a 'pos' and 'end' columns.")
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

  if(use_transcript){
    if( !("Feature" %in% colnames(truth_set))){
      stop("truth_set requires a 'Feature' which indicates the 'ensembl_transcript_id' column.")
    }
  }

  if(use_transcript){

    need_colmuns <- c("SampleName","Feature","SYMBOL")
    check_columns <- sum(!(need_colmuns %in% colnames(variants)))

    if(check_columns > 0){
      missing <- need_colmuns[!(need_colmuns %in% colnames(variants))]
      stop(paste0("Check requirements for column names of 'variants' The following columns are missing: ",missing))
    }
  }else{
    need_colmuns <- c("SampleName","SYMBOL")
    check_columns <- sum(!(need_colmuns %in% colnames(variants)))

    if(check_columns > 0){
      missing <- need_colmuns[!(need_colmuns %in% colnames(variants))]
      stop(paste0("Check requirements for column names of 'variants' The following columns are missing: ",missing))
    }
  }

  # Create GRanges of truth set
  df <- truth_set %>% dplyr::select(-chrom,-pos) # locus is left to indicate the position in the truth set
  truth_set_GR <- GRanges(seqnames = truth_set$chrom,
                          IRanges(start = start_var,
                                  end = end_var),mcols=df)

  # as a final result I want to obtain a data frame with ncol >= than ncol(truth_set) that for every variant in the truth set it gives me
  # the match with a variant in the downsampled set and all the information next to it.

  # for each variant in the truth set
  match_truth_set <- lapply(1:length(truth_set_GR), function(index_indel){

    # 1. Select one indel
    truth_indel <- truth_set_GR[index_indel]

    # 2. Reduce search space using only information in truth set
    samplename <- as.character(mcols(truth_indel)$mcols.SampleName)
    gene <- as.character(mcols(truth_indel)$mcols.SYMBOL)

    if(use_transcript){

      transcript <- as.character(mcols(truth_indel)$mcols.Feature)

      callset_indel <- variants %>% filter(SampleName %in% samplename,
                                                SYMBOL %in% gene,
                                                Feature %in% transcript)
    } else {

      callset_indel <- variants %>% filter(SampleName %in% samplename,
                                           SYMBOL %in% gene)

    }

    # The caller reports only the pos of an INDEL
    df <- callset_indel %>% dplyr::select(-chrom,-pos)
    callset_GR <- GRanges(seqnames = callset_indel$chrom,
                            IRanges(start = as.numeric(as.character(callset_indel$pos)),
                                    end = as.numeric(as.character(callset_indel$pos))),mcols=df)

    over <- findOverlaps(query = truth_indel,subject = callset_GR,maxgap = maxgap)

    if(length(over) < 1 ){

      df_indel <- data.frame(truth_indel)
      colnames(df_indel) <- stringr::str_remove(colnames(df_indel),"mcols.")

    } else{

      df_indel <- cbind(data.frame(truth_indel),data.frame(mcols(callset_GR[subjectHits(over)]))) # I bind the positionfrom the truth set with the mcols from the callset, Location from the variants set should preserve the position
      # found by a caller
      colnames(df_indel) <- stringr::str_remove(colnames(df_indel),"mcols.")

    }

    return(df_indel)

  })

  df_indel_merge <- do.call(bind_rows,match_truth_set)  %>%
    dplyr::rename(chrom = seqnames,
           pos = start) %>%
    dplyr::select(-width,-strand)
  return(df_indel_merge)

}

