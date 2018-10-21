############
## Parse PON
############

create_vaf <- function(vaf){
  as.numeric(substr(vaf,1,nchar(vaf)-1))/100
}

# Extract VAF from VarScan output
parse_vaf_varscan <- function(freq){
  freq <- gsub("%","",freq)
  freq <- as.numeric(as.character(freq))/100
  return(freq)
}

vaf_to_number <- function(freq){
  freq <- as.numeric(as.character(freq))
  return(freq)
}

# Nset
n_normals <- function(set){
  length(gregexpr("-", set)[[1]])+1
}

# path_to_pon is the link to the file created with combinevariants using GATK

parse_pon <- function(path_to_pon,
                      caller){



  pon <- fread(path_to_pon,header = TRUE)
  pon <- data.frame(pon)


  # Extract VAF for every variant called by the caller
  # dtermine number of normal samples where eah variant was called
  # determine the minimum VAF that each variant was called with

  if(caller %in% "varscan"){

    FREQcol <- grep("FREQ",colnames(pon))
    pon_FREQ <- pon[,FREQcol]
    pon_VAF <- apply(pon_FREQ,2,function(x) parse_vaf_varscan(x))
    minVAF <- apply(pon_VAF,1,min,na.rm=TRUE)
    combine_vaf <- cbind(pon[,c("CHROM","POS","REF","ALT","set")],minVAF)
  }

  if(caller %in% c("mutect","vardict")){
    FREQcol <- grep(".AF",colnames(pon))
    pon_FREQ <- pon[,FREQcol]
    pon_VAF <- apply(pon_FREQ,2,function(x) vaf_to_number(x))
    minVAF <- apply(pon_VAF,1,min,na.rm=TRUE)
    combine_vaf <- cbind(pon[,c("CHROM","POS","REF","ALT","set")],minVAF)
  }

  combine_vaf$nsam <- sapply(combine_vaf$set,function(entry_set) n_normals(entry_set))

  return(combine_vaf)
}
