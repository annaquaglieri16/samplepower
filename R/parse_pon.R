############
## Parse PON
############

create_vaf <- function(vaf){
  as.numeric(substr(vaf,1,nchar(vaf)-1))/100
}

# Extract VAF from VarScan output
VAFvarscan <- function(freq){
  freq <- gsub("%","",freq)
  freq <- as.numeric(as.character(freq))/100
  return(freq)
}

VAFother <- function(freq){
  freq <- as.numeric(as.character(freq))
  return(freq)
}


parsePON <- function(linkPON,caller="vardict"){
  pon <- fread(linkPON,header = TRUE)
  pon <- data.frame(pon)


  # Extract VAF for every variant called by the caller
  # dtermine number of normal samples where eah variant was called
  # determine the minimum VAF that each variant was called with

  if(caller %in% "varscan"){
    FREQcol <- grep("FREQ",colnames(pon))
    pon_FREQ <- pon[,FREQcol]
    pon_VAF <- apply(pon_FREQ,2,function(x) VAFvarscan(x))
    minVAF <- apply(pon_VAF,1,min,na.rm=TRUE)
    combine_vaf <- cbind(pon[,c("CHROM","POS","REF","ALT","set")],minVAF)
  }

  if(caller %in% c("mutect","vardict")){
    FREQcol <- grep(".AF",colnames(pon))
    pon_FREQ <- pon[,FREQcol]
    pon_VAF <- apply(pon_FREQ,2,function(x) VAFother(x))
    minVAF <- apply(pon_VAF,1,min,na.rm=TRUE)
    combine_vaf <- cbind(pon[,c("CHROM","POS","REF","ALT","set")],minVAF)
  }

  # Nset
  Nsam <- function(set){
    length(gregexpr("SRX", set)[[1]])
  }
  combine_vaf$nsam <- sapply(combine_vaf$set,function(entry_set) Nsam(entry_set))
  combine_vaf$nsam_prop <- combine_vaf$nsam/17

  return(combine_vaf)
}

