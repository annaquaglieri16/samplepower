#' Parse a VCF input to return a data frame with standardised fields across callers to use for caller comparison. It only works with germline calls and for VCF from the following callers: GATK3 MuTect2, VarScan2 and VarDict.
#' @param vcf_path path to where the `.vcf` file for one sample is saved.
#' @param sample_name character. Sample name of the current `vcf` file.
#' @param caller character. One of `mutect`, `vardict` or `varscan`.


parse_vcf_output <- function(vcf_path, sample_name, caller) {

  vcf <- VariantAnnotation::readVcf(vcf_path)

  # Check if is comes from somatic calls
  if(ncol(VariantAnnotation::geno(vcf)$GT) > 1){

    stop("parse_vcf_output wasn't implemented for somatic calls")

  } else {

    if(caller == "varscan"){

      vcf_df <- data.frame(data.frame(ranges(vcf)),
                           genotype= VariantAnnotation::geno(vcf)$GT[,1],
                           #SampleName = VariantAnnotation::info(vcf)$SAMPLE,
                           #qual = VariantAnnotation::qual(vcf),
                           filter = VariantAnnotation::filt(vcf),
                           ref_base_quality = VariantAnnotation::geno(vcf)$RBQ[,1],
                           alt_base_quality = VariantAnnotation::geno(vcf)$ABQ[,1],
                           tot_depth = VariantAnnotation::geno(vcf)$DP[,1],
                           freq = VariantAnnotation::geno(vcf)$FREQ[,1],
                           ref_depth = VariantAnnotation::geno(vcf)$RD[,1],
                           alt_depth = VariantAnnotation::geno(vcf)$AD[,1],
                           ref_forw = VariantAnnotation::geno(vcf)$RDF[,1],
                           ref_rev = VariantAnnotation::geno(vcf)$RDR[,1],
                           alt_forw = VariantAnnotation::geno(vcf)$ADF[,1],
                           alt_rev = VariantAnnotation::geno(vcf)$ADR[,1]) %>%

        tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
        tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
        tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
        dplyr::mutate(caller="varscan",
                      Location = gsub(":","_",Location),
                      SampleName = sample_name) %>%
        dplyr::mutate(VAF = parse_vaf_varscan(freq),
                      qual = (ref_base_quality + alt_base_quality)/2) %>%  # mean of alt/ref base qualitites
        dplyr::select(-freq)

    }

    if(caller == "mutect"){

      vcf_df <- data.frame(data.frame(ranges(vcf)),
                           genotype= VariantAnnotation::geno(vcf)$GT[,1],
                           #SampleName = VariantAnnotation::info(vcf)$SAMPLE,
                           #qual = VariantAnnotation::qual(vcf),
                           filter = VariantAnnotation::filt(vcf),
                           base_quality = VariantAnnotation::geno(vcf)$QSS[,1],
                           alleles_depth = VariantAnnotation::geno(vcf)$AD[,1],
                           VAF = VariantAnnotation::geno(vcf)$AF[,1],
                           alt_forw = VariantAnnotation::geno(vcf)$ALT_F1R2[,1],
                           alt_rev = VariantAnnotation::geno(vcf)$ALT_F2R1[,1],
                           ref_forw = VariantAnnotation::geno(vcf)$REF_F1R2[,1],
                           ref_rev = VariantAnnotation::geno(vcf)$REF_F2R1[,1]) %>%

        tidyr::separate(alleles_depth, into = c("ref_depth","alt_depth"),sep=" ",remove=TRUE) %>%
        dplyr::mutate(ref_depth =  as.numeric(as.character(ref_depth)),
                      alt_depth = as.numeric(as.character(alt_depth)),
                      tot_depth = ref_depth + alt_depth) %>%

        tidyr::separate(base_quality,into = c("ref_base_quality","alt_base_quality"),sep = ",",remove=TRUE) %>%
        dplyr::mutate(ref_base_quality =  as.numeric(as.character(ref_base_quality)),
                      alt_base_quality = as.numeric(as.character(alt_base_quality))) %>%
        dplyr::mutate(qual_ref = ifelse(ref_depth == 0,0,ref_base_quality/ref_depth),
                      qual_alt= ifelse(alt_depth == 0,0,alt_base_quality/alt_depth)) %>%

        dplyr::mutate(qual = qual_ref + qual_alt/2) %>%
        dplyr::select(-qual_ref,-qual_alt) %>%

        tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
        tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
        tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
        dplyr::mutate(caller="mutect2",
                      Location = gsub(":","_",Location))

    }

    if(caller == "vardict"){

      vcf_df <- data.frame(data.frame(ranges(vcf)),
                           genotype= VariantAnnotation::geno(vcf)$GT[,1],
                           #SampleName = VariantAnnotation::info(vcf)$SAMPLE,
                           qual = VariantAnnotation::info(vcf)$QUAL,
                           filter = VariantAnnotation::filt(vcf),
                           DP = VariantAnnotation::info(vcf)$DP,
                           VAF = VariantAnnotation::info(vcf)$AF,
                           ADJVAF_ADJ_indels = VariantAnnotation::info(vcf)$ADJAF,
                           VD = VariantAnnotation::info(vcf)$VD,
                           REFBIAS = VariantAnnotation::info(vcf)$REFBIAS,
                           VARBIAS = VariantAnnotation::info(vcf)$VARBIAS) %>%

        tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
        tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
        tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
        tidyr::separate(REFBIAS,into=c("ref_forw","ref_rev"),sep=":") %>%
        tidyr::separate(VARBIAS,into=c("alt_forw","alt_rev"),sep=":") %>%
        dplyr::mutate(ref_depth = DP - VD) %>%
        dplyr::rename(tot_depth = DP,
                      alt_depth = VD) %>%
        dplyr::select(-start) %>%
        dplyr::mutate(caller="vardict",
                      Location = gsub(":","_",Location))

    }

  }

    vcf_df <- vcf_df %>% dplyr::select(Location,caller,chrom,pos,end,ref,alt,qual,filter,
                                        genotype,tot_depth,VAF,ref_depth,
                                        alt_depth,ref_forw,ref_rev,alt_forw,alt_rev,everything())


    return(vcf_df)

}
