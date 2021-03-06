% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_fields.R
\name{extract_fields}
\alias{extract_fields}
\title{Extract specific fields from `variant_files` that will be needed for the sensitivity analysis}
\usage{
extract_fields(variants, label, use_transcript)
}
\arguments{
\item{variants}{dataframe containing all the variants in the cohort. This is created within the `variants_power()` function.}

\item{label}{character. Any label to be assigned to the current run. This is passed from the argument `down_label` when calling the `variants_power()` function.}

\item{use_transcript}{logical. Whether transcript information is kept in the sensitivity analysis. This will be passed from `variants_power()` and it will be equivalent to `!TCGA`.}
}
\description{
Extract specific fields from `variant_files` that will be needed for the sensitivity analysis
}
\details{
It is expected that variants were annotated with the Variant Effect Predictor (VEP) https://asia.ensembl.org/info/docs/tools/vep/index.html. This function works with VEP 0.89
}
