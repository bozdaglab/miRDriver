\name{gatherCisGenes}
\alias{gatherCisGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gathering of Cis genes for each trans gene
}
\description{
This function collects all the DE cis genes and all miRNAs for each trans gene. This function uses the output from the *makingCisAndTransGene()* function as well as a *mirna_bedfile* and a *gistic_bedfile*. The gistic_bedfile is described in the makingCisAndTransGene function description and can be made using the *getRegionWiseGistic()* function. The *mirna_bedfile* is described in the example datasets section and needs to be provided by the user. The output from this function is in the folders gene_all_DEcis and gene_all_mirna. In these folders there is a folder for each trans gene where all the miRNA's or DE cis genes are listed. As with all the functions the mirdirectory needs to be consistent with the other functions.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ncore}{
An integer value that represents the number of cores that the user wants to use (default = 2)
} \item{mirna_bedfile}{
A dataframe in the format of a bedfile for all of the micro RNA
} \item{gistic_bedfile}{
A dataframe in the format of a bedfile for the gistic results
} \item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}\author{
Banabithi Bose
}\examples{

load(
  file.path(
    tempdir(),
    "mirDriverFold",
    "miRDriver_Step1",
    "GisticResults",
    "gistic_bedfile.rda"
  )
)

load(system.file("data", "miRDriverData.RData", package = "miRDriver"))

gatherCisGenes(
  ncore = 1,
  mirna_bedfile = mirna_bedfile,
  gistic_bedfile = region,
  mirdirectory = tempdir()
)
}\value{
Returns nothing
}
