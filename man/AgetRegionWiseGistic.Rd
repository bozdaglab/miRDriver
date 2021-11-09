\name{getRegionWiseGistic}
\alias{getRegionWiseGistic}
\title{
produce regionwise Gistic files i.e. peak wise files
}
\description{
This function used the lesions file to find the aberrations for each peak. The lesions file, as stated in the example data sets section, can be created with the GISTIC2.0 module on genepattern.com (https://cloud.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00125:6.15.28).If one uses GISTIC2.0 to make the lesions file, then the lesions file will be labeled all_lesions.conf_90.txt. As one can see in the example code, this file is used in the function as the argument gisticfile. The gistic_bedfile argument should be either TRUE or FALSE, representing if the user would like the function to make a gistic_bedfile. The mirdirectory argument needs to be the file path that the user wants the results to be in. The mirdirectory needs to be the same for all functions in the same analysis.}
\arguments{
    \item{gisticfile}{
a character string for the memory location of the gistic file or amplification/deletion peaks with patients file in .txt
} \item{gistic_bedfile}{
a logical value for whether or not one wants a gistic_bedfile
}\item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}\author{
Banabithi Bose
}\examples{

gisticfile <-
  system.file("extdata", "Dummy.gistic.lesions.txt", package = "miRDriver")

getRegionWiseGistic(
  gisticfile = gisticfile,
  gistic_bedfile = TRUE,
  mirdirectory = tempdir()
)
}\value{
Returns nothing
}