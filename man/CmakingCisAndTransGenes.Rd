\name{makingCisAndTransGenes}
\alias{makingCisAndTransGenes}
\title{
Making Cis and Trans Genes
}
\description{
This function finds the trans and cis genes for each peak. The trans and cis genes will be stored in R dataframe objects located in the file path ~/mirdirectory/miRDriver_Step2/Trans_Cis_Files/peak/. In the file path the mirdirectory is the mirdirectory specified in the mirdirectory argument and should be the same as the mirdirectory used in the previous functions. The peak in the file path represents the sundirectories for each peak which will follow the convention "Amplification Peak 1". The gistic_bedfile argument needs to be set to a R dataframe object that follows the convention of having four columns, "chromosome", "start position", "end position", and "peak". The chromosome column will have the chromosome number in the form "chr1". The start and end positions will be integers representing the start and end positions on the chromosome. The peak column will have a string representing the peak following the convention of "Amplification Peak 1". The gistic_bedfile can be made in the getRegionWiseGisticFunction.}
\arguments{
    \item{Genes_coord_bedfile}{
A dataframe that fits the format of a bedfile for gene coordinates
} \item{gistic_bedfile}{
A dataframe that fits the format of a bedfile for the gistic results with amplification/deletion peaks with patients
}\item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}
\author{
Banabithi Bose
}
\examples{

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

makingCisAndTransGenes(
  Genes_coord_bedfile = Genes_coord,
  gistic_bedfile = region,
  mirdirectory = tempdir()
)

}\value{
Returns nothing
}