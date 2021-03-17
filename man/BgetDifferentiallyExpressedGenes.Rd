\name{getDifferentiallyExpressedGenes}
\alias{getDifferentiallyExpressedGenes}
\title{
Find Differentially Expressed Genes
}
\description{
This function finds the differentially expressed genes based off of the peak wise files generated in the getRegionWiseGistic function and RNASeq count data. This function uses parallel computation by setting the argument value of *ncore*; the default is *ncore = 2*. The argument *RNAcount* needs to be set to the file path for the RNASeq count data. The argument *mirdirectory* needs to be the same as in all previous functions.
}\arguments{
\item{ncore}{
An integer value that represents the number of cores that the user wants to use (default = 2)
} \item{RNACount}{
a R dataframe of RNA seq count
}\item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}\author{
Banabithi Bose
}\examples{


load(system.file("data", "miRDriverData.RData", package = "miRDriver"))

getDifferentiallyExpressedGenes(ncore = 1, RNACount = Dummy_RnaSeq_Count_Data,
                                mirdirectory = tempdir())
}\value{
Returns nothing
}
