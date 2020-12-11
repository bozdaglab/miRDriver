\name{getDifferentiallyExpressedGenes}
\alias{getDifferentiallyExpressedGenes}
\title{
Find Differentially Expressed Genes
}
\description{
This function finds the differentially expressed genes based off of peak wise files generated in the getRegionWiseGistic function and RNASeq count data. This function can be run in parallel or on one core. The user sets the number of cores in the ncore argument or uses the default of two. Whether or not the code is run in parallel or not is set with the argument parallel. To run in parallel the argument is set to TRUE, to run linearly the argument is set to FALSE. The argumen RNAcount needs to be set to the file path for the RNASeq coutn data. The peakwise files should be in the correct folder to be used with the code,  "~/mirdirectory/mirDriver/miRDriver_Step1/Regoinwise_Gistic_Files". In the example file path mirdirectory is whatever the user has set the mirdirectory to. If the user has run the getRegionWiseGistic function and set the mirdirectory to the same thing in both functions then the peakwise files are in the correct place.}
\usage{
getDifferentiallyExpressedGenes(ncore=2, RNACount,parallel=c('TRUE','FALSE'), 
mirdirectory="~/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ncore}{
assuming parallel core use, an integer value for the number of cores user wants to use
} \item{RNACount}{
a R dataframe of RNA seq count
} \item{parallel}{
a logic value for whether or not one wants to run the code in parallel
}\item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}\author{
Banabithi Bose
}\examples{
print("This function finds the differentially expressed genes.")
##The following code is a valid example of how to use the function but is commented as it
##takes a very long time to run
##ncore<- 2
##load("~/mirDriverFold/EXAMPLE_DATASETS/RnaSeq_count_data.Rda")
##getDifferentiallyExpressedGenes(ncore, RNACount=RnaSeq_count_data,
##parallel = TRUE, mirdirectory="~/")

}\value{
Returns nothing
}
