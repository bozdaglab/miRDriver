\name{getTransGenePredictorFile}
\alias{getTransGenePredictorFile}
\title{
Creation of the trans gene predictor file
}
\description{
This function creates a TransGenePredictorFile for each trans gene based on the mirna and DEcis files from the output of the gatherCisGenes function and the user provided dataframes methylationData, RSeqData, CNVData, TFData, mirnaData. The five dataframes needed for this function are described in the example datasets section. The TransGenePredictorFiles will contain information about the possible predictors for each patient based on the data provided. The TransGenePredictorFiles will be saved in the file path ~/mirdirectory/mirDriverFold;miRDriver_Step4/TransGenePredFile/. As always mirdirectory in the file path represents the mirdirectory chosen for the analysis and should be the same as the mirdirectory used in all previous functions.
}
\usage{
getTransGenePredictorFile(ncores=2, methylationData, RSeqData, CNVData, 
TFData,mirnaData,parallel=c('TRUE','FALSE'), mirdirectory = "~/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ncores}{
assuming parallel use, the number of cores the user wants to use
} \item{methylationData}{
a dataframe for the methylation data
} \item{RSeqData}{
a dataframe for the RNA sequence data
} \item{CNVData}{
a dataframe for the copy number varation data
} \item{TFData}{
a dataframe for the transcription factor data
} \item{mirnaData}{
a dataframe for the mirna data
} \item{parallel}{
a logic value for whether or not the user wants to run the code in parallel
}\item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}
\author{
Banabithi Bose
}
\examples{
print("This function creates the trans gene predictor file.")
##The following code is a valid example of how to use the function but is commented as it
##takes a very long time to run
##ncores<- 2
##load("~/mirDriverFold/EXAMPLE_DATASETS/methylationData.Rda")
##load("~/mirDriverFold/EXAMPLE_DATASETS/RSeqData.Rda")
##load("~/mirDriverFold/EXAMPLE_DATASETS/CNVData.Rda")
##load("~/mirDriverFold/EXAMPLE_DATASETS/TFData.Rda")
##load("~/mirDriverFold/EXAMPLE_DATASETS/mirnaData.Rda")
##getTransGenePredictorFile(ncores, methylationData, RSeqData, CNVData,TFData,mirnaData      
## ,(parallel = TRUE), mirdirectory="~/")
}\value{
Returns nothing
}