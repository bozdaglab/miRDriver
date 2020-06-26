\name{getTransGenePredictorFile}
\alias{getTransGenePredictorFile}
\title{
Creation of the trans gene predictor file
}
\description{
This function creates the trans gene predictor file which is a dataframe with all of the contributors to the expression including methylation data, RNA sequence data, copy number varaition data, transcription factor data and mirna data.
}
\usage{
getTransGenePredictorFile(ncores=2, methylationData, RSeqData, CNVData, TFData,
                         mirnaData,parallel=c('TRUE','FALSE'))
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
}
}
\author{
Banabithi Bose
}
\examples{
ncores<- 2
load("~/mirDriverFold/EXAMPLE_DATASETS/methylationData.Rda")
load("~/mirDriverFold/EXAMPLE_DATASETS/RSeqData.Rda")
load("~/mirDriverFold/EXAMPLE_DATASETS/CNVData.Rda")
load("~/mirDriverFold/EXAMPLE_DATASETS/TFData.Rda")
load("~/mirDriverFold/EXAMPLE_DATASETS/mirnaData.Rda")
  \donttest{
    getTransGenePredictorFile(ncores, methylationData, RSeqData, CNVData,
                          TFData,mirnaData,(parallel = TRUE))
  }
}\value{
Returns nothing
}
