\name{dummy_TFData}
\alias{dummy_TFData}
\docType{data}
\title{
dummy_TFData}
\description{
A R dataframe object that contains gene and transcription factor(TF) interactions with confidence score. The first column labeled as "V1", contains the genes, the second column labeled as "TF", contains the transcription factors and the third column labeled as "Confidence", contains the confidence level of the interaction.}
\format{
653 observations of 3 variables
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(miRDriverData)
dummy_TFData[1:3,]
}
\keyword{datasets}