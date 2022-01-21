\name{Dorothea_data}
\alias{Dorothea_data}
\docType{data}
\title{
  Dorothea_data}
\description{
  A R dataframe object that contains gene and transcription factor(TF) interactions with confidence score. The first column labeled as "V1", contains the genes, the second column labeled as "TF", contains the transcription factors and the third column labeled as "Confidence", contains the confidence level of the interaction.}
\format{
  247531 observations of 3 variables
}
\source{
  DoRothEA gene set resource (https://saezlab.github.io/dorothea/)
}
\examples{
  data(miRDriverData)
  dummy_TFData[1:3,]
}
\keyword{datasets}