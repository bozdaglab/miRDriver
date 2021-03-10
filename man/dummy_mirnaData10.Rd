\name{dummy_mirnaData10}
\alias{dummy_mirnaData10}
\docType{data}
\title{
dummy_mirnaData10}
\description{
A R dataframe object. It contains the miRNA expression data for the samples. The first column is labeled as "sample" and contains the sample IDs. The rest of the columns are the miRNAs with expression values.
}
\format{
A data frame with 300 observations for 10 samples. The first column contains the sample IDs and the rest of the columns are the observations.
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(miRDriverData)
dummy_mirnaData10[1:3, 1:3]
}
\keyword{datasets}
