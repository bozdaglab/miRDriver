\name{dummy_methylationData10}
\alias{dummy_methylationData10}
\docType{data}
\title{
dummy_methylationData10}
\description{
A R dataframe object. It contains the gene centric DNA methylation values for each sample. The first column labeled as "sample" contains the sample IDs. The rest of the columns are the genes with methylation values.
}
\format{
A data frame with 300 observations for 10 samples. The first column contains the sample IDs and the rest of the columns are the observations.
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(miRDriverData)
dummy_methylationData10[1:3, 1:3]
}
\keyword{datasets}

