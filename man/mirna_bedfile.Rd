\name{mirna_bedfile}
\alias{mirna_bedfile}
\docType{data}
\title{
x}
\description{
This is an R dataframe object where the data represents the genomic position of the miRNAs. The rows in this dataset are the chromosomal regions that start with "chr" followed by chromosome numbers as integers 1-22.. The columns are "chromosome number", "start position", "end position", and "miRNA ID".}
\format{
2024 observations of 4 variables
}
\source{
TCGA Database (https://portal.gdc.cancer.gov/)
}
\examples{
data(miRDriverData)
## maybe str(miRDriverData) ; plot(miRDriverData) ...
}
\keyword{datasets}
