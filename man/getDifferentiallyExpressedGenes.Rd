\name{getDifferentiallyExpressedGenes}
\alias{getDifferentiallyExpressedGenes}
\title{
Find Differentially Expressed Genes
}
\description{
This function finds the differentially expressed genes.
}
\usage{
getDifferentiallyExpressedGenes(ncore=2, RNACount,parallel=c('TRUE','FALSE'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncore}{
  assuming parallel core use, an integer value for the number of cores user wants to use
} \item{RNACount}{
  a R dataframe of RNA seq count
} \item{parallel}{
  a logic value for whether or not one wants to run the code in parallel
}
}\author{
Banabithi Bose
}\examples{
ncore<- 2
load("~/mirDriverFold/EXAMPLE_DATASETS/RnaSeq_count_data.Rda")
  \donttest{
    getDifferentiallyExpressedGenes(ncore, RNACount=RnaSeq_count_data,parallel = TRUE)
  }
}\value{
Returns nothing
}

