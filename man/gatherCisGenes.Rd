\name{gatherCisGenes}
\alias{gatherCisGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gathering of Cis genes for each trans gene
}
\description{
This function gathers all of the cis genes that go to the same trans gene.
}
\usage{
gatherCisGenes(ncores=2,mirna_bedfile,gistic_bedfile,parallel=c('TRUE','FALSE'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncores}{
  Assuming parallel, and integer value that represents the number of cores that the user wants to use
} \item{mirna_bedfile}{
  A dataframe in the format of a bedfile for all of the micro RNA
} \item{gistic_bedfile}{
  A dataframe in the format of a bedfile for the gistic results
} \item{parallel}{
  A logic value for whether or not one wants to run the code in parallel
}
}\author{
Banabithi Bose
}\examples{
load("~/mirDriverFold/EXAMPLE_DATASETS/mirna_bedfile.Rda")
load("~/mirDriverFold/EXAMPLE_DATASETS/gistic_bedfile.rda")
ncores<- 2
  \donttest{
    gatherCisGenes(ncores,mirna_bedfile=mirna_bedfile,gistic_bedfile=FALSE,parallel        =TRUE)
  }
}\value{
Returns nothing
}
