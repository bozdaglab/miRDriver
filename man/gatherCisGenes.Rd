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
A logic value for whether or not one wants to run the code in parallel, if running code on         windows, must be false
}
}\author{
Banabithi Bose
}\examples{
print("This function gathers all of the cis genes that go to the same trans gene.")
##The following code is a valid example of how to use the function but is commented as it
##takes a very long time to run
##load("~/mirDriverFold/EXAMPLE_DATASETS/mirna_bedfile.rda")
##load("~/mirDriverFold/EXAMPLE_DATASETS/gistic_bedfile.rda")
##ncores<- 2
##gatherCisGenes(ncores,mirna_bedfile=mirna_bedfile,gistic_bedfile=FALSE,parallel =FALSE)
}\value{
Returns nothing
}