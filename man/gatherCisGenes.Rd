\name{gatherCisGenes}
\alias{gatherCisGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gathering of Cis genes for each trans gene
}
\description{
This function collects all the differentially expressed cis genes and all miRNA for each trans gene. This function uses the output from the makingCisAndTransGene function as well as a mirna_bedfile and a gistic_bedfile. The gistic_bedfile is described in the makingCisAndTransGene function description and can be made in the getRegionWiseGistic function. The mirna_bedfile is described in the example datasets section and needs to be provided by the user. The output from this function is in the folders gene_all_DEcis and gene_all_mirna. In these folders there is a folder for each trans gene where all the miRNA's or DEcis genes are listed. As with all the functions the mirdirectory needs to be consistent with the other functions.
}
\usage{
gatherCisGenes(ncores=2,mirna_bedfile,gistic_bedfile,
parallel=c('TRUE','FALSE'), mirdirectory="~/")
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
}\item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}\author{
Banabithi Bose
}\examples{
print("This function gathers all of the cis genes that go to the same trans gene.")
##The following code is a valid example of how to use the function but is commented as it
##takes a very long time to run
##load("~/mirDriverFold/EXAMPLE_DATASETS/mirna_bedfile.rda")
##load("~/mirDriverFold/EXAMPLE_DATASETS/gistic_bedfile.rda")
##ncores<- 2
##gatherCisGenes(ncores,mirna_bedfile=mirna_bedfile,gistic_bedfile=FALSE,
##parallel =FALSE, mirdirectory="~/")
}\value{
Returns nothing
}