\name{makingCisAndTransGenes}
\alias{makingCisAndTransGenes}
\title{
Making Cis and Trans Genes
}
\description{
This function makes the cis and trans genes for the differentially expressed genes
}
\usage{
makingCisAndTransGenes(Genes_coord_bedfile,gistic_bedfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Genes_coord_bedfile}{
A dataframe that fits the format of a bedfile for gene coordinates
} \item{gistic_bedfile}{
A dataframe that fits the format of a bedfile for the gistic results with amplification/deletion peaks with patients
}
}
\author{
Banabithi Bose
}
\examples{
print("This function makes the cis and trans genes for the differentially expressed genes")
##The following code is a valid example of how to use the function but is commented as it
##takes a very long time to run 
##load("~/mirDriverFold/EXAMPLE_DATASETS/hg38_genes_coord_bedfile.rda")
##load("~/mirDriverFold/EXAMPLE_DATASETS/gistic_bedfile.rda")
##makingCisAndTransGenes(Genes_coord_bedfile=Genes_coord, gistic_bedfile=region)

}\value{
Returns nothing
}