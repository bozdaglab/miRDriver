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
load("~/mirDriverFold/EXAMPLE_DATASETS/hg38_genes_coord_bedfile.rda")
load("~/mirDriverFold/EXAMPLE_DATASETS/gistic_bedfile.rda")
  \donttest{
    makingCisAndTransGenes(Genes_coord_bedfile=Genes_coord, gistic_bedfile=region)
  }
}\value{
Returns nothing
}
