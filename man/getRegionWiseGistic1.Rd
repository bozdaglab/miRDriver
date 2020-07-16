\name{getRegionWiseGistic1}
\alias{getRegionWiseGistic1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
run GISTIC 2.0 and produce regionwise Gistic files i.e. aberration peak wise files
}
\description{
getRegionWiseGistic1 will run GISTIC 2.0 and will produce regionwise Gistic files i.e. peak wise files
}
\usage{
getRegionWiseGistic1(humanGene,copyNumber,conf.lev,username,
                    password,gistic_bedfile=c('TRUE','FALSE'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{humanGene}{
a character string that is a link to a matlab file for the humanGene link of choice
} \item{copyNumber}{
a character string that represents a link to a copy number text file stored in user's Genepattern account
} \item{conf.lev}{
a double value that represents the confidence level of GISTIC 2.0
} \item{username}{
a character string that is the username to user's Genepattern account
} \item{password}{
a character string that is the password to user's Genepattern account
} \item{gistic_bedfile}{
a logic value of  either true or false representing whether or not one wants to produce a bedfile from GISTIC 2.0 output
}
}\author{
Banabithi Bose
}\examples{
print("getRegionWiseGistic1 will run GISTIC 2.0 and will produce
regionwise Gistic files i.e. peak wise files")
##The following code is a valid example of how to use the function but is commented as it
##takes a very long time to run
##humanGene<- "Human_Hg38.UCSC.add_miR.160920.refgene.mat"
##copyNumber<- "https://cloud.genepattern.org/gp/users/miRDriver/ovary_brca_cnv_masked.txt"
##conf.lev<-0.90
##username <- "miRDriver"
##password <- "1234"
##getRegionWiseGistic1(humanGene,copyNumber,conf.lev,username,password,gistic_bedfile=TRUE)
}\value{
Returns nothing
}