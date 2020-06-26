\name{getRegionWiseGistic2}
\alias{getRegionWiseGistic2}
\title{
produce regionwise Gistic files i.e. peak wise files
}
\description{
getRegionWiseGistic2 will produce regionwise Gistic files i.e. aberration peak wise files
}
\usage{
getRegionWiseGistic2(gisticfile,gistic_bedfile)
}
\arguments{
  \item{gisticfile}{
  a character string for the memory location of the gistic file or amplification/deletion peaks with patients file in .txt
} \item{gistic_bedfile}{
  a logical value for whether or not one wants a gistic_bedfile
}
}\author{
Banabithi Bose
}\examples{
gisticfile<-"~/mirDriverFold/EXAMPLE_DATASETS/all_lesions.conf_90.txt"
  \donttest{
    getRegionWiseGistic2(gisticfile=gisticfile,gistic_bedfile=TRUE)
  }
}\value{
Returns nothing
}
