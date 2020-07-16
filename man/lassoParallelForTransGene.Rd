\name{lassoParallelForTransGene}
\alias{lassoParallelForTransGene}
\title{
Lasso Regression For All Predictors
}
\description{
This function runs a lasso regression on all of the predictors for a data set.
}
\usage{
lassoParallelForTransGene(ncores = 2,numCounter = 100,Nfolds = 10,
nonZeroPercent = 70 ,parallel=c('TRUE','FALSE'))
}
\arguments{
    \item{ncores}{
assuming the code in run in parallel, an integer for the number of cores one wants to use
} \item{numCounter}{
an integer for the number of times the lasso regression is run
} \item{Nfolds}{
an integer for the number of folds one wants for their lasso regression
} \item{nonZeroPercent}{
an integer for the percentage of appearance as a non-zero LASSO coefficient for a mirna to be selected
} \item{parallel}{
a logic value for whether or not the user wants to run the code in parallel, if the code is run on windows, must be false
}
}\author{
Banabithi Bose
}\examples{
print("This function runs a lasso regression on all of the predictors for a data set.")
##The following code is a valid example of how to use the function but is commented as it
##takes a very long time to run
##ncores=2
##Nfolds<-3
##numCounter = 20
##nonZeroPercent = 10
##x<-lassoParallelForTransGene(ncores = 70,numCounter = 3,Nfolds = 3,nonZeroPercent = 10     
## ,parallel = TRUE)

}\value{
Returns nothing
}