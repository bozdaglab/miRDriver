\name{lassoParallelForTransGene}
\alias{lassoParallelForTransGene}
\title{
Lasso Regression For All Predictors
}
\description{
This function runs a lasso regression to find what contributes to the expression of the trans gene. This type of regression eliminate non-contributors by settting their coefficient to zero. If the coeffiecient of a miRNA is non-zero after the regression finishes then it is a contributor. However, to be accurate the regression is set to run multiple times (the user specifies as an argument), with a certain percentage of times needed with the coefficient to be non-zero. (the user specifies this as an argument as well). The user also needs to provide the number of folds they want their regression to use. The regression is run on the TransGenePredictorFiles created from the getTransPredictorFile function. The output from this function is found in the file path ~/mirdirectory/mirDriverFold/miRDriver_RESULTS/. The mirdirectory in the file path represents the directory chosen by the user and specified in the argument, mirdirectory. The argument mirdirectory needs to be set to the same directory as the mirdirectory in the previous steps. }
\arguments{
    \item{ncore}{
An integer value that represents the number of cores that the user wants to use (default = 2)
} \item{numCounter}{
an integer for the number of times the lasso regression is run
} \item{Nfolds}{
an integer for the number of folds one wants for their lasso regression
} \item{nonZeroPercent}{
an integer for the percentage of appearance as a non-zero LASSO coefficient for a mirna to be selected
}\item{mirdirectory}{
a character string for the file path of the directory one wants the results in}
}\author{
Banabithi Bose
}\examples{
lassoParallelForTransGene(
  ncore = 1,
  numCounter = 10,
  Nfolds = 3,
  nonZeroPercent = 20,
  mirdirectory = tempdir()
)
}\value{
Returns nothing
}