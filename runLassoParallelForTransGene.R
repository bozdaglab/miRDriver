######################################
# File: runLassoParallelForTransGene.R
# Created By: Banabithi Bose
# Date Created: 15th Jan 2019
# Objective: This code will run LASSO using files from the folder 
#            OvarianLASSO/Pipe2_TransGenePredFile.
#            The code will run in parallel in 70 cores.
#            This code will require helper file runLassoParallelForTransGeneHelper.R
#            which is also available in the same location.
######################################

source("runLassoParallelForTransGeneHelper.R")
numCounter = 100 #This counter determines the number of times Lasso is run.
#################
# STEP 1: Make parallel cluster and export the variables to all the clusters
#################
print("cluster making started")
cl <- makeCluster(mc <- getOption("cl.cores", 70),type = "FORK")
clusterEvalQ(cl, library('glmnet'))
clusterEvalQ(cl, library(janitor))
invisible(clusterEvalQ(cl, library(stringr)))
invisible(clusterEvalQ(cl, library(plyr)))
invisible(clusterEvalQ(cl, library(data.table)))
invisible(clusterEvalQ(cl, library(readr)))
invisible(clusterEvalQ(cl, library(gsubfn)))
invisible(clusterEvalQ(cl, library(parallel)))
invisible(clusterEvalQ(cl, library(sqldf)))
invisible(clusterEvalQ(cl, library(reshape)))
invisible(clusterEvalQ(cl, library(reshape2)))
clusterEvalQ(cl, library(pbapply))

print("cluster export started")
parallel::clusterExport(cl=cl,c("numCounter","processLassoCtr"),envir=environment())#envir needed to be correct, export from global aswell as local
print("cluster making done")

#################
# STEP 2: Call the helper file function to run Lasso.
#################
fileNames<-list.files("OvarianLASSO/PipeTransGeneDEPredFile/", pattern="*.Rda") ##we have 10494 GenePred files in this folder
parLapplyLB(cl,fileNames,processUpDownGeneFiles,processLassoCtr)
stopCluster(cl)
print("cluster stopped")

##########End Of Code###############