######################################
# File: createTransGenePredictorFile.R
# Created By: Banabithi Bose
# Date Created: 15th Jan 2019
# Objective: This code will create the Trans Gene Predictor Files.
#            The code will run in parallel in 70 cores.
#            This code will require helper file createTransGenePredictorFileHelper.R
#            which is also available in the same location.
######################################

source("createTransGenePredictorFileHelper.R")
#################
# STEP 1: Load all the necessary dataframes
#################
load("OvarianLASSO/DATAFRAMES/OV_All_methylation.Rda")
load("OvarianLASSO/DATAFRAMES/RSeq_RPKM.Rda")
load("OvarianLASSO/DATAFRAMES/OV_LASSO_CNV.Rda")
load("OvarianLASSO/DATAFRAMES/TF_in_Rseq.Rda")
load("OvarianLASSO/DATAFRAMES/OV_filtered_mirna_0.01_30_gistic.Rda")
#################


#################
# STEP 2: Make parallel cluster and export the dataframes to all the clusters
#################
print("cluster making started")
cl <- makeCluster(mc <- getOption("cl.cores", 70))
invisible(clusterEvalQ(cl, library('glmnet')))
invisible(clusterEvalQ(cl, library(janitor)))
invisible(clusterEvalQ(cl, library(stringr)))
invisible(clusterEvalQ(cl, library(plyr)))
invisible(clusterEvalQ(cl, library(data.table)))
invisible(clusterEvalQ(cl, library(readr)))
invisible(clusterEvalQ(cl, library(gsubfn)))
invisible(clusterEvalQ(cl, library(parallel)))
invisible(clusterEvalQ(cl, library(sqldf)))
invisible(clusterEvalQ(cl, library(reshape)))
invisible(clusterEvalQ(cl, library(reshape2)))
invisible(clusterEvalQ(cl, library(pbapply)))
print("cluster export started")
parallel::clusterExport(cl=cl,c("RSeq_RPKM","OV_LASSO_CNV","OV_All_methylation","TF_in_Rseq","OV_filtered_mirna_0.01_30_gistic"),envir=environment())#envir needed to be correct, export from global aswell as local
print("cluster export done")
#################

#################
# STEP 3: Call the helper file function to create the Trans Gene Predictor Files.
#################

GenefileNames<-list.files("OvarianLASSO/gene_all_mirna", pattern="*.Rda")
parLapply(cl,GenefileNames,processEachGene)
stopCluster(cl)

##########End Of Code###############

