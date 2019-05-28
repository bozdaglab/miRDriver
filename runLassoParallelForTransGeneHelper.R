######################################
# File: runLassoParallelForTransGene.R
# Created By: Banabithi Bose
# Date Created: 15th Jan 2019
# Objective: This code is the helper file to run LASSO using files from the folder 
#            OvarianLASSO/PipeTransGeneDEPredFile.
#            This code will be called from the runLassoParallelForTransGene.R
#            which is also available in the same location.
######################################
.libPaths()
gc(reset=T)
library(parallel)
library('glmnet')
library(janitor)
library(stringr)
library(plyr)
library(data.table)
library(readr)
library(gsubfn)
library(parallel)
library(sqldf) 
library(reshape)
library(reshape2)
library(pbapply)

processLassoCtr<- function(numCounter,Lx,Ly,gene_id){
  tryCatch(
    {
      cvfit = cv.glmnet(Lx, Ly,alpha=1,standardize=T,standardize.response=F)
      Coeff<-data.frame("","","")
      colnames(Coeff)<-c("coefficients_min","genes_min","predictors_min")
      Coeff_N<-data.frame("","","")
      colnames(Coeff_N)<-c("coefficients_1se","genes_1se","predictors_1se")
      #time<-system.time(cvfit) #didn't work
      coef_data<-coef(cvfit, s = "lambda.min")[-1,]
      x<-data.frame(coef_data)
      x$Gene_ID<-gene_id
      x$Predictors<-row.names(x)
      x<-data.table(x)
      colnames(x)<-c("coefficients_min","genes_min","predictors_min")
      Coeff<- rbind(Coeff,x)
      Coeff<-Coeff[-1,]
      coef_data_N<-coef(cvfit, s = "lambda.1se")[-1,]
      N<-data.frame(coef_data_N)
      N$Gene_ID<-gene_id
      N$Predictors<-row.names(N)
      N<-data.table(N)
      colnames(N)<-c("coefficients_1se","genes_1se","predictors_1se")
      Coeff_N<- rbind(Coeff_N,N)
      Coeff_N<-Coeff_N[-1,]
      Tot_Coeff<-data.frame(Coeff,Coeff_N)
      return(Tot_Coeff)
    }, error = function(error_condition) {
      cat(paste0(gene_id," : ",error_condition),file="OvarianLASSO/ERROR_FILES/LassoErrors.txt",sep="\n",append=TRUE)
      return()
    }, finally={
    }
  )
  gc()
}

processUpDownGeneFiles <- function(fileName,processLassoCtr){
  load(paste0("OvarianLASSO/PipeTransGeneDEPredFile/",fileName, sep = ""))
  ## Sample size varies from 1052 to 1057 by observation
  gene_id<-str_sub(fileName, length(fileName),-5)
  if (colnames(GenePred)[1]=="sample")
  {
    GenePred<-GenePred[,-1]
  }
  else
  {
    GenePred<-GenePred
  }
  
  GenePred<-subset(GenePred, select=which(!duplicated(names(GenePred))))##checking on duplicated colnames
  if (ncol(GenePred)>2)
  {
    Lx<-data.matrix(GenePred[,-1,drop=T])
    Ly<-data.matrix(GenePred[,1])
    
    j <- 1:numCounter
    
    Tot_Coeff_2<-lapply(j,processLassoCtr,Lx,Ly,gene_id)
    Tot_Coeff<-rbindlist(Tot_Coeff_2, use.names=TRUE, fill=TRUE)
    min_Coeff<-Tot_Coeff[,1:3]
    fse_Coeff<-Tot_Coeff[,4:6]
    save(min_Coeff,file=paste0("OvarianLASSO/LassoMinCoeff/",gene_id,".Rda"))
    
  }
  else
  {
    save(GenePred,file=paste0("OvarianLASSO/ERROR_FILES/lasso_error_genes/df_",gene_id,".Rda"))##this will give us genes on which lasso didn't run
  }
  
  rm(list=ls())
  gc()
}
