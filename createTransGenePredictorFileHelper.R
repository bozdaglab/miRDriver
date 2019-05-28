######################################
# File: createTransGenePredictorFileHelper.R
# Created By: Banabithi Bose
# Date Created: 15th May 2019
# Objective: This code is the helper file to create the Trans Gene Predictor Files.
#            This code will be called from the createTransGenePredictorFile.R
#            which is also available in the same location.
######################################
gc(reset=T)
.libPaths()
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
processEachGene <- function(GeneFile){
  tryCatch(
    {
      load(paste0("OvarianLASSO/gene_all_mirna/",GeneFile))
      colnames(GM)<-"miRNA_ID"
      mirna<-data.frame(colnames(OV_filtered_mirna_0.01_30_gistic[,-1]))
      colnames(mirna)<-"miRNA"
      mirna$miRNA<-as.character(mirna$miRNA)
      
      common<- intersect(GM$miRNA_ID,mirna$miRNA )
      GeneMiRNA1<-OV_filtered_mirna_0.01_30_gistic[,-1][, ..common]
      
      if (nrow(GeneMiRNA1)==0){
        print(paste0(GeneFile," this gene does not have any filtered mirna for lasso"))
      }
      else { 
        sample<-data.table(OV_filtered_mirna_0.01_30_gistic$sample)
        colnames(sample)<-"sample"
        GeneMiRNA <- cbind(sample,GeneMiRNA1)
        GeneMiRNA$sample<-str_sub(GeneMiRNA$sample,-40,-13)
        
        ## Merging each gene file with cis gene and rnaseq
        load(paste0("OvarianLASSO/gene_all_DEcis/",GeneFile))
        colnames(TC)<-"cis_gene"
        all_genes<-data.frame(colnames(RSeq_RPKM[,-1]))
        colnames(all_genes)<-"gene"
        common1 <- intersect(TC$cis_gene,all_genes$gene)
        CisGene1<-RSeq_RPKM[,-1][, ..common1]
        sample1<-data.table(gsub("[.]","-",RSeq_RPKM$sample))
        colnames(sample1)<-"sample"
        CisGene <- cbind(sample1,CisGene1)
        CisGene$sample<-str_sub(CisGene$sample,-40,-13)
        
        
        ##Merging each genes and Methylation
        Geneid<-str_sub(GeneFile,1,15)
        colnames(OV_All_methylation)<-str_sub(colnames(OV_All_methylation),-40,-13)
        m_sample<-data.table(colnames(OV_All_methylation))
        colnames(m_sample)<-"sample"
        
        V1<-OV_All_methylation[Geneid,]
        x<-V1[rowSums(is.na(V1[ , 1:623])) == 0, ]##13414 unique genes
        if (nrow(x)==0){
          print("this gene does not have beta value")
          V1<-m_sample
        }
        else{
          V1<-melt(V1)
          colnames(V1)<-c("sample","beta")
        }
        
        ##Merging each genes and RnaSeq
        RnaSeq<-data.frame(RSeq_RPKM)
        RnaSeq$sample<-gsub("[.]","-",RnaSeq$sample)
        Y<-RSeq_RPKM[,..Geneid]
        Y$sample<-RnaSeq$sample
        colnames(Y)<-c("rseq","sample")
        Y$sample <- str_sub(Y$sample,-40,-13)
        
        ##Merging genes and Gene Centric CNV
        colnames(OV_LASSO_CNV)<-str_sub(colnames(OV_LASSO_CNV),-40,-13)
        V2<-OV_LASSO_CNV[Geneid,]
        V2<-melt(V2)
        colnames(V2)<-c("sample","cnv")
        
        #Merging Upgenes and TF with expression
        V3<-TF_in_Rseq[TF_in_Rseq$ensembl_gene_id==Geneid,]
        tf_sample<-data.table(RSeq_RPKM$sample)
        colnames(tf_sample)<-"sample"
        tf_sample$sample<-gsub("[.]","-",str_sub(tf_sample$sample,-40,-13))
        
        if (nrow(V3)==0){
          print("this gene does not have TF")
          V3<-tf_sample
        }
        else { 
          tf<-V3$TF
          V3<-RSeq_RPKM[, ..tf]
          V3$sample<-tf_sample$sample
          V3$sample<-gsub("[.]","-",V3$sample)
        }
        
        ## Gene Mirna data
        if (ncol(GeneMiRNA)==0 ){
          V4 = data.frame(GeneMiRNA$sample)
          colnames(V4)=="sample"
        }
        else{
          V4<-GeneMiRNA
        }
        ## Cis gene data
        if (ncol(CisGene)==0 ){
          V5 = data.frame(GeneMiRNA$sample)
          colnames(V5)="sample"
        }
        else{
          V5<-CisGene
        }
        
        ### Making Gene Predictors file
        Y$sample<-as.character(Y$sample)
        V1$sample<-as.character(V1$sample)
        V2$sample<-as.character(V2$sample)
        V3$sample<-as.character(V3$sample)
        V4$sample<-as.character(V4$sample)
        V5$sample<-as.character(V5$sample)
        Y<-data.frame(unique(Y))
        V1<-data.frame(unique(V1))
        V2<-data.frame(unique(V2))
        V3<-data.frame(unique(V3))
        V4<-data.frame(unique(V4))
        V5<-data.frame(unique(V5))
        dfs<-list(Y,V1,V2,V3,V4,V5)
        df<-Reduce(function(...) merge(..., by = "sample", x.all = TRUE), dfs) #1052 unique samples
        df <- remove_empty_cols(df)
        GenePred<-df[complete.cases(df), ]
        save(GenePred,file=paste0("OvarianLASSO/PipeTransGeneDEPredFile/",Geneid,".Rda"))
      }
    }, error = function(error_condition) {
      cat(paste0(GeneFile," : ",error_condition),file="OvarianLASSO/ERROR_FILES/transGenePredictorError.txt",sep="\n",append=TRUE)
    }, finally={
    }
  )
}
##########End Of Code###############