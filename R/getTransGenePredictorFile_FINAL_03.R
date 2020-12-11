## miRDriver STEP 4
## Code 1
## Created By: Banabithi Bose
## Date Created: 5/1/2019
getTransGenePredictorFile <-function(ncores = 2,methylationData,RSeqData,CNVData,TFData,mirnaData,parallel = c('TRUE', 'FALSE'),mirdirectory="~/") {
    
    
    if (is(methylationData,"SummarizedExperiment")){
        methylationData<-assay(methylationData)
    }else {
        methylationData<-methylationData
    }


    if (is(RSeqData,"SummarizedExperiment")){
        RSeqData<-assay(RSeqData)
    }else {
        RSeqData<-RSeqData
    }

    if (is(CNVData,"SummarizedExperiment")){
        CNVData<-assay(CNVData)
    }else {
        CNVData<-CNVData
    }

    if (is(TFData,"SummarizedExperiment")){
        TFData<-assay(TFData)
    }else {
        TFData<-TFData
    }

    
    if (is(mirnaData,"SummarizedExperiment")){
        mirnaData<-assay(mirnaData)
    }else {
        mirnaData<-mirnaData
    }


    ##..Geneid <-..common<- 
    ..tf <- processEachGene <- NULL
    dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step4"))
    dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step4/TransGenePredFile"))
    #GeneFile<-GeneFileNames[1]
    processEachGene <- function(GeneFile) {
        tryCatch({
            methylation <- methylationData
            RSeq <- RSeqData
            CNV <- CNVData
            TF <- TFData
            MG <- mirnaData
            load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_mirna/",GeneFile))
            colnames(GM) <- "miRNA_ID"
            mirna <-data.frame(colnames(MG[, -1]))
            colnames(mirna) <- "miRNA"
            mirna$miRNA <- as.character(mirna$miRNA)
            common <- intersect(GM$miRNA_ID, mirna$miRNA)
            GeneMirnaData <- MG[, ..common]
            if (nrow(GeneMirnaData) == 0) {
                #print(paste0(GeneFile, " no mirna, LASSO will not run for this trans gene"))
            } else
            {
                GeneMirnaData <- cbind(MG[, 1, drop = FALSE], MG[, ..common])
                MG<-data.table(MG)
                #MG<-data.frame(MG)
                GeneMirnaData1 <- MG[, ..common]
                #GeneMirnaData1 <- MG[, common]
                ## Merging each gene file with cis gene and rnaseq
                load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_DEcis/",GeneFile))
                colnames(TC) <- "cis_gene"
                all_genes <- data.frame(colnames(RSeq[, -1]))
                colnames(all_genes) <- "gene"
                common <- intersect(TC$cis_gene, all_genes$gene)
                
                CisGene <- RSeq[, common]
                if (ncol(CisGene) == 0) {
                    #print(paste0(GeneFile, " no cis genes for this trans gene"))
                    CisGeneData <- RSeq[, 1, drop = FALSE]
                } else {
                    CisGeneData <- cbind(RSeq[, 1, drop = FALSE], RSeq[, common])
                }
                ## Finding Methylation beta values for each trans gene
                Geneid <- str_sub(GeneFile, end = -5)
                if (is.element(Geneid, colnames(methylation))) {
                    MethGene <- methylation[, ..Geneid]
                    MethGene <- na.omit(MethGene)
                    if (nrow(MethGene) == 0) {
                        #print(paste0(GeneFile, " no methylation values for this trans gene"))
                        MethGeneData <- methylation[, 1, drop = FALSE]
                    } else{
                        MethGeneData <-cbind(methylation[, 1, drop = FALSE], methylation[, ..Geneid])
                        colnames(MethGeneData) <- c("sample", "beta")
                    }
                } else{
                    MethGeneData <- methylation[, 1, drop = FALSE]
                }
                ##Merging trans genes and Gene Centric CNV
                if (is.element(Geneid, colnames(CNV))) {
                    cnvGene <- CNV[, ..Geneid]
                    cnvGene <- na.omit(cnvGene)
                    if (nrow(cnvGene) == 0) {
                        #print(paste0(GeneFile, " no cnv values for this trans gene"))
                        CNVGeneData <- CNV[, 1, drop = FALSE]
                    } else{
                        CNVGeneData <- cbind(CNV[, 1, drop = FALSE], CNV[, ..Geneid])
                        colnames(CNVGeneData) <- c("sample", "cnv")
                    }
                } else{
                    CNVGeneData <- CNV[, 1, drop = FALSE]
                }
                #Merging transgene and TF with expression
                TF<-data.frame(TF)
                tfGene <- TF[TF$gene_id == Geneid, ]
                if (nrow(tfGene) == 0) {
                    #print("this trans gene does not have TF")
                    tfGeneData <- RSeq[, 1, drop = FALSE]
                } else{
                    tf <- tfGene$TF
                    if (is.element(tf, colnames(RSeq))) {
                        tfGeneData <- cbind(RSeq[, 1, drop = FALSE], RSeq[, ..tf])
                    } else{
                        tfGeneData <- RSeq[, 1, drop = FALSE]
                    }
                }
                if (is.element(Geneid, colnames(RSeq))) {
                    Y <-cbind(RSeq[, 1, drop = FALSE], RSeq[, Geneid])
                } else{
                    #print("this transgene isn't present in rnaseq data")
                }
                V1 <- GeneMirnaData
                V2 <- CisGeneData
                V3 <- MethGeneData
                V4 <- CNVGeneData
                V5 <- tfGeneData
                ## Making Gene Predictors file
                Y$sample <- as.character(Y$sample)
                V1$sample <- as.character(V1$sample)
                V2$sample <- as.character(V2$sample)
                V3$sample <- as.character(V3$sample)
                V4$sample <- as.character(V4$sample)
                V5$sample <- as.character(V5$sample)
                Y <- data.frame(unique(Y))
                colnames(Y)<-c("sample",Geneid)
                V1 <- data.frame(unique(V1))
                V2 <- data.frame(unique(V2))
                V3 <- data.frame(unique(V3))
                V4 <- data.frame(unique(V4))
                V5 <- data.frame(unique(V5))
                dfs <- list(Y, V1, V2, V3, V4, V5)
                df <-Reduce(function(...)merge(..., by = "sample", x.all = TRUE), dfs)
                #df <- remove_empty(df)
                df <- remove_empty(df, c("rows", "cols"))
                GenePred <- df[complete.cases(df),]
                save(GenePred,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step4/TransGenePredFile/",Geneid,".Rda"))
            }
        }, error = function(error_condition) {
            cat(paste0(GeneFile, " : ", error_condition),file = paste0(mirdirectory,"/mirDriverFold/ERROR_FILES/Step4_transGenePredictorError.txt"),sep = "\n",append = TRUE)
        }, finally = {
        })
    }
    ##########End of function###############
    if (parallel == TRUE) {
        ##Make parallel cluster and export the dataframes to all the clusters
        #################
        #print("cluster making started")
        if(.Platform$OS.type == "windows"){
            cl <- makeCluster(mc <- getOption("cl.cores", ncores), type = "PSOCK")
        }else{
            cl <- makeCluster(mc <- getOption("cl.cores", ncores), type = "FORK")
        }
        invisible(clusterEvalQ(cl, library(glmnet)))
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
        #print("cluster export started")
        parallel::clusterExport(cl = cl,c("methylationData","RSeqData","CNVData","TFData","mirnaData"),envir = environment())
        #print("cluster export done")
        #################
        GeneFileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_mirna/"), pattern = "*.Rda")
        parLapply(cl, GeneFileNames, processEachGene)
        stopCluster(cl)
    }
    if (parallel == FALSE) {
        GeneFileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_mirna/"), pattern = "*.Rda")
        lapply(GeneFileNames, processEachGene)
    }
    return(paste0("All predictors are gathered for trans genes in mirDriverFold/miRDriver_Step4/TransGenePredFile"))
}
