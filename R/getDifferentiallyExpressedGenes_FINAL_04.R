## miRDriver STEP 2
## Code 1
## Created By: Banabithi Bose
## Date Created: 5/1/2019
getDifferentiallyExpressedGenes <-function(ncore = 2,RNACount, parallel = c('TRUE', 'FALSE'),mirdirectory="~/") {
    if ((is(RNACount, "data.frame")) || (is(RNACount,"matrix"))){
        RNACount<-RNACount
    }else if (is(RNACount,"SummarizedExperiment")){
        RNACount<-assay(RNACount)
    }else{
        #print("")
    }
    
    ##Create all directories needed
    dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2"))
    dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/UpDownGeneLibrary"))
    dir.create(paste0(mirdirectory,"/mirDriverFold/ERROR_FILES"))
    ##DELETING -CN values.txt FILES
    #fileNames_CN <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files/",pattern = "*CN values.txt"))
    #if(length(fileNames_CN) != 0){
    # file.remove(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files/"),fileNames_CN)
    #}
    fileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files/"),pattern = "*.txt")
    Amplification1 <-read.table(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files/",fileNames[1]),sep = ",",header = TRUE)
    Amplification1[, 1] <- str_trim(Amplification1[, 1])
    
    ##Getting common samples
    colnames(RNACount)<- str_sub(colnames(RNACount),1,16)
    X <- data.table(colnames(RNACount))
    X$V1 <- gsub("[.]", "-", X$V1)
    X$V1 <- as.character(X$V1)
    Y <- data.table(Amplification1[, 1])
    Y$V1 <- gsub("[.]", "-", Y$V1)
    Y$V1 <- as.character(Y$V1)
    Y$V1 <- str_sub(Y$V1,1,16)
    commonSample <- data.table(intersect(X$V1, Y$V1))
    ##This common sample is common between gistic samples i.e. cnv samples and rnaseq data.
    colnames(RNACount) <- gsub("[.]", "-", colnames(RNACount))
    Z <- RNACount[, commonSample$V1]
    processEdgeR <- function(fileName) {
        tryCatch({
            fileNameOnly <- str_sub(fileName, length(fileName), -5)
            AmplDelTables <- read.table(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files/",fileName),sep = ",",header = TRUE)
            AmplDelTables[, 1] <- str_trim(AmplDelTables[, 1])
            AmplDelTables$X <- str_sub(AmplDelTables$X,1,16)
            AmplDelTables <- unique(AmplDelTables)
            AmplDelTables$X <- gsub("\\.", "-", AmplDelTables$X)
            AmplDelTables <- AmplDelTables[!duplicated(AmplDelTables$X), ]
            rownames(AmplDelTables) <- AmplDelTables$X
            AmplDelTables <- data.frame(AmplDelTables[, -1])
            AmplDelTables <- AmplDelTables[commonSample$V1, , drop = FALSE]
            AmplDelTables <- AmplDelTables[order(rownames(AmplDelTables)), ]
            Aberr <- factor(AmplDelTables[, 2])
            y <- DGEList(counts = Z, group = Aberr)
            keep <- rowSums(cpm(y) > 1) >= 1
            y <- y[keep, , keep.lib.sizes = FALSE]
            y <- calcNormFactors(y)
            aberation_type <- AmplDelTables$Aberration
            design <- model.matrix( ~ aberation_type)
            y <- estimateDisp(y, design, robust = TRUE)
            fit <- glmQLFit(y, design, robust = TRUE)
            qlf <- glmQLFTest(fit)
            allgenelist <- qlf$table
            allgenelist <-cbind(Row.Names = rownames(allgenelist), allgenelist)
            rownames(allgenelist) <- NULL
            ##print("DE analysis has started")
            UPDEgene <-sqldf("select * from allgenelist where logFC >= 1 and PValue <= 0.05")
            DNDEgene <-sqldf("select * from allgenelist where logFC <= -1  and PValue <= 0.05")
            write.table(UPDEgene,paste(mirdirectory,"/mirDriverFold/miRDriver_Step2/UpDownGeneLibrary/UPgene_",fileNameOnly,".txt",sep = ""))
            write.table(DNDEgene,paste(mirdirectory,"/mirDriverFold/miRDriver_Step2/UpDownGeneLibrary/DOWNgene_",fileNameOnly,".txt",sep = ""))
            ##print("DE analysis has continued")
        }, error = function(error_condition) {
            cat(paste0(fileName, " : ", error_condition),file = paste0(mirdirectory,"/mirDriverFold/ERROR_FILES/Step2_DEanalysisError.txt"),sep = "\n",append = TRUE)
        }, finally = {
        })
    }
    if (parallel == TRUE) {
        if(.Platform$OS.type == "windows"){
            cl <- makeCluster(mc <- getOption("cl.cores", ncore), type = "PSOCK")
        }else{
            cl <- makeCluster(mc <- getOption("cl.cores", ncore), type = "FORK")
        }
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
        invisible(clusterEvalQ(cl, library(edgeR)))
        invisible(clusterEvalQ(cl, library('statmod')))
        invisible(clusterEvalQ(cl, library(pbapply)))
        #print("cluster export started")
        parallel::clusterExport(cl = cl, c("commonSample", "Z"), envir = environment())
        #print("cluster making done, work started")
        fileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files"),pattern = "*.txt")
        parLapplyLB(cl, fileNames, processEdgeR)
        stopCluster(cl)
        #print("cluster stopped")
    }
    if (parallel == FALSE) {
        fileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files"),pattern = "*.txt")
        lapply(fileNames, processEdgeR)
    }
    return(paste0("Differential analysis has been done, output files in folder miRDriver_Step2/UpDownGeneLibrary"))
}
