## miRDriver STEP 2
## Code 1
## Created By: Banabithi Bose
## Date Created: 5/1/2019

getDifferentiallyExpressedGenes <-
  function(ncore = 2,
           RNACount,
           mirdirectory = "~", DEgene = "All") {
    if ((is(RNACount, "data.frame")) || (is(RNACount, "matrix"))) {
      RNACount <- RNACount
    } else{
      RNACount <- assay(RNACount)
    }

    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step2",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step2",
          fsep = .Platform$file.sep
        )
      )
    }
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step2",
      "UpDownGeneLibrary",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step2",
          "UpDownGeneLibrary",
          fsep = .Platform$file.sep
        )
      )
    }
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "ERROR_FILES",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "ERROR_FILES",
          fsep = .Platform$file.sep
        )
      )
    }
    
    
    fileNames <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step1",
          "Regionwise_Gistic_Files"
        ),
        pattern = "*.txt"
      )
    Amplification1 <-
      read.table(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step1",
          "Regionwise_Gistic_Files",
          fileNames[1]
        ),
        sep = ",",
        header = TRUE
      )
    Amplification1[, 1] <- str_trim(Amplification1[, 1])
    
    
    colnames(RNACount) <- str_sub(colnames(RNACount), 1, 16)
    X <- data.table(colnames(RNACount))
    X$V1 <- gsub("[.]", "-", X$V1)
    X$V1 <- as.character(X$V1)
    Y <- data.table(Amplification1[, 1])
    Y$V1 <- gsub("[.]", "-", Y$V1)
    Y$V1 <- as.character(Y$V1)
    Y$V1 <- str_sub(Y$V1, 1, 16)
    commonSample <- data.table(intersect(X$V1, Y$V1))
    
    
    colnames(RNACount) <- gsub("[.]", "-", colnames(RNACount))
    Z <- RNACount[, commonSample$V1]
    processEdgeR <- function(fileName) {
      tryCatch({
        fileNameOnly <- str_sub(fileName, length(fileName),-5)
        AmplDelTables <-
          read.table(
            file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step1",
              "Regionwise_Gistic_Files",
              fileName
            ),
            sep = ",",
            header = TRUE
          )
        AmplDelTables[, 1] <- str_trim(AmplDelTables[, 1])
        AmplDelTables$X <- str_sub(AmplDelTables$X, 1, 16)
        AmplDelTables <- unique(AmplDelTables)
        AmplDelTables$X <- gsub("\\.", "-", AmplDelTables$X)
        AmplDelTables <- AmplDelTables[!duplicated(AmplDelTables$X),]
        rownames(AmplDelTables) <- AmplDelTables$X
        AmplDelTables <- data.frame(AmplDelTables[,-1])
        AmplDelTables <-
          AmplDelTables[commonSample$V1, , drop = FALSE]
        AmplDelTables <-
          AmplDelTables[order(rownames(AmplDelTables)),]
        Aberr <- factor(AmplDelTables[, 2])
        y <- DGEList(counts = Z, group = Aberr)
        keep <- rowSums(cpm(y) > 1) >= 1
        y <- y[keep, , keep.lib.sizes = FALSE]
        y <- calcNormFactors(y)
        aberation_type <- AmplDelTables$Aberration
        design <- model.matrix(~ aberation_type)
        y <- estimateDisp(y, design, robust = TRUE)
        fit <- glmQLFit(y, design, robust = TRUE)
        qlf <- glmQLFTest(fit)
        allgenelist <- qlf$table
        allgenelist <-
          cbind(Row.Names = rownames(allgenelist), allgenelist)
        rownames(allgenelist) <- NULL
        UPDEgene <-
          sqldf("select * from allgenelist where logFC >= 1 and PValue <= 0.05")
        DNDEgene <-
          sqldf("select * from allgenelist where logFC <= -1  and PValue <= 0.05")
        
        ##Giving the option for UP/Down Genes
        if (DEgene == "UP"||DEgene == "All"){
        write.table(
          UPDEgene,
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "UpDownGeneLibrary",
            paste0("UPgene_", fileNameOnly, ".txt")
          )
        )}
        if (DEgene == "Down"||DEgene == "All"){
        write.table(
          DNDEgene,
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "UpDownGeneLibrary",
            paste0("DOWNgene_", fileNameOnly, ".txt")
          )
        )}
      }, error = function(error_condition) {
        cat(
          paste0(fileName, " : ", error_condition),
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "ERROR_FILES",
            "Step2_DEanalysisError.txt",
            fsep = .Platform$file.sep
          ),
          sep = "\n",
          append = TRUE
        )
      }, finally = {
        
      })
    }
    if (ncore > 1) {
      cl <-
        makeCluster(mc <- getOption("cl.cores", ncore))
      invisible(clusterEvalQ(cl, library('glmnet')))
      invisible(clusterEvalQ(cl, library(janitor)))
      invisible(clusterEvalQ(cl, library(stringr)))
      invisible(clusterEvalQ(cl, library(data.table)))
      invisible(clusterEvalQ(cl, library(gsubfn)))
      invisible(clusterEvalQ(cl, library(parallel)))
      invisible(clusterEvalQ(cl, library(sqldf)))
      invisible(clusterEvalQ(cl, library(reshape)))
      invisible(clusterEvalQ(cl, library(reshape2)))
      invisible(clusterEvalQ(cl, library(edgeR)))
      invisible(clusterEvalQ(cl, library('statmod')))
      parallel::clusterExport(cl = cl, c("commonSample", "Z"), envir = environment())
      fileNames <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step1",
            "Regionwise_Gistic_Files"
          ),
          pattern = "*.txt"
        )
      parLapplyLB(cl, fileNames, processEdgeR)
      stopCluster(cl)
    } else {
      fileNames <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step1",
            "Regionwise_Gistic_Files"
          ),
          pattern = "*.txt"
        )
      lapply(fileNames, processEdgeR)
    }
    return(
      paste0(
        "Differential analysis has been done, output files in folder miRDriver_Step2/UpDownGeneLibrary"
      )
    )
  }