## miRDriver STEP 4
## Code 2
## Created By: Banabithi Bose
## Date Created: 5/1/2019
lassoParallelForTransGene <-
  function(ncore = 1,
           numCounter = 100,
           Nfolds = 10,
           nonZeroPercent = 70,
           mirdirectory = "~") {
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_RESULTS",
      "LassoMinCoeff"
    ))) {
      dir.create(file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_RESULTS",
        "LassoMinCoeff"
      ))
    } 

    processLassoCtr <- function(numCounter, Lx, Ly, gene_id) {
      tryCatch({
        cvfit <-
          suppressWarnings(
            cv.glmnet(
              Lx,
              Ly,
              alpha = 1,
              standardize = TRUE,
              standardize.response = FALSE,
              nfolds = Nfolds
            )
          )
        Coeff <- data.frame("", "", "")
        colnames(Coeff) <-
          c("coefficients_min", "genes_min", "predictors_min")
        Coeff_N <- data.frame("", "", "")
        colnames(Coeff_N) <-
          c("coefficients_1se", "genes_1se", "predictors_1se")
        coef_data <- coef(cvfit, s = "lambda.min")[-1,]
        x <- data.frame(coef_data)
        x$Gene_ID <- gene_id
        x$Predictors <- row.names(x)
        x <- data.table(x)
        colnames(x) <-
          c("coefficients_min", "genes_min", "predictors_min")
        Coeff <- rbind(Coeff, x)
        Coeff <- Coeff[-1,]
        coef_data_N <- coef(cvfit, s = "lambda.1se")[-1,]
        N <- data.frame(coef_data_N)
        N$Gene_ID <- gene_id
        N$Predictors <- row.names(N)
        N <- data.table(N)
        colnames(N) <-
          c("coefficients_1se", "genes_1se", "predictors_1se")
        Coeff_N <- rbind(Coeff_N, N)
        Coeff_N <- Coeff_N[-1,]
        Tot_Coeff <- data.frame(Coeff, Coeff_N)
        return(Tot_Coeff)
      }, error = function(error_condition) {
        cat(
          paste0(gene_id, " : ", error_condition),
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "ERROR_FILES",
            "Step4_LassoErrors.txt"
          ),
          sep = "\n",
          append = TRUE
        )
        return()
      }, finally = {
        
      })
      return()
    }
    processUpDownGeneFiles <- function(fileName, processLassoCtr) {
      load(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step4",
          "TransGenePredFile",
          fileName
        )
      )
      gene_id <- str_sub(fileName, length(fileName),-5)
      if (colnames(GenePred)[1] == "sample")
      {
        GenePred <- GenePred[,-1]
      } else
      {
        GenePred <- GenePred
      }
      GenePred <-
        subset(GenePred, select = which(!duplicated(names(GenePred))))
      if (ncol(GenePred) > 2)
      {
        Lx <- data.matrix(GenePred[,-1, drop = TRUE])
        Ly <- data.matrix(GenePred[, 1])
        j <- seq_len(numCounter)
        Tot_Coeff_2 <- lapply(j, processLassoCtr, Lx, Ly, gene_id)
        Tot_Coeff <-
          rbindlist(Tot_Coeff_2, use.names = TRUE, fill = TRUE)
        min_Coeff <- Tot_Coeff[, c(1, 2, 3)]
        fse_Coeff <- Tot_Coeff[, c(4, 5, 6)]
        save(
          min_Coeff,
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_RESULTS",
            "LassoMinCoeff",
            paste0(gene_id, ".Rda")
          )
        )
      } else
      {
        save(
          GenePred,
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "ERROR_FILES",
            "Step4_Lasso_",
            paste0(gene_id, ".Rda")
          )
        )
      }
      rm(list = ls())
      gc()
    }
    if (ncore > 1) {
      cl <-
        makeCluster(mc <- getOption("cl.cores", ncore))
      
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
      invisible(clusterEvalQ(cl, library(pbapply)))
      parallel::clusterExport(
        cl = cl,
        c("numCounter", "Nfolds", "processLassoCtr"),
        envir = environment()
      )
      fileNames <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step4",
            "TransGenePredFile"
          ),
          pattern = "*.Rda"
        )
      parLapplyLB(cl, fileNames, processUpDownGeneFiles, processLassoCtr)
      stopCluster(cl)
    } else {
      fileNames <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step4",
            "TransGenePredFile"
          ),
          pattern = "*.Rda"
        )
      lapply(fileNames, processUpDownGeneFiles, processLassoCtr)
    }
    processMirGene <- function(fileNames, nonZeroTimes) {
      load(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_RESULTS",
          "LassoMinCoeff",
          fileNames
        )
      )
      gene_id <- str_sub(fileNames, end = -5)
      nx <-
        sqldf(
          "select genes_min,predictors_min,count(coefficients_min) coeff_min_count, avg(coefficients_min)avg_coeff from min_Coeff where coefficients_min <> 0 group by genes_min,predictors_min "
        )
      selected_features <-
        sqldf(paste0("select * from nx where coeff_min_count >=", nonZeroTimes))
      if (nrow(selected_features) > 0)
      {
        selected_features$genes_min <- gene_id
        gene_with_mirana <-
          sqldf("select * from selected_features where predictors_min like '%hsa%'")
        d1 <- data.table(gene_with_mirana[,-3])
        d1 <- unique(d1)
        colnames(d1) <- c("genes", "miRNA", "av_coeff")
        Mir_Gen <- data.table(d1)
        return(Mir_Gen)
      } else
      {
        return()
      }
    }
    nonZeroTimes <- round(numCounter * (nonZeroPercent / 100))
    LASSOfileNames <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_RESULTS",
          "LassoMinCoeff"
        ),
        pattern = "*.Rda"
      )
    Final_gene_mirna_network <-
      rbindlist(lapply(LASSOfileNames, processMirGene, nonZeroTimes))
    Final_gene_mirna_network <- unique(Final_gene_mirna_network)
    save(
      Final_gene_mirna_network,
      file = file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_RESULTS",
        "Final_gene_mirna_network.Rda"
      )
    )
    return(Final_gene_mirna_network)
    
  }
