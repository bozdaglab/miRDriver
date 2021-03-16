## miRDriver STEP 3
## Code 1
## Created By: Banabithi Bose
## Date Created: 5/1/2019

gatherCisGenes <-
  function(ncore = 2,
           mirna_bedfile,
           gistic_bedfile,
           mirdirectory = "~") {
    if ((is(mirna_bedfile, "data.frame")) ||
        (is(mirna_bedfile, "matrix"))) {
      mirna_bedfile <- mirna_bedfile
    } else{
      mirna_bedfile <- assay(mirna_bedfile)
    }
    if ((is(gistic_bedfile, "data.frame")) ||
        (is(gistic_bedfile, "matrix"))) {
      gistic_bedfile <- gistic_bedfile
    } else {
      gistic_bedfile <- assay(gistic_bedfile)
    }
    
    tc_de <- NULL
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          fsep = .Platform$file.sep
        )
      )
    } 
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      "transgene_w_DEcisgene",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "transgene_w_DEcisgene",
          fsep = .Platform$file.sep
        )
      )
    } 
    
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      "DE_GCIS",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "DE_GCIS",
          fsep = .Platform$file.sep
        )
      )
    } 
    
    
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      "miRNA_by_GisticRegion",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "miRNA_by_GisticRegion",
          fsep = .Platform$file.sep
        )
      )
    }
    
    
    
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      "Trans_Mirna_Files",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "Trans_Mirna_Files",
          fsep = .Platform$file.sep
        )
      )
    }   
    
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      "GMIRNA",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "GMIRNA",
          fsep = .Platform$file.sep
        )
      )
    } 
    
    
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      "gene_all_mirna",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "gene_all_mirna",
          fsep = .Platform$file.sep
        )
      )
    }
    
    
    
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step3",
      "gene_all_DEcis",
      fsep = .Platform$file.sep
    ))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "gene_all_DEcis",
          fsep = .Platform$file.sep
        )
      )
    }
    
    
    processGeneFiles <- function(fileName) {
      cis_genes <- NULL
      tryCatch({
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fileName,
            "trans_genes.Rda",
            fsep = .Platform$file.sep
          )
        )
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fileName,
            "cis_genes.Rda",
            fsep = .Platform$file.sep
          )
        )
        trans_genes <- data.table(trans_genes)
        trans_genes <- trans_genes[trans_genes$gene_id != "",]
        if (nrow(cis_genes) != 0) {
          for (i in seq_along(trans_genes$gene_id))
          {
            tc_de <- cbind(trans_genes[i], t(cis_genes))
            fname <- paste0(trans_genes$gene_id[i], "_", fileName)
            save(
              tc_de,
              file = file.path(
                mirdirectory,
                "mirDriverFold",
                "miRDriver_Step3",
                "transgene_w_DEcisgene",
                paste0(fname, ".Rda")
              )
            )
          }
        }
        return()
      }, error = function(error_condition) {
        cat(
          paste0(fileName, " : ", error_condition),
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "ERROR_FILES",
            "Step3_TransGene_w_all_cisgenes_from_all_regions.txt",
            fsep = .Platform$file.sep
          ),
          sep = "\n",
          append = TRUE
        )
        return()
      }, finally = {
        
      })
    }
    if (ncore > 2) {
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
      clusterExport(cl = cl, "processGeneFiles", envir = environment())
      
      fileNames <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fsep = .Platform$file.sep
          )
        )
      parLapply(cl, fileNames, processGeneFiles)
      stopCluster(cl)
    } else{
      fileNames <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fsep = .Platform$file.sep
          )
        )
      lapply(fileNames, processGeneFiles)
    }
    
    fileList <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "transgene_w_DEcisgene",
          fsep = .Platform$file.sep
        ),
        pattern = "*.Rda"
      )
    files <- data.table(fileList)
    files <-
      data.table(str_split_fixed(files$fileList, "_", 2))[, 1, drop = FALSE]
    files <- unique(files)
    files <- data.table(unique(files))
    for (i in seq_along(files$V1)) {
      newdir <-
        file.path(mirdirectory,
                  "mirDriverFold",
                  "miRDriver_Step3",
                  "DE_GCIS",
                  files$V1[i])
      
      if (!file.exists(newdir)) {
        dir.create(newdir)
      } 
      
      
      
    }
    current_folder <-
      file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_Step3",
        "transgene_w_DEcisgene",
        fsep = .Platform$file.sep
      )
    fileList <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "transgene_w_DEcisgene",
          fsep = .Platform$file.sep
        ),
        pattern = "*.Rda"
      )
    
    for (i in seq_along(fileList)) {
      newlocation = unique(data.table(str_split_fixed(fileList[i], "_", 2))[, 1])$V1
      newfolder <-
        file.path(mirdirectory,
                  "mirDriverFold",
                  "miRDriver_Step3",
                  "DE_GCIS",
                  newlocation)
      load(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "transgene_w_DEcisgene",
          fileList[i]
        )
      )
      save(tc_de, file = file.path(newfolder, fileList[i]))
    }
    
    processGM <- function(GeneFile) {
      tc_de <- NULL
      tryCatch({
        file_list1 <-
          list.files(
            file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step3",
              "DE_GCIS",
              GeneFile
            )
          )
        for (file in file_list1) {
          load(
            file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step3",
              "DE_GCIS",
              GeneFile,
              file
            )
          )
          tc <- data.table(t(tc_de))
          colnames(tc) <- as.character(unlist(tc[1,]))
          tc <- tc[-1,]
          if (!exists("TC"))
          {
            TC <- tc
          } else
          {
            temp_TC <- tc
            TC <- rbind(TC, temp_TC)
            rm(temp_TC)
          }
        }
        TC <- unique(TC)
        save(
          TC,
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step3",
            "gene_all_DEcis",
            paste0(GeneFile, ".Rda")
          )
        )
      }, error = function(error_condition) {
        cat(
          paste0(GeneFile, " : ", error_condition),
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "ERROR_FILES",
            "Step3_TransGene_w_all_cisgenes_from_all_regions.txt",
            fsep = .Platform$file.sep
          ),
          sep = "\n",
          append = TRUE
        )
        return()
      }, finally = {
        
      })
    }
    if (ncore > 2) {
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
      clusterExport(cl = cl, "processGM", envir = environment())
      GeneFileList <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step3",
            "DE_GCIS",
            fsep = .Platform$file.sep
          )
        )
      parLapply(cl, GeneFileList, processGM)
      stopCluster(cl)
    } else{
      GeneFileList <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step3",
            "DE_GCIS",
            fsep = .Platform$file.sep
          )
        )
      lapply(GeneFileList, processGM)
    }
    GeneFileList1 <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "gene_all_DEcis",
          fsep = .Platform$file.sep
        )
      )
    
    region <- gistic_bedfile
    
    region$strand <- "+"
    region <- region[, c(1, 5, 2, 3, 4)]
    region$chr <- paste0("chr", region$chr)
    colnames(region) <-
      c("seqnames", "strand", "start", "end", "peak")
    query <- region %>% as_granges()
    
    
    mirna_bedfile$strand <- "+"
    mirna_bedfile <- mirna_bedfile[, c(1, 5, 2, 3, 4)]
    mirna_bedfile$chr <- paste0("chr", mirna_bedfile$chr)
    colnames(mirna_bedfile) <-
      c("seqnames", "strand", "start", "end", "mirna")
    subject <- mirna_bedfile %>% as_granges()
    intersect_rng <- join_overlap_intersect(query, subject)
    overlapping_pairs <- as.data.frame(intersect_rng)
    
    overlapping_pairs <- overlapping_pairs[, c(6, 7)]
    colnames(overlapping_pairs) <- c("region.PEAK", "mirna")
    
    z <-
      data.table(str_split_fixed(overlapping_pairs$region.PEAK, ",", 2))
    overlapping_pairs$region.PEAK <- z$V1
    colnames(overlapping_pairs) <- c("regionPEAK", "miRNA")
    
    mirna_in_gistic_region <-
      unique(sqldf("SELECT regionPEAK, miRNA from overlapping_pairs order by 1"))
    for (i in seq_len(nrow(mirna_in_gistic_region)))
    {
      myfile <-
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "miRNA_by_GisticRegion",
          paste0(mirna_in_gistic_region[i, 1], ".txt")
        )
      fileConn <- file(myfile, open = "a")
      writeLines(paste(mirna_in_gistic_region[i, 2]), fileConn)
      close(fileConn)
    }
    tcga_mirnas_count_by_gistic_region <-
      sqldf(
        "select regionPEAK, count(regionPEAK) from mirna_in_gistic_region group by regionPEAK"
      )
    save(
      tcga_mirnas_count_by_gistic_region,
      file = file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_RESULTS",
        "mirnas_count_by_gistic_region.Rda"
      )
    )
    
    fileNamesMiRNA <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "miRNA_by_GisticRegion"
        ),
        pattern = "*.txt"
      )
    
    processUpDownGeneFiles <- function(fileNames) {
      f <- NULL
      file_list1 <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step3",
            "DE_GCIS",
            fileNames
          )
        )
      gene <- fileNames
      for (i in seq_along(file_list1)) {
        mirnafilename <- str_split_fixed(file_list1[i], "_", 2)[, 2]
        mirnafilename <- str_replace(mirnafilename, ".Rda", ".txt")
        if (is.element(mirnafilename, fileNamesMiRNA))
        {
          RegionMiRna <-
            read.csv(
              file.path(
                mirdirectory,
                "mirDriverFold",
                "miRDriver_Step3",
                "miRNA_by_GisticRegion",
                mirnafilename
              ),
              sep = ";",
              header = FALSE
            )
          Lasso_RegionMiRna <- data.table(RegionMiRna)
          colnames(Lasso_RegionMiRna) <- "miRNA_ID"
          gm <- data.table(cbind(gene, t(Lasso_RegionMiRna)))
          fname <- paste0(gene, "_", mirnafilename)
          fname <- str_replace(fname, ".txt", "")
          save(
            gm,
            file = file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step3",
              "Trans_Mirna_Files",
              paste0(fname, ".Rda")
            )
          )
        }
      }
    }
    
    fileList <-
      list.files(file.path(mirdirectory, "mirDriverFold", "miRDriver_Step3", "DE_GCIS"))
    files <- unique(data.table(fileList))
    lapply(files$fileList, processUpDownGeneFiles)
    
    fileList <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "Trans_Mirna_Files"
        ),
        pattern = "*.Rda"
      )
    files <- data.table(fileList)
    files <-
      data.table(str_split_fixed(files$fileList, "_", 2))[, 1, drop = FALSE]
    files <- unique(files)
    for (i in seq_along(files$V1)) {
      newdir <-
        file.path(mirdirectory,
                  "mirDriverFold",
                  "miRDriver_Step3",
                  "GMIRNA",
                  files$V1[i])
      if (!file.exists(newdir)) {
        dir.create(newdir)
      } 
    }
    for (i in seq_along(fileList)) {
      newlocation = unique(data.table(str_split_fixed(fileList[i], "_", 2))[, 1])$V1
      current_folder <-
        file.path(mirdirectory,
                  "mirDriverFold",
                  "miRDriver_Step3",
                  "Trans_Mirna_Files")
      new_folder <-
        file.path(mirdirectory,
                  "mirDriverFold",
                  "miRDriver_Step3",
                  "GMIRNA",
                  newlocation)
      file.copy(file.path(current_folder, fileList[i]), new_folder)
    }
    
    processGM <- function(GeneFile) {
      tryCatch({
        file_list1 <-
          list.files(
            file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step3",
              "GMIRNA",
              GeneFile
            )
          )
        for (file in file_list1) {
          load(
            file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step3",
              "GMIRNA",
              GeneFile,
              file
            )
          )
          gm <- data.table(t(gm))
          colnames(gm) <- as.character(unlist(gm[1,]))
          gm <- gm[-1,]
          if (!exists("GM"))
          {
            GM <- gm
          } else
          {
            temp_GM <- gm
            GM <- rbind(GM, temp_GM)
            rm(temp_GM)
          }
        }
        GM <- unique(GM)
        save(
          GM,
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step3",
            "gene_all_mirna",
            paste0(GeneFile, ".Rda")
          )
        )
      }, error = function(error_condition) {
        cat(
          paste0(fileName, " : ", error_condition),
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "ERROR_FILES",
            "Step3_TransGene_w_all_mirnagenes_from_all_regions.txt",
            fsep = .Platform$file.sep
          ),
          sep = "\n",
          append = TRUE
        )
        return()
      }, finally = {
        
      })
    }
    if (ncore > 2) {
      cl <-
        makeCluster(mc <- getOption("cl.cores", ncore))
      
      clusterEvalQ(cl, library('glmnet'))
      clusterEvalQ(cl, library(janitor))
      clusterEvalQ(cl, library(stringr))
      clusterEvalQ(cl, library(plyr))
      clusterEvalQ(cl, library(data.table))
      clusterEvalQ(cl, library(readr))
      clusterEvalQ(cl, library(gsubfn))
      clusterEvalQ(cl, library(parallel))
      clusterEvalQ(cl, library(sqldf))
      clusterEvalQ(cl, library(reshape))
      clusterEvalQ(cl, library(reshape2))
      clusterEvalQ(cl, library(pbapply))
      clusterExport(cl = cl, "processGM", envir = environment())
      GeneFileList <-
        list.files(file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "GMIRNA"
        ))
      parLapply(cl, GeneFileList, processGM)
      stopCluster(cl)
    } else{
      GeneFileList <-
        list.files(file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step3",
          "GMIRNA"
        ))
      lapply(GeneFileList, processGM)
    }
    GeneFileList1 <-
      list.files(file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_Step3",
        "gene_all_mirna"
      ))
    return(
      paste0(
        "mirDriverFold/miRDriver_Step3/gene_all_mirna and gene_all_DEcis are filled with cis genes and cis mirnas"
      )
    )
    return(
      paste0(
        "gene_all_mirna and gene_all_DEcis folders are filled in mirDriverFold/miRDriver_Step3"
      )
    )
  }
