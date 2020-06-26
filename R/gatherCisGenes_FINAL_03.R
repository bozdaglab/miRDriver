## miRDriver STEP 2
## Code 4
## Created By: Banabithi Bose
## Date Created: .........
gatherCisGenes <-
  function(ncores = 2,
           mirna_bedfile,
           gistic_bedfile,
           parallel = c('TRUE', 'FALSE')) {
    dir.create("~/mirDriverFold/miRDriver_Step3")
    dir.create("~/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene")
    dir.create("~/mirDriverFold/miRDriver_Step3/DE_GCIS")
    dir.create("~/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion")
    dir.create("~/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files")
    dir.create("~/mirDriverFold/miRDriver_Step3/GMIRNA")
    dir.create("~/mirDriverFold/miRDriver_Step3/gene_all_mirna")
    dir.create("~/mirDriverFold/miRDriver_Step3/gene_all_DEcis")

    processGeneFiles <- function(fileName) {
      cis_genes <- NULL
      tryCatch({
        load(
          paste0(
            "~/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",
            fileName,
            "/trans_genes.Rda"
          )
        )
        load(
          paste0(
            "~/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",
            fileName,
            "/cis_genes.Rda"
          )
        )

        trans_genes <- data.table(trans_genes)
        trans_genes$gene_id
        trans_genes <- trans_genes[trans_genes$gene_id != "", ]

        if (nrow(cis_genes) == 0) {
          print("this region does not have any cis genes")
        } else{
          for (i in 1:length(trans_genes$gene_id))
          {
            tc_de <- cbind(trans_genes[i], t(cis_genes))
            fname <- paste0(trans_genes$gene_id[i], "_", fileName)

            save(
              tc_de,
              file = paste0(
                "~/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/",
                fname,
                ".Rda"
              )
            )
          }
        }
        return()
      }
      , error = function(error_condition) {
        cat(
          paste0(fileName, " : ", error_condition),
          file = "~/mirDriverFold/ERROR_FILES/Step3_TransGene_w_all_cisgenes_from_all_regions.txt",
          sep = "\n",
          append = TRUE
        )
        return()
      }, finally = {

      })
    }
    if (parallel == TRUE) {
      cl <- makeCluster(mc <- getOption("cl.cores", ncores))
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
      ###### MAIN PROGRAM############################################################
      fileNames <-
        list.files("~/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/")
      parLapply(cl, fileNames, processGeneFiles)
      stopCluster(cl)
    }
    if (parallel == FALSE) {
      fileNames <-
        list.files("~/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/")
      lapply(fileNames, processGeneFiles)
    }
    #################### END ##############################################

    fileNames1 <-
      list.files("~/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/")##27374

    ###### CODE 2
    ## Putting same gene files from different regions in a common folder DE_GCIS for each trans gene
    ### First Create the gene folders
    ###CODE START###
    fileList <-
      list.files("~/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/",
                 pattern = "*.Rda")
    files <- data.table(fileList)
    files <- str_sub(files$fileList, 1, 15)
    files <- unique(files)
    files <- data.table(unique(files))

    for (i in 1:length(files$V1)) {
      newdir <- paste0("~/mirDriverFold/miRDriver_Step3/DE_GCIS/", files$V1[i])
      dir.create(newdir)
    }

    #################
    fileList <-
      list.files("~/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/",
                 pattern = "*.Rda")
    for (i in 1:length(fileList)) {
      newlocation = str_sub(fileList[i], 1, 15)

      current_folder <-
        "~/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/"
      new_folder <-
        paste0("~/mirDriverFold/miRDriver_Step3/DE_GCIS/", newlocation)
      file.copy(file.path(current_folder, fileList[i]), new_folder)

    }

    # # ##The file numbers in each newlocation says for each gene how many gistic regions are there.
    # # ### CODE END ####

    ##Code 3
    processGM <- function(GeneFile) {
      tc_de <- NULL
      tryCatch({
        file_list1 <-
          list.files(paste0("~/mirDriverFold/miRDriver_Step3/DE_GCIS/", GeneFile))

        for (file in file_list1) {
          load(paste0(
            "~/mirDriverFold/miRDriver_Step3/DE_GCIS/",
            GeneFile,
            "/",
            file
          ))
          tc <- data.table(t(tc_de))
          colnames(tc) <- as.character(unlist(tc[1, ]))

          tc <- tc[-1, ]
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
        save(TC,
             file = paste0(
               "~/mirDriverFold/miRDriver_Step3/gene_all_DEcis/",
               GeneFile,
               ".Rda"
             ))


      }, error = function(error_condition) {
        cat(
          paste0(GeneFile, " : ", error_condition),
          file = "~/mirDriverFold/ERROR_FILES/Step3_TransGene_w_all_cisgenes_from_all_regions.txt",
          sep = "\n",
          append = TRUE
        )
        return()
      }, finally = {

      })
    }
    if (parallel == TRUE) {
      cl <- makeCluster(mc <- getOption("cl.cores", ncores))
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
        list.files("~/mirDriverFold/miRDriver_Step3/DE_GCIS") ##3634
      parLapply(cl, GeneFileList, processGM)
      stopCluster(cl)
    }
    if (parallel == FALSE) {
      GeneFileList <-
        list.files("~/mirDriverFold/miRDriver_Step3/DE_GCIS") ##3634
      lapply(GeneFileList, processGM)
    }
    GeneFileList1 <-
      list.files("~/mirDriverFold/miRDriver_Step3/gene_all_DEcis") ##3634

    ##Code 4

    ###This code is for merging TCGA mirna expression data with Gistic regions to get miRNAs
    ### PIPELINE2_OV_LASSO/MiRNA_DIRECTORIES/miRNA_by_GisticRegion has miRNA list, this list is nothing but cis miRNA
    ### then use the miRNAexpression data in LASSO_DATAFRAMES
    #Using Gistic region bed file

    if (gistic_bedfile == FALSE) {
      load("~/mirDriverFold/EXAMPLE_DATASETS/gistic_bedfile.rda")
      region$chr <- paste0("chr", region$chr)
      region <- bedr.merge.region(region)##64 valid regions
      region <- bedr.sort.region(region)
    } else{
      region <- gistic_bedfile
      region$chr <- paste0("chr", region$chr)
      region <- bedr.merge.region(region)##64 valid regions
      region <- bedr.sort.region(region)

    }

    mirna_bedfile <- bedr.sort.region(mirna_bedfile)
    overlapping_pairs <- bedr.join.region(mirna_bedfile,
                                          region,
                                          report.n.overlap = TRUE,
                                          check.chr = FALSE)


    colnames(overlapping_pairs) <- c(
      'miRNA.CHROM',
      'miRNA.POS',
      'miRNA.END',
      'miRNA',
      'region.CHROM',
      'region.POS',
      'region.END',
      'region.PEAK',
      'Overlap'
    )


    ########  making files with mirnas for each gistic regions

    z <-
      data.table(str_split_fixed(overlapping_pairs$region.PEAK, ",", 2))
    overlapping_pairs$region.PEAK <- z$V1
    colnames(overlapping_pairs) <- c(
      'miCHROM',
      'miPOS',
      'miEND',
      'miRNA',
      'regionCHROM',
      'regionPOS',
      'regionEND',
      'regionPEAK',
      'Overlap'
    )

    mirna_in_gistic_region <-
      unique(sqldf("SELECT regionPEAK, miRNA from overlapping_pairs order by 1"))
    length(unique(mirna_in_gistic_region$regionPEAK))

    for (i in 1:(nrow(mirna_in_gistic_region)))
    {
      myfile <-
        file.path(
          paste0(
            "~/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion/",
            mirna_in_gistic_region[i, 1],
            ".txt"
          )
        )
      fileConn <- file(myfile, open = "a")
      writeLines(paste(mirna_in_gistic_region[i, 2]), fileConn)
      close(fileConn)

    }

    tcga_mirnas_count_by_gistic_region <-
      sqldf(
        "select regionPEAK, count(regionPEAK) from mirna_in_gistic_region group by regionPEAK"
      )
    save(tcga_mirnas_count_by_gistic_region, file = "~/mirDriverFold/miRDriver_RESULTS/mirnas_count_by_gistic_region.Rda")

    max(tcga_mirnas_count_by_gistic_region$`count(regionPEAK)`)
    min(tcga_mirnas_count_by_gistic_region$`count(regionPEAK)`)

    ###################Reading region with miRNAs list files###########

    ##Code 5

    ## This will use miRNA_by_GisticRegion and EDGER_DIRECTORIES/up_down_gene files to create gene with all miRNA files.


    # CODE 5:
    fileNamesMiRNA <-
      list.files("~/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion/",
                 pattern = "*.txt")##56 files

    processUpDownGeneFiles <- function(fileNames) {
      file_list1 <-
        list.files(paste0("~/mirDriverFold/miRDriver_Step3/DE_GCIS/", fileNames))
      gene <- fileNames
      for (i in 1:length(file_list1)) {
        mirnafilename <- str_sub(file_list1[i], 17)
        mirnafilename <- str_replace(mirnafilename, ".Rda", ".txt")

        if (is.element(mirnafilename, fileNamesMiRNA))
        {
          RegionMiRna <-
            read.csv(
              paste0(
                "~/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion/",
                mirnafilename,
                sep = ""
              ),
              sep = ";",
              header = FALSE
            )
          Lasso_RegionMiRna <- data.table(RegionMiRna)
          colnames(Lasso_RegionMiRna) <- "miRNA_ID"
          gm <- cbind(gene, t(Lasso_RegionMiRna))
          fname <- paste0(gene, "_", mirnafilename)

          save(
            gm,
            file = paste0(
              "~/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/",
              fname,
              ".Rda"
            )
          )

        }
        else
        {
          print(paste0(mirnafilename, " this region has no mirna for ", gene))### think about it
        }

      }

      return()
    }

    ##### MAIN PROGRAM############################################################

    fileList <-
      list.files("~/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/",
                 pattern = "*.Rda")
    files <- data.table(fileList)
    files <- str_sub(files$fileList, 1, 15)
    files <- unique(files)
    lapply(files, processUpDownGeneFiles)
    # ################### END ##############################################

    fileListmir <-
      list.files("~/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/",
                 pattern = "*.Rda")#29056 files
    # ###### CODE 6
    ##The file numbers in each newlocation says for each gene how many gistic regions are there.
    # ## Putting same gene files from different regions in a common folder
    # ### First Create the gene folders
    # # ##CODE START###
    # # ##Created 12042 unique gene folders
    fileList <-
      list.files("~/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/",
                 pattern = "*.Rda")
    files <- data.table(fileList)
    files <- str_sub(files$fileList, 1, 15)
    files <- unique(files)
    files <- data.table(unique(files))
    for (i in 1:length(files$V1)) {
      newdir <- paste0("~/mirDriverFold/miRDriver_Step3/GMIRNA/", files$V1[i])
      dir.create(newdir)
    }
    # #
    # # ### Putting files in gene folders
    # # ################
    fileList <-
      list.files("~/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/",
                 pattern = "*.Rda")
    for (i in 1:length(fileList)) {
      newlocation = str_sub(fileList[i], 1, 15)

      current_folder <-
        "~/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/"
      new_folder <-
        paste0("~/mirDriverFold/miRDriver_Step3/GMIRNA/", newlocation)
      file.copy(file.path(current_folder, fileList[i]), new_folder)
    }

    ### CODE END ####

    ###CODE 7
    ##Merging all files in PIPELINE2_OV_LASSO/GMIRNA , each gene folder Parallel Code

    processGM <- function(GeneFile) {
      tryCatch({
        file_list1 <-
          list.files(paste0("~/mirDriverFold/miRDriver_Step3/GMIRNA/", GeneFile))

        for (file in file_list1) {
          load(paste0(
            "~/mirDriverFold/miRDriver_Step3/GMIRNA/",
            GeneFile,
            "/",
            file
          ))
          gm <- data.table(t(gm))
          colnames(gm) <- as.character(unlist(gm[1, ]))

          gm <- gm[-1, ]
          if (!exists("GM"))
          {
            GM <- gm

          }
          else
          {
            temp_GM <- gm
            GM <- rbind(GM, temp_GM)
            rm(temp_GM)

          }

        }
        GM <- unique(GM)
        save(GM,
             file = paste0(
               "~/mirDriverFold/miRDriver_Step3/gene_all_mirna/",
               GeneFile,
               ".Rda"
             ))


      }, error = function(error_condition) {
        cat(
          paste0(GeneFile, " : ", error_condition),
          file = "~/mirDriverFold/ERROR_FILES/Step3_TransGene_w_all_mirnagenes_from_all_regions.txt",
          sep = "\n",
          append = TRUE
        )
        return()
      }, finally = {

      })
    }

    if (parallel == TRUE) {
      cl <- makeCluster(mc <- getOption("cl.cores", ncores))
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
        list.files("~/mirDriverFold/miRDriver_Step3/GMIRNA") ##3489
      parLapply(cl, GeneFileList, processGM)
      stopCluster(cl)
    }
    if (parallel == FALSE) {
      GeneFileList <-
        list.files("~/mirDriverFold/miRDriver_Step3/GMIRNA") ##3489
      lapply(GeneFileList, processGM)
    }
    ###End
    GeneFileList1 <-
      list.files("~/mirDriverFold/miRDriver_Step3/gene_all_mirna/") ##3489

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
