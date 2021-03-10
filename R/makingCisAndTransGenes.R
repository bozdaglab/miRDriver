## miRDriver STEP 2
## Code 2
## Created By: Banabithi Bose
## Date Created: 5/1/2019
makingCisAndTransGenes <-
  function(Genes_coord_bedfile,
           gistic_bedfile,
           mirdirectory = "~") {
    if ((is(Genes_coord_bedfile, "data.frame")) ||
        (is(Genes_coord_bedfile, "matrix"))) {
      Genes_coord_bedfile <- Genes_coord_bedfile
    } else {
      Genes_coord_bedfile <- assay(Genes_coord_bedfile)
    }
    
    if ((is(gistic_bedfile, "data.frame")) ||
        (is(gistic_bedfile, "matrix"))) {
      gistic_bedfile <- gistic_bedfile
    } else{
      gistic_bedfile <- assay(gistic_bedfile)
    }
    
    dir.create(
      file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_RESULTS",
        fsep = .Platform$file.sep
      )
    )
    dir.create(
      file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_Step2",
        "CIS_TRANS_FILES",
        fsep = .Platform$file.sep
      )
    )
    dir.create(
      file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_Step2",
        "Trans_Cis_Files",
        fsep = .Platform$file.sep
      )
    )
    
    fileNames <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step2",
          "UpDownGeneLibrary"
        ),
        pattern = "*.txt"
      )
    Count_Cis_Trans <- data.frame("", "", "")
    colnames(Count_Cis_Trans) <-
      c("Gistic_Region", "cis_genes_count", "trans_genes_count")
    i <- as.integer(seq_along(fileNames))
    MakecistransFile <- function(i) {
      upgene <-
        read.table(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "UpDownGeneLibrary",
            fileNames[i]
          ),
          sep = " ",
          header = TRUE
        )
      genes <- data.frame(upgene$Row.Names)
      colnames(genes) <- "gene_id"
      Gistic_Region1 <- gsub(".*_", "", fileNames[i])
      Gistic_Region <-
        str_sub(Gistic_Region1, 1, str_length(Gistic_Region1) - 4)
      if ((gistic_bedfile == FALSE)[1]) {
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step1",
            "GisticResults",
            "gistic_bedfile.rda"
          )
        )
        region$chr <- gsub("[chr]", "", region$chr)
        Region <- region[region$peak == Gistic_Region,]
      } else{
        Region <- gistic_bedfile[gistic_bedfile$peak == Gistic_Region,]
      }
      if (nrow(Region) != 0) {
        Region$strand <- "+"
        Region <- Region[, c(1, 5, 2, 3, 4)]
        Region$chr <- paste0("chr", Region$chr)
        colnames(Region) <-
          c("seqnames", "strand", "start", "end", "peak")
        query <- Region %>% as_granges()
        
        Genes_coord_bedfile$strand <- "+"
        Genes_coord_bedfile <- Genes_coord_bedfile[, c(1, 5, 2, 3, 4)]
        Genes_coord_bedfile$chr <-
          paste0("chr", Genes_coord_bedfile$chr)
        colnames(Genes_coord_bedfile) <-
          c("seqnames", "strand", "start", "end", "gene")
        subject <- Genes_coord_bedfile %>% as_granges()
        intersect_rng <-
          suppressWarnings(join_overlap_intersect(query, subject))
        Gistic_Genes <- as.data.frame(intersect_rng)
        
        cis_genes <- data.frame(Gistic_Genes$gene)
        if (nrow(cis_genes) == 0) {
          trans_genes <- genes
          colnames(trans_genes) <- "gene_id"
          save(
            trans_genes,
            file = file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step2",
              "CIS_TRANS_FILES",
              paste0("trans_", fileNames[i], ".Rda")
            )
          )
        } else{
          colnames(cis_genes) <- "gene_id"
          save(
            cis_genes,
            file = paste0(
              "~/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/cis_",
              fileNames[i],
              ".Rda"
            )
          )
          trans_genes <- data.frame(genes[-c(cis_genes$gene_id),])
          colnames(trans_genes) <- "gene_id"
          save(
            trans_genes,
            file = file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step2",
              "CIS_TRANS_FILES",
              paste0("trans_", fileNames[i], ".Rda")
            )
          )
        }
        Count_Cis_Trans1 <-
          data.table(cbind(fileNames[i], nrow(unique(cis_genes)), nrow(unique(trans_genes))))
        colnames(Count_Cis_Trans1) <-
          c("Gistic_Region",
            "cis_genes_count",
            "trans_genes_count")
        return(Count_Cis_Trans1)
      }
    }
    Count_Cis_Trans <- rbindlist(lapply(i, MakecistransFile))
    Count_Cis_Trans <- Count_Cis_Trans[-1,]
    Count_Cis_Trans <- unique(Count_Cis_Trans)
    save(
      Count_Cis_Trans,
      file = file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_RESULTS",
        "Count_Cis_Trans.Rda"
      )
    )
    r <- Count_Cis_Trans[Count_Cis_Trans$trans_genes_count == 0,]
    fileNames1 <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step2",
          "CIS_TRANS_FILES"
        ),
        pattern = "*.Rda"
      )
    fileNames2 <- r$Gistic_Region
    for (i in seq_along(fileNames1)) {
      fileNamesD <- gsub("(cis_)", "", fileNames1[i])
      fileNamesD <- gsub("(trans_)", "", fileNamesD)
      fileNamesD <- str_sub(fileNamesD, end = -5)
      if (is.element(fileNamesD, fileNames2)) {
        file.remove(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "CIS_TRANS_FILES",
            fileNames1[i]
          )
        )
      } else
      {
        
      }
    }
    processUpDownGeneFiles <- function(fileNameUpDn) {
      trans_genes <- NULL
      RegionNameOnly <-
        str_sub(gsub(".*_", "", fileNameUpDn), length(fileNameUpDn),-9)
      if (str_sub(fileNameUpDn, 1, 3) == "cis") {
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "CIS_TRANS_FILES",
            fileNameUpDn
          )
        )
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "CIS_TRANS_FILES",
            gsub(".*cis", "trans", fileNameUpDn)
          )
        )
      } else if (str_sub(fileNameUpDn, 1, 5) == "trans")
      {
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "CIS_TRANS_FILES",
            fileNameUpDn
          )
        )
        cis_file_name <- paste0(gsub(".*trans", "cis", fileNameUpDn))
        if (is.element(cis_file_name, fileNames))
        {
          load(
            file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step2",
              "CIS_TRANS_FILES",
              gsub(".*cis", "trans", fileNameUpDn)
            )
          )
        } else
        {
          cis_genes <- data.frame("")
          colnames(cis_genes) <- "gene_id"
        }
      }
      transfile <- str_sub(fileNameUpDn,-80,-9)
      transfile <- gsub("cis_", "trans_", transfile)
      if (!dir.exists(
        paste0(
          mirdirectory,
          "/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",
          RegionNameOnly
        )
      )) {
        dir.create(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            RegionNameOnly
          )
        )
      }
      save(
        cis_genes,
        file = file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step2",
          "Trans_Cis_Files",
          RegionNameOnly,
          "cis_genes.Rda"
        )
      )
      save(
        trans_genes,
        file = file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step2",
          "Trans_Cis_Files",
          RegionNameOnly,
          paste0(transfile, ".Rda")
        )
      )
    }
    fileNames <-
      list.files(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step2",
          "CIS_TRANS_FILES"
        ),
        pattern = "*.Rda"
      )
    lapply(fileNames, processUpDownGeneFiles)
    
    processUpDownGeneFiles2 <- function(fileName) {
      trans_gene <- NULL
      genefiles <- NULL
      genefiles <-
        list.files(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fileName
          )
        )
      transfiles <- paste0("trans_DOWNgene_", fileName, ".Rda")
      transfiles_new <- paste0("trans_UPgene_", fileName, ".Rda")
      if (is.element(transfiles, genefiles)) {
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fileName,
            transfiles
          )
        )
        transgene1 <- trans_genes
        colnames(transgene1) <- "gene_id"
        if (is.element(transfiles_new, genefiles)) {
          load(
            file.path(
              mirdirectory,
              "mirDriverFold",
              "miRDriver_Step2",
              "Trans_Cis_Files",
              fileName,
              transfiles_new
            )
          )
          transgene2 <- trans_genes
          colnames(transgene2) <- "gene_id"
        } else {
          transgene2 <- data.frame("")
          colnames(transgene2) <- "gene_id"
        }
        trans_genes <- rbind(transgene1, transgene2)
        save(
          trans_genes,
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fileName,
            "trans_genes.Rda"
          )
        )
      } else{
        load(
          file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fileName,
            transfiles_new
          )
        )
        transgene2 <- trans_genes
        transgene1 <- data.frame("")
        colnames(transgene2) <- "gene_id"
        colnames(transgene1) <- "gene_id"
        trans_genes <- rbind(transgene1, transgene2)
        save(
          trans_genes,
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step2",
            "Trans_Cis_Files",
            fileName,
            "trans_genes.Rda"
          )
        )
      }
    }
    folderNames <-
      list.files(file.path(
        mirdirectory,
        "mirDriverFold",
        "miRDriver_Step2",
        "Trans_Cis_Files"
      ))
    lapply(folderNames, processUpDownGeneFiles2)
    return(
      paste0(
        "Trans genes and cis genes files produced in folder mirDriverFold/miRDriver_Step2/Trans_Cis_Files"
      )
    )
  }
