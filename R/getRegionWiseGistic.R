## miRDriver STEP 1
## Code 1
## Created By: Banabithi Bose
## Date Created: 5/1/2019
getRegionWiseGistic <-
  function(gisticfile,
           gistic_bedfile,
           mirdirectory = "~", Aber = "All") {
    if (!file.exists(file.path(mirdirectory, "mirDriverFold"))) {
      dir.create(file.path(mirdirectory, "mirDriverFold"))
      
    }
    if (!file.exists(file.path(mirdirectory, "mirDriverFold", "miRDriver_Step1"))) {
      dir.create(file.path(mirdirectory, "mirDriverFold", "miRDriver_Step1"))
      
    }
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step1",
      "GisticResults"))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step1",
          "GisticResults",
          fsep = .Platform$file.sep
        ))
    }
    if (!file.exists(file.path(
      mirdirectory,
      "mirDriverFold",
      "miRDriver_Step1",
      "Regionwise_Gistic_Files"))) {
      dir.create(
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step1",
          "Regionwise_Gistic_Files",
          fsep = .Platform$file.sep
        ))
    }
    
    
    Gistic_lesions <-
      read.table(gisticfile, sep = "\t", header = TRUE)
    
    ## Adding the options for Amp/Del/All peaks
    if (Aber == "Amp"){
      Gistic_lesions <- Gistic_lesions[grep("Amp", Gistic_lesions$Unique.Name, ignore.case = FALSE, perl = FALSE, value = FALSE,
                                            fixed = FALSE, useBytes = FALSE, invert = FALSE),]
    }else if (Aber == "Del"){
      Gistic_lesions <- Gistic_lesions[grep("Del", Gistic_lesions$Unique.Name, ignore.case = FALSE, perl = FALSE, value = FALSE,
                                            fixed = FALSE, useBytes = FALSE, invert = FALSE),]
      }
    Gistic_lesions <-
      Gistic_lesions[-grep("CN values", Gistic_lesions$Unique.Name),]
    Gistic_lesions <-  Gistic_lesions[seq_len(nrow(Gistic_lesions)),]
    indx <-
      which(Gistic_lesions == '1' |
              Gistic_lesions == '2', arr.ind = TRUE)
    val <-
      tapply(names(Gistic_lesions)[indx[, 2]], indx[, 1], FUN = toString)
    Gistic_lesions$Aberrated_patients <-
      ifelse(seq_len(nrow(Gistic_lesions)) %in% names(val), val, NA)
    indx1 <- which(Gistic_lesions == '0', arr.ind = TRUE)
    val1 <-
      tapply(names(Gistic_lesions)[indx1[, 2]], indx1[, 1], FUN = toString)
    Gistic_lesions$NonAberrated_patients <-
      ifelse(seq_len(nrow(Gistic_lesions)) %in% names(val1), val1, NA)
    BR <-
      Gistic_lesions[, c(1,
                         2,
                         3,
                         4,
                         5,
                         6,
                         7,
                         8,
                         9,
                         ncol(Gistic_lesions) - 1,
                         ncol(Gistic_lesions))]
    for (j in seq_len(nrow(BR))) {
      BR_Peak <- as.character(BR$Unique.Name[j])
      BR_Aberrated <-
        as.vector(strsplit(c(as.vector(BR$Aberrated[j])), ",")[[1]])
      BR_NonAberrated <-
        as.vector(strsplit(c(as.vector(
          BR$NonAberrated[j]
        )), ",")[[1]])
      myfile <-
        file.path(
          mirdirectory,
          "mirDriverFold",
          "miRDriver_Step1",
          "Regionwise_Gistic_Files",
          paste0(BR_Peak, ".txt")
        )
      fileConn <- file(myfile, open = "a")
      writeLines(paste(" ", ",Patient,", "Aberration"), fileConn)
      for (i in seq_along(BR_Aberrated))
        writeLines(paste(BR_Aberrated[i], ",1,", "aberrated"), fileConn)
      for (i in seq_along(BR_NonAberrated))
        writeLines(paste(BR_NonAberrated[i], ",1,", "non-aberrated"), fileConn)
      close(fileConn)
      if (gistic_bedfile == TRUE) {
        gistic_region1 <- Gistic_lesions[, c(1, 5)]
        x <- str_split_fixed(gistic_region1$Region.Limits, ":", 2)
        x <- data.table(x)
        y <- data.table(str_split_fixed(x$V2, "-", 2))
        z <- data.table(str_split_fixed(y$V2, "p", 2))
        z$V1 <- str_sub(z$V1,-1000,-2)
        Gistic_Bedfile <-
          cbind(x[, 1], y[, 1], z[, 1], gistic_region1[, 1])
        colnames(Gistic_Bedfile) <- c("chr", "start", "end", "peak")
        Gistic_Bedfile$chr <- as.character(Gistic_Bedfile$chr)
        Gistic_Bedfile$start <- as.integer(Gistic_Bedfile$start)
        Gistic_Bedfile$end <- as.integer(Gistic_Bedfile$end)
        Gistic_Bedfile <- na.omit(Gistic_Bedfile)
        Gistic_Bedfile <- data.frame(Gistic_Bedfile)
        region <- Gistic_Bedfile
        save(
          region,
          file = file.path(
            mirdirectory,
            "mirDriverFold",
            "miRDriver_Step1",
            "GisticResults",
            "gistic_bedfile.rda"
          )
        )
      }
    }
    if (gistic_bedfile == TRUE) {
      return(
        paste0(
          "gistic_bedfile.rda is created in mirdirectory under 'mirDriverFold/miRDriver_Step1/GisticResults' and regionwise gistic files are ready in folder 'miRDriver_Step1/Regionwise_Gistic_Files'"
        )
      )
    } else {
      return(list(
        paste0(
          "regionwise gistic files are ready in mirdirectory under 'mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files'"
        )
      ))
    }
  }
