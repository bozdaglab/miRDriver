## miRDriver STEP 1
## Code 1
## Created By: Banabithi Bose
## Date Created: 5/1/2019
getRegionWiseGistic1 <- function(humanGene,copyNumber,conf.lev,username,password,gistic_bedfile = c('TRUE', 'FALSE')){
    gisticfile <- NULL
    dir.create("~/mirDriverFold/miRDriver_Step1")
    dir.create("~/mirDriverFold/miRDriver_Step1/GisticResults")
    dir.create("~/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files")
    servername <- "https://cloud.genepattern.org/gp"
    gp.client <- gp.login(servername, username, password)
    file <-paste0("ftp://gpftp.broadinstitute.org/module_support_files/GISTIC/refgene/",humanGene)
    GISTIC_2.0.result <-run.analysis(gp.client,"urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00125:6.15.28",refgene.file = file,seg.file = copyNumber,confidence.level = conf.lev)
    ##Download the results from GISTIC 2.0
    setwd("~/")
    GISTIC_2.0.out.files <-job.result.download.files(GISTIC_2.0.result,"mirDriverFold/miRDriver_Step1/GisticResults/")
    ##Getting lesions file in miRDriver_Step1/GisticResults
    lesion_file <-list.files("~/mirDriverFold/miRDriver_Step1/GisticResults", pattern = "all_lesions")
    ##Processing region wise gistic files
    Gistic_lesions <-read.table(paste0("~/mirDriverFold/miRDriver_Step1/GisticResults/", lesion_file),sep = "\t",header = TRUE)
    Gistic_lesions <- data.table(Gistic_lesions)
    indx <-which(Gistic_lesions == '1' | Gistic_lesions == '2', arr.ind = TRUE)
    val <-tapply(names(Gistic_lesions)[indx[, 2]], indx[, 1], FUN = toString)
    Gistic_lesions$Aberrated_patients <-ifelse(seq_len(nrow(Gistic_lesions)) %in% names(val), val, NA)
    indx1 <- which(Gistic_lesions == '0', arr.ind = TRUE)
    val1 <-tapply(names(Gistic_lesions)[indx1[, 2]], indx1[, 1], FUN = toString)
    Gistic_lesions$NonAberrated_patients <-ifelse(seq_len(nrow(Gistic_lesions)) %in% names(val1), val1, NA)
    BR <-Gistic_lesions[, c(1,2,3,4,5,6,7,8,9, ncol(Gistic_lesions) - 1),(ncol(Gistic_lesions))]
    for (j in seq_len(nrow(BR))) {
        BR_Peak <- as.character(BR$Unique.Name[j])
        BR_Aberrated <-as.vector(strsplit(c(as.vector(BR$Aberrated[j])), ",")[[1]])
        BR_NonAberrated <-as.vector(strsplit(c(as.vector(BR$NonAberrated[j])), ",")[[1]])
        myfile <-file.path(paste0("~/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files/",BR_Peak,".txt"))
        fileConn <- file(myfile, open = "a")
        writeLines(paste(" ", ",Patient,", "Aberration"), fileConn)
        for (i in seq_along(BR_Aberrated))
            writeLines(paste(BR_Aberrated[i], ",1,", "aberrated"), fileConn)
        for (i in seq_along(BR_NonAberrated))
            writeLines(paste(BR_NonAberrated[i], ",1,", "non-aberrated"), fileConn)
        close(fileConn)
        gisticfile<-"~/mirDriverFold/miRDriver_Step1/GisticResults/all_lesions.conf_90.txt"
        if (gistic_bedfile == TRUE) {
            gistic_region <- read.table(gisticfile, sep = "\t", header = TRUE)
            gistic_region <-gistic_region[-grep("CN values", gistic_region$Unique.Name), ]
            gistic_region1 <- gistic_region[, c(1, 5)]
            x <- str_split_fixed(gistic_region1$Region.Limits, ":", 2)
            x <- data.table(x)
            y <- data.table(str_split_fixed(x$V2, "-", 2))
            z <- data.table(str_split_fixed(y$V2, "p", 2))
            z$V1 <- str_sub(z$V1, -1000, -2)
            Gistic_Bedfile <- cbind(x[, 1], y[, 1], z[, 1], gistic_region1[, 1])
            colnames(Gistic_Bedfile) <- c("chr", "start", "end", "peak")
            Gistic_Bedfile$chr <- as.character(Gistic_Bedfile$chr)
            Gistic_Bedfile$start <- as.integer(Gistic_Bedfile$start)
            Gistic_Bedfile$end <- as.integer(Gistic_Bedfile$end)
            Gistic_Bedfile <- na.omit(Gistic_Bedfile)
            Gistic_Bedfile <-data.frame(Gistic_Bedfile)
            region <- bedr.sort.region(Gistic_Bedfile)
            region$chr <- gsub("chr", "", region$chr)
            is.valid.region(region)
            save(region, file = "~/mirDriverFold/miRDriver_Step1/GisticResults/gistic_bedfile.rda")
        }
        if (gistic_bedfile == FALSE) {
            print("no gistic bedfile created")
        }
    }
    if (gistic_bedfile == TRUE) {
        return(paste0("gistic_bedfile.rda is created in folder 'miRDriver_Step1/GisticResults' and regionwise gistic files are ready in folder 'miRDriver_Step1/Regionwise_Gistic_Files'"))
    }
    if (gistic_bedfile == FALSE) {
        return(list(paste0("regionwise gistic files are ready in folder 'mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files'")))
    }
}
getRegionWiseGistic2 <-function(gisticfile,gistic_bedfile = c('TRUE', 'FALSE')) {
    dir.create("~/mirDriverFold")
    dir.create("~/mirDriverFold/miRDriver_Step1")
    dir.create("~/mirDriverFold/miRDriver_Step1/GisticResults")
    dir.create("~/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files")
    Gistic_lesions <- read.table(gisticfile, sep = "\t", header = TRUE)
    Gistic_lesions <-Gistic_lesions[-grep("CN values", Gistic_lesions$Unique.Name), ]
    Gistic_lesions <-  Gistic_lesions[seq_len(nrow(Gistic_lesions)), ]
    indx <-which(Gistic_lesions == '1' | Gistic_lesions == '2', arr.ind = TRUE)
    val <-tapply(names(Gistic_lesions)[indx[, 2]], indx[, 1], FUN = toString)
    Gistic_lesions$Aberrated_patients <-ifelse(seq_len(nrow(Gistic_lesions)) %in% names(val), val, NA)
    indx1 <- which(Gistic_lesions == '0', arr.ind = TRUE)
    val1 <-tapply(names(Gistic_lesions)[indx1[, 2]], indx1[, 1], FUN = toString)
    Gistic_lesions$NonAberrated_patients <-ifelse(seq_len(nrow(Gistic_lesions)) %in% names(val1), val1, NA)
    BR <-Gistic_lesions[, c(1,2,3,4,5,6,7,8,9, ncol(Gistic_lesions) - 1),(ncol(Gistic_lesions))]
    for (j in seq_len(nrow(BR))) {
        BR_Peak <- as.character(BR$Unique.Name[j])
        BR_Aberrated <-as.vector(strsplit(c(as.vector(BR$Aberrated[j])), ",")[[1]])
        BR_NonAberrated <-as.vector(strsplit(c(as.vector(BR$NonAberrated[j])), ",")[[1]])
        myfile <-file.path(paste0("~/mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files/",BR_Peak,".txt"))
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
            z$V1 <- str_sub(z$V1, -1000, -2)
            Gistic_Bedfile <- cbind(x[, 1], y[, 1], z[, 1], gistic_region1[, 1])
            colnames(Gistic_Bedfile) <- c("chr", "start", "end", "peak")
            Gistic_Bedfile$chr <- as.character(Gistic_Bedfile$chr)
            Gistic_Bedfile$start <- as.integer(Gistic_Bedfile$start)
            Gistic_Bedfile$end <- as.integer(Gistic_Bedfile$end)
            Gistic_Bedfile <- na.omit(Gistic_Bedfile)
            Gistic_Bedfile <-data.frame(Gistic_Bedfile)
            region <- bedr.sort.region(Gistic_Bedfile)
            is.valid.region(region)
            save(region, file = "~/mirDriverFold/miRDriver_Step1/GisticResults/gistic_bedfile.rda")
        }
        if (gistic_bedfile == FALSE) {
            print("no gistic bedfile created")
        }
    }
    if (gistic_bedfile == TRUE) {
        return(paste0("gistic_bedfile.rda is created in folder 'miRDriver_Step1/GisticResults' and regionwise gistic files are ready in folder 'miRDriver_Step1/Regionwise_Gistic_Files'"))
    }
    if (gistic_bedfile == FALSE) {
        return(list(paste0("regionwise gistic files are ready in folder 'mirDriverFold/miRDriver_Step1/Regionwise_Gistic_Files'")))
    }
}