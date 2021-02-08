## miRDriver STEP 2
## Code 2
## Created By: Banabithi Bose
## Date Created: 5/1/2019
makingCisAndTransGenes <-function(Genes_coord_bedfile, gistic_bedfile, mirdirectory="~/") {
    if ((is(Genes_coord_bedfile, "data.frame")) || (is(Genes_coord_bedfile,"matrix"))){
        Genes_coord_bedfile<-Genes_coord_bedfile
    }else if (is(Genes_coord_bedfile,"SummarizedExperiment")){
        Genes_coord_bedfile<-assay(Genes_coord_bedfile)
    }else{
        #print("")
    }

    if ((is(gistic_bedfile, "data.frame")) || (is(gistic_bedfile,"matrix"))){
        gistic_bedfile<-gistic_bedfile
    }else if (is(gistic_bedfile,"SummarizedExperiment")){
        gistic_bedfile<-assay(gistic_bedfile)
    }else{
        #print("")
    }


    ##making necessary directories
    dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_RESULTS"))
    dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES"))
    dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files"))
    ############## Reading up/down genes files for a gistic region
    fileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/UpDownGeneLibrary/"),pattern = "*.txt")
    Count_Cis_Trans <- data.frame("", "", "")
    colnames(Count_Cis_Trans) <-c("Gistic_Region", "cis_genes_count", "trans_genes_count")
    i <- as.integer(seq_along(fileNames))
    MakecistransFile <- function(i) {
        #print(i)
        upgene <-read.table(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/UpDownGeneLibrary/",fileNames[i]),sep = " ",header = TRUE)
        genes <- data.frame(upgene$Row.Names)
        colnames(genes) <- "gene_id"
        Gistic_Region1 <- gsub(".*_", "", fileNames[i])
        Gistic_Region <-str_sub(Gistic_Region1, 1, str_length(Gistic_Region1) - 4)
        if ((gistic_bedfile == FALSE)[1] ){
            load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step1/GisticResults/gistic_bedfile.rda"))
            region$chr <- gsub("[chr]", "", region$chr)
            Region = region[region$peak == Gistic_Region, ]
        } else{
            Region = gistic_bedfile[gistic_bedfile$peak == Gistic_Region, ]
        }
        if (nrow(Region) == 0) {
            #print("this gistic region is not valid")
        } else {
            ##New join code
            Region$strand<-"+"
            Region<-Region[,c(1,5,2,3,4)]
            Region$chr<-paste0("chr", Region$chr)
            colnames(Region)<-c("seqnames", "strand", "start", "end", "peak")
            query<-Region%>%as_granges()
            
            Genes_coord_bedfile$strand<-"+"
            Genes_coord_bedfile<-Genes_coord_bedfile[,c(1,5,2,3,4)]
            Genes_coord_bedfile$chr<-paste0("chr", Genes_coord_bedfile$chr)
            colnames(Genes_coord_bedfile)<-c("seqnames", "strand", "start", "end", "gene")
            subject<-Genes_coord_bedfile%>%as_granges()
            intersect_rng <- suppressWarnings(join_overlap_intersect(query, subject))
            Gistic_Genes<-as.data.frame(intersect_rng)
            ##end new join code
            
            cis_genes <-data.frame(Gistic_Genes$gene)
            if (nrow(cis_genes) == 0) {
                #print("This region has no cis genes")
                trans_genes <- genes
                colnames(trans_genes) <- "gene_id"
                save(trans_genes,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/trans_",fileNames[i],".Rda"))
            } else{
                colnames(cis_genes) <- "gene_id"
                save(cis_genes,file = paste0("~/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/cis_",fileNames[i],".Rda"))
                trans_genes <- data.frame(genes[-c(cis_genes$gene_id), ])
                colnames(trans_genes) <- "gene_id"
                save(trans_genes,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/trans_",fileNames[i],".Rda"))
            }
        Count_Cis_Trans1 <-data.table(cbind(fileNames[i], nrow(unique(cis_genes)), nrow(unique(trans_genes))))
        colnames(Count_Cis_Trans1) <-c("Gistic_Region","cis_genes_count","trans_genes_count")
        return(Count_Cis_Trans1)
        }
    }   
    Count_Cis_Trans <- rbindlist(lapply(i, MakecistransFile))
    Count_Cis_Trans <- Count_Cis_Trans[-1, ]
    Count_Cis_Trans <- unique(Count_Cis_Trans)
    save(Count_Cis_Trans, file = paste0(mirdirectory,"/mirDriverFold/miRDriver_RESULTS/Count_Cis_Trans.Rda"))
    r <- Count_Cis_Trans[Count_Cis_Trans$trans_genes_count == 0, ]
    #### Removing regions from CIS_TRANS_FILES where there are no trans genes
    fileNames1 <-list.files("~/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/", pattern = "*.Rda")
    fileNames2 <- r$Gistic_Region
    for (i in seq_along(fileNames1)){
        fileNamesD <- gsub("(cis_)", "", fileNames1[i])
        fileNamesD <- gsub("(trans_)", "", fileNamesD)
        fileNamesD <- str_sub(fileNamesD, end = -5)
        if (is.element(fileNamesD, fileNames2)) {
            file.remove(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/",fileNames1[i]))
        } else
        {
        }
    }
    ##Merging script portion of the function
    processUpDownGeneFiles <- function(fileNameUpDn) {
        trans_genes <- NULL
        RegionNameOnly <-str_sub(gsub(".*_", "", fileNameUpDn), length(fileNameUpDn), -9)
        if (str_sub(fileNameUpDn, 1, 3) == "cis") {
            load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/",fileNameUpDn))
            load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/",gsub(".*cis", "trans", fileNameUpDn)))
        } else if (str_sub(fileNameUpDn, 1, 5) == "trans")
        {
            load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/",fileNameUpDn))
            cis_file_name <- paste0(gsub(".*trans", "cis", fileNameUpDn))
            if (is.element(cis_file_name, fileNames))
            {
                load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES/",gsub(".*trans", "cis", fileNameUpDn)))
            } else
            {
                cis_genes <- data.frame("")
                colnames(cis_genes) <- "gene_id"
            }
        }
        transfile <- str_sub(fileNameUpDn, -80, -9)
        transfile <- gsub("cis_", "trans_", transfile)
        if(!dir.exists(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",RegionNameOnly))){
            dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",RegionNameOnly))
        }
        save(cis_genes,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",RegionNameOnly,"/cis_genes.Rda"))
        save(trans_genes,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",RegionNameOnly,"/",transfile,".Rda"))
    }
    ##### MAIN PROGRAM ############################################################
    fileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/CIS_TRANS_FILES"), pattern = "*.Rda")
    lapply(fileNames, processUpDownGeneFiles)
    ################### END #######################################################
    ## MERGING trans_UPgene and trans_DOWNgene data
    processUpDownGeneFiles2 <- function(fileName) {
        trans_gene <- NULL
        genefiles <- NULL
        genefiles<-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/"))
        transfiles <- paste0("trans_DOWNgene_", fileName, ".Rda")
        transfiles_new <- paste0("trans_UPgene_", fileName, ".Rda")
        if (is.element(transfiles, genefiles)) {
            load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/",transfiles))
            transgene1 <- trans_genes
            colnames(transgene1) <- "gene_id"
            if (is.element(transfiles_new, genefiles)) {
                load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/",transfiles_new))
                transgene2 <- trans_genes
                colnames(transgene2) <- "gene_id"
            } else {
                transgene2 <- data.frame("")
                colnames(transgene2) <- "gene_id"
            }
            trans_genes <- rbind(transgene1, transgene2)
            save(trans_genes,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/trans_genes.Rda"))
        } else{
            load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/",transfiles_new))
            transgene2 <- trans_genes
            transgene1 <- data.frame("")
            colnames(transgene2) <- "gene_id"
            colnames(transgene1) <- "gene_id"
            trans_genes <- rbind(transgene1, transgene2)
            save(trans_genes,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/trans_genes.Rda"))
        }
    }
    folderNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/"))
    lapply(folderNames, processUpDownGeneFiles2)
    return(paste0("Trans genes and cis genes files produced in folder mirDriverFold/miRDriver_Step2/Trans_Cis_Files"))
}
