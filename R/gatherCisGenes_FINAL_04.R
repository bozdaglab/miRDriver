## miRDriver STEP 3
## Code 1
## Created By: Banabithi Bose
## Date Created: 5/1/2019

gatherCisGenes <-function(ncores = 2,mirna_bedfile,gistic_bedfile,parallel = c('TRUE', 'FALSE'),mirdirectory="~/") {



    if ((is(mirna_bedfile, "data.frame")) || (is(mirna_bedfile,"matrix"))){
        mirna_bedfile<-mirna_bedfile
    }else if (is(mirna_bedfile,"SummarizedExperiment")){
        mirna_bedfile<-assay(mirna_bedfile)
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

        tc_de<-NULL

        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3"))
        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene"))
        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS"))
        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion"))
        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files"))
        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/GMIRNA"))
        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_mirna"))
        dir.create(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_DEcis"))
        processGeneFiles <- function(fileName) {
                cis_genes <- NULL
                tryCatch({
                        load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/trans_genes.Rda"))
                        load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/",fileName,"/cis_genes.Rda"))
                        trans_genes <- data.table(trans_genes)
                        trans_genes <- trans_genes[trans_genes$gene_id != "", ]
                        if (nrow(cis_genes) == 0) {
                                #print("this region does not have any cis genes")
                        } else{
                                for (i in seq_along(trans_genes$gene_id))
                                {
                                        tc_de <- cbind(trans_genes[i], t(cis_genes))
                                        fname <- paste0(trans_genes$gene_id[i], "_", fileName)
                                        ## save(tc_de,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/",fname,".Rda"))
                                        save(tc_de,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/",fname,".Rda"))
                                        
                                }
                        }
                        return()
                }, error = function(error_condition) {
                        cat(paste0(fileName, " : ", error_condition),file = paste0(mirdirectory,"/mirDriverFold/ERROR_FILES/Step3_TransGene_w_all_cisgenes_from_all_regions.txt"),sep = "\n",append = TRUE)
                        return()
                }, finally = {
                })
        }
        if (parallel == TRUE) {
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
                invisible(clusterEvalQ(cl, library(pbapply)))
                clusterExport(cl = cl, "processGeneFiles", envir = environment())
                ###### MAIN PROGRAM############################################################
                fileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/"))
                parLapply(cl, fileNames, processGeneFiles)
                stopCluster(cl)
        }
        if (parallel == FALSE) {
                fileNames <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step2/Trans_Cis_Files/"))
                lapply(fileNames, processGeneFiles)
        }
        #################### END ##############################################
        ###### CODE 2
        ## Putting same gene files from different regions in a common folder DE_GCIS for each trans gene
        ## First Create the gene folders
        ###CODE START###
        #print("second code starts")
        fileList <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/"),pattern = "*.Rda")
        files <- data.table(fileList)
        files <- data.table(str_split_fixed(files$fileList,"_", 2))[,1,drop=FALSE]
        files <- unique(files)
        files <- data.table(unique(files))
        for (i in seq_along(files$V1)) {
                newdir <- paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS/", files$V1[i])
                dir.create(newdir)
        }
        #################
        current_folder <-paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene")
        fileList <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/"),pattern = "*.Rda")
        for (i in seq_along(fileList)) {
                newlocation = unique(data.table(str_split_fixed(fileList[i],"_", 2))[,1])$V1
                new_folder <-paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS/", newlocation)
                
                load(paste0(mirdirectory, "/mirDriverFold/miRDriver_Step3/transgene_w_DEcisgene/",fileList[i]))
                save(tc_de, file=paste0(new_folder, "/", fileList[i]))
                #defaultW <- getOption("warn") 
                #options(warn = -1)
                
                #file.move(file.path(current_folder, fileList[i]), new_folder, overwrite = FALSE)
                
                #options(warn = defaultW)
        }
        ##The file numbers in each newlocation says for each gene how many gistic regions are there.
        ### CODE END ####
        ##Code 3
        #print("code 3 start")
        processGM <- function(GeneFile) {
                tc_de <- NULL
                tryCatch({
                        file_list1 <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS/", GeneFile))
                        for (file in file_list1) {
                                load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS/",GeneFile,"/",file))
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
                        save(TC,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_DEcis/",GeneFile,".Rda"))
                }, error = function(error_condition) {
                        cat(paste0(GeneFile, " : ", error_condition),file = paste0(mirdirectory,"/mirDriverFold/ERROR_FILES/Step3_TransGene_w_all_cisgenes_from_all_regions.txt"),sep = "\n",append = TRUE)
                        return()
                }, finally = {
                })
        }
        if (parallel == TRUE) {
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
                invisible(clusterEvalQ(cl, library(pbapply)))
                clusterExport(cl = cl, "processGM", envir = environment())
                ### MAIN FUNCTION###
                GeneFileList <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS"))
                parLapply(cl, GeneFileList, processGM)
                stopCluster(cl)
        }
        if (parallel == FALSE) {
                GeneFileList <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS"))
                lapply(GeneFileList, processGM)
        }
        GeneFileList1 <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_DEcis"))
        ##Code 4
        ##This code is for merging TCGA mirna expression data with Gistic regions to get miRNAs
        ##Using Gistic region bed file
        #print("code 4 starts")
        if ((gistic_bedfile == FALSE) [1,1]) {
                load(paste0(mirdirectory,"/mirDriverFold/EXAMPLE_DATASETS/gistic_bedfile.rda"))
        } else{
                region <- gistic_bedfile
        }
        #new bedr replacement
        region$strand<-"+"
        region<-region[,c(1,5,2,3,4)]
        region$chr<-paste0("chr", region$chr)
        colnames(region)<-c("seqnames", "strand", "start", "end", "peak")
        query<-region%>%as_granges()
        
        
        mirna_bedfile$strand<-"+"
        mirna_bedfile<-mirna_bedfile[,c(1,5,2,3,4)]
        mirna_bedfile$chr<-paste0("chr", mirna_bedfile$chr)
        colnames(mirna_bedfile)<-c("seqnames", "strand", "start", "end", "mirna")
        subject<-mirna_bedfile%>%as_granges()
        intersect_rng <- join_overlap_intersect(query, subject)
        overlapping_pairs <-as.data.frame(intersect_rng)
        #end
        
        overlapping_pairs<-overlapping_pairs[,c(6,7)]
        colnames(overlapping_pairs)<-c("region.PEAK", "mirna")
        
        z <-data.table(str_split_fixed(overlapping_pairs$region.PEAK, ",", 2))
        overlapping_pairs$region.PEAK <- z$V1
        colnames(overlapping_pairs)<-c("regionPEAK", "miRNA")
        
        ########  making files with mirnas for each gistic regions
        mirna_in_gistic_region <-unique(sqldf("SELECT regionPEAK, miRNA from overlapping_pairs order by 1"))
        for (i in seq_len(nrow(mirna_in_gistic_region)))
        {
                myfile <-file.path(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion/",mirna_in_gistic_region[i, 1],".txt"))
                fileConn <- file(myfile, open = "a")
                writeLines(paste(mirna_in_gistic_region[i, 2]), fileConn)
                close(fileConn)
        }
        tcga_mirnas_count_by_gistic_region <-sqldf("select regionPEAK, count(regionPEAK) from mirna_in_gistic_region group by regionPEAK")
        save(tcga_mirnas_count_by_gistic_region, file = paste0(mirdirectory,"/mirDriverFold/miRDriver_RESULTS/mirnas_count_by_gistic_region.Rda"))
        ###################Reading region with miRNAs list files##########
        ## This will use miRNA_by_GisticRegion and EDGER_DIRECTORIES/up_down_gene files to create gene with all miRNA files.
        ## CODE 5:
        #print("code5 starts")
        fileNamesMiRNA <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion/"),pattern = "*.txt")
        processUpDownGeneFiles <- function(fileNames) {
                f<- NULL
                file_list1 <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS/", fileNames))
                gene <- fileNames
        for (i in seq_along(file_list1)) {
                mirnafilename <- str_split_fixed(file_list1[i], "_", 2)[,2]
                mirnafilename <- str_replace(mirnafilename, ".Rda", ".txt")
                if (is.element(mirnafilename, fileNamesMiRNA))
                {
                        RegionMiRna <-read.csv(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/miRNA_by_GisticRegion/",mirnafilename,sep = ""),sep = ";",header = FALSE)
                        Lasso_RegionMiRna <- data.table(RegionMiRna)
                        colnames(Lasso_RegionMiRna) <- "miRNA_ID"
                        gm <- data.table(cbind(gene, t(Lasso_RegionMiRna)))
                        fname <- paste0(gene, "_", mirnafilename)
                        fname <- str_replace(fname, ".txt","")
                        save( gm,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/",fname,".Rda")) 
        
                }else {
                        #print(paste0(mirnafilename, " this region has no mirna for ", gene))
                }
            }
        }
        ##### MAIN PROGRAM############################################################
        ######################
        fileList <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/DE_GCIS/"))
        files <- unique(data.table(fileList))
        lapply(files$fileList, processUpDownGeneFiles) 
        ################### END ##############################################
        ###### CODE 6
        ##The file numbers in each newlocation says for each gene how many gistic regions are there.
        ## Putting same gene files from different regions in a common folder
        ## First Create the gene folders
        ## CODE START###
        fileList <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/"),pattern = "*.Rda")
        files <- data.table(fileList)
        files <- data.table(str_split_fixed(files$fileList, "_", 2))[,1,drop=FALSE]
        files <- unique(files)
        for (i in seq_along(files$V1)) {
                newdir <- paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/GMIRNA/", files$V1[i])
                dir.create(newdir)
        }
        ## Putting files in gene folders
        for (i in seq_along(fileList)) {
                newlocation = unique(data.table(str_split_fixed(fileList[i],"_", 2))[,1])$V1
                current_folder <-paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/Trans_Mirna_Files/")
                new_folder <-paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/GMIRNA/", newlocation)
                file.copy(file.path(current_folder, fileList[i]), new_folder)
        }
        ### CODE END ####
        ###CODE 7
        ###Merging all files in PIPELINE2_OV_LASSO/GMIRNA , each gene folder Parallel Code
        #print("code 7 start")
        processGM <- function(GeneFile) {
                tryCatch({
                        file_list1 <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/GMIRNA/", GeneFile))
                        for (file in file_list1) {
                                load(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/GMIRNA/",GeneFile,"/",file))
                                gm <- data.table(t(gm))
                                colnames(gm) <- as.character(unlist(gm[1, ]))
                                gm <- gm[-1, ]
                                if (!exists("GM"))
                                {
                                GM <- gm
                                }else
                                {
                                temp_GM <- gm
                                GM <- rbind(GM, temp_GM)
                                rm(temp_GM)
                                }
                        }
                        GM <- unique(GM)
                        save(GM,file = paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_mirna/",GeneFile,".Rda"))
                }, error = function(error_condition) {
                        cat(paste0(GeneFile, " : ", error_condition),file = paste0(mirdirectory,"/mirDriverFold/ERROR_FILES/Step3_TransGene_w_all_mirnagenes_from_all_regions.txt"),sep = "\n",append = TRUE)
                        return()
                }, finally = {
                })
        }
        if (parallel == TRUE) {
                if(.Platform$OS.type == "windows"){
                    cl <- makeCluster(mc <- getOption("cl.cores", ncores), type = "PSOCK")
                }else{
                    cl <- makeCluster(mc <- getOption("cl.cores", ncores), type = "FORK")
                }
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
                ########### MAIN FUNCTION######################
                GeneFileList <-list.files(paste0(mirdirectory,"mirDriverFold/miRDriver_Step3/GMIRNA"))
                parLapply(cl, GeneFileList, processGM)
                stopCluster(cl)
        }
        if (parallel == FALSE) {
                GeneFileList <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/GMIRNA"))
                lapply(GeneFileList, processGM)
        }
        ###End
        GeneFileList1 <-list.files(paste0(mirdirectory,"/mirDriverFold/miRDriver_Step3/gene_all_mirna/"))
        return(paste0("mirDriverFold/miRDriver_Step3/gene_all_mirna and gene_all_DEcis are filled with cis genes and cis mirnas"))
        return(paste0("gene_all_mirna and gene_all_DEcis folders are filled in mirDriverFold/miRDriver_Step3"))
        #options(warn = defaultW)
}


