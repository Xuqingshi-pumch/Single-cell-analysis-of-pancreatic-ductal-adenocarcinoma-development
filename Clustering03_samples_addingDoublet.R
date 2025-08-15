#### prepare core&environment ####
rm(list=ls())
sink("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.txt", append=TRUE, split=TRUE )
pdf("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.pdf") #
library(future)
plan("multicore", workers = 55) 
options(future.globals.maxSize= 900000*1024^2)
future.seed = NULL
nbrOfWorkers()

#### prepare library ####
    suppressPackageStartupMessages(library(Seurat))
    suppressPackageStartupMessages(library(future))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(clustree)) 
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(dplyr)) 
    suppressPackageStartupMessages(library(data.table)) 
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(Matrix))
    suppressPackageStartupMessages(library(patchwork))
    suppressPackageStartupMessages(library(sctransform))
    library(viridis)
    library(maps)
    library(fields)
    library(spam)
    suppressPackageStartupMessages(library(DoubletFinder))
    suppressPackageStartupMessages(library(harmony))
    #suppressPackageStartupMessages(library(clusterProfiler))
    #suppressPackageStartupMessages(library(org.Hs.eg.db))
    options(stringsAsFactors = F)
    filesDir<-'/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/'
    setwd(filesDir)
    getwd()

    if(!dir.exists(paste0(filesDir, "Results/02Clustering/SaveData/"))){dir.create(paste0(filesDir, "Results/02Clustering/SaveData/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Figures/"))){dir.create(paste0(filesDir, "Results/02Clustering/Figures/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Tables/"))){dir.create(paste0(filesDir, "Results/02Clustering/Tables/"), recursive = TRUE)}


    CN1_23filter_PSCT_20pc.1res_DF <- readRDS(paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_PSCT_bothDF_230426.rds"))



    namelist <- names(CN1_23filter_PSCT_20pc.1res_DF)

#### rearrange metadata column before merging ####
    CN1_23filter_adddoublet_persample_DF <- lapply(namelist,function(samplename){
        samplename
        sc_persampledata <- CN1_23filter_PSCT_20pc.1res_DF[[samplename]]
        #### add "doublet" column before merging ####
        # DF_hi_lo                  DF_hi_lo_20pc_5res
        # Doublet_High_Confidience  Doublet_High_Confidience
        # Doublet_Low_Confidience       Doublet_Low_Zero_Confidience
        # Doublet_Zero_Confidience
        ##DF没有修改        
        sc_persampledata@meta.data[,"Doublet"] <- sc_persampledata@meta.data[,"DF_hi_lo_20pc_0.1res"]
            sc_persampledata@meta.data$Doublet[which(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res == "Doublet_Low_Zero_Confidience" & 
                                                  sc_persampledata@meta.data[, "DF_hi_lo_20pc_5res"] == "Doublet_Low_Zero_Confidience")] <- "Doublet_Zero_Confidience"


            sc_persampledata@meta.data$Doublet[which(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res == "Doublet_Low_Zero_Confidience" & 
                                                  sc_persampledata@meta.data[, "DF_hi_lo_20pc_5res"] == "Doublet_High_Confidience")] <- "Doublet_Low_Confidience" #Low
            
            sc_persampledata@meta.data$Doublet[which(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res == "Doublet_High_Confidience" & 
                                                  sc_persampledata@meta.data[, "DF_hi_lo_20pc_5res"] == "Doublet_High_Confidience")] <- "Doublet_High_Confidience"

            sc_persampledata@meta.data$Doublet[which(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res == "Doublet_High_Confidience" & 
                                                  sc_persampledata@meta.data[, "DF_hi_lo_20pc_5res"] == "Doublet_Low_Zero_Confidience")] <- "Doublet_Low_Confidience"
            
        head(sc_persampledata)
        return(sc_persampledata)
    })
    
    names(CN1_23filter_adddoublet_persample_DF)
    names(CN1_23filter_adddoublet_persample_DF) <- namelist

    names(CN1_23filter_adddoublet_persample_DF)
    saveRDS(CN1_23filter_adddoublet_persample_DF,paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_adddoublet_persample_DF_230426.rds"))


rm(list = ls())

