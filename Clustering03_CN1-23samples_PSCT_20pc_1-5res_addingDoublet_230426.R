## 02Clustering ####
##Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.R
#!/bin/bash
#SBATCH -J  Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426
#SBATCH -D /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3 
#SBATCH -p compute
#SBATCH --output=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.out         
#SBATCH --error=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.err                  
#SBATCH --nodes=1                          
#SBATCH --ntasks-per-node=56 
#SBATCH --mem=500000    
#SBATCH --nodelist=n01


#**** 这个在MyeloidAnalysis05_CN1UMCN2_annoting_filting3_resfinding 文件后面

# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results


# ###mkdir /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/

# module load anaconda_global/anaconda_22
# source activate /share/home/sigrid/biosoft_sigrid/anaconda3/envs/R4
# R CMD BATCH --no-restore /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.R  /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.Rout
# conda deactivate
# module unload anaconda_global/anaconda_22 

# Error in paste0(filesDir, "Results/02Clustering/Figures/sceList2_MReRNA/",  : 
#   argument is missing, with no default
# Calls: res_persample -> dir.exists -> paste0
# Execution halted

#### prepare core&environment ####
rm(list=ls())
sink("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.txt", append=TRUE, split=TRUE )
pdf("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426/Clustering03_CN1-23samples_PSCT_20pc_1-5res_addingDoublet_230426.pdf") #
library(future)
plan("multicore", workers = 55) 
options(future.globals.maxSize= 900000*1024^2)
future.seed = NULL
nbrOfWorkers()


#### 9 samples ####

# save(sceList2_SCT1_1, file = paste0(filesDir, "Results/02Clustering/SaveData/sceList2_SCT1_1.RData"))
# save(sceList2_SCT1_3, file = paste0(filesDir, "Results/02Clustering/SaveData/sceList2_SCT1_3.RData"))
# save(sceList2_SCT1_3, file = paste0(filesDir, "Results/02Clustering/SaveData/sceList2_SCT1_3.RData"))

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



    #load(paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_adddoublet_persample_DF.RData"))
    #CN1_23filter_adddoublet_persample_DF

        

#### extract doubletfinder parameters and create doublet excel ####
    #先在一个dataframe中添加信息
    sceList1_DF_df <- as.data.frame(matrix(nrow=0,ncol=16)) 
    sceList1_DF_df
    ncol(sceList1_DF_df)
    colname <- c("samplename", "pos_QC1", "manpc", "manres",
                              "pc20.res1_Doublet", "pc20.res1_Singlet", "pc20.res1_doublet_rate", "pc20.res1_Singlet_rate", 
                              "pc20.res5_Doublet", "pc20.res5_Singlet", "pc20.res5_doublet_rate", "pc20.res5_Singlet_rate",
                              "Doublet_Doublet", "Doublet_Singlet", "Doublet_doublet_rate", "Doublet_Singlet_rate")
    names(sceList1_DF_df) <- colname
    length(colname)



    # sceList1_DF_df <- data.frame()


    # colnames(sceList1_DF_df) <- c("samplename", "pos_QC1", "manpc", "manres",
    #                               "pc20.res1_Doublet", "pc20.res1_Singlet", "pc20.res1_doublet_rate", "pc20.res1_Singlet_rate", 
    #                               "pc20.res5_Doublet", "pc20.res5_Singlet", "pc20.res5_doublet_rate", "pc20.res5_Singlet_rate",
    #                               "Doublet_Doublet", "Doublet_Singlet", "Doublet_doublet_rate", "Doublet_Singlet_rate")
    #### load pc&res_df  ####
        # npc_res_persample_xlsx<-readxl::read_xlsx(paste0(filesDir, "Results/02Clustering/pre_annotation_persample/sceList1_SCT_npc_res_persample.xlsx")) 
        # npc_res_persample_xlsx
        # npc_res_persample_df<-as.data.frame(npc_res_persample_xlsx)
        # rm(npc_res_persample_xlsx)
        # npc_res_persample_df
        namelist <- names(CN1_23filter_adddoublet_persample_DF)
    #### extract doubletfinder parameters and create doublet excel ####
    for (samplename in namelist) {
        print(samplename)

        sc_persampledata <- CN1_23filter_adddoublet_persample_DF[[samplename]]
        colnames(sc_persampledata@meta.data)
        print(paste0(samplename, "_colnames:_", colnames(sc_persampledata@meta.data)))
        manpc <- 20
        manres <- 0.1

        # manpc <- npc_res_persample_df[npc_res_persample_df$samplename == samplename,"npc_num"]
        # manres <- npc_res_persample_df[npc_res_persample_df$samplename == samplename,"res_num"]

        #### pANN_0.25_0.29_418 manpcres ####

            table(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res) 
            pc20.res1_DoubletConfidience_df <- as.data.frame(table(sc_persampledata@meta.data[ ,"DF_hi_lo_20pc_0.1res"]))
            pc20.res1_DoubletConfidience_df
            pc20.res1_Doublet <- as.character(pc20.res1_DoubletConfidience_df[pc20.res1_DoubletConfidience_df$Var1=="Doublet_High_Confidience",2]) 
            pc20.res1_Doublet
            pc20.res1_Singlet <- as.character(pc20.res1_DoubletConfidience_df[pc20.res1_DoubletConfidience_df$Var1=="Doublet_Low_Zero_Confidience",2])
            pc20.res1_Singlet
            pc20.res1_doublet_rate <- round(as.numeric(pc20.res1_Doublet)/ncol(sc_persampledata),2)
            pc20.res1_Singlet_rate <- round(as.numeric(pc20.res1_Singlet)/ncol(sc_persampledata),2)
            rm(pc20.res1_DoubletConfidience_df)


        #### pANN_0.25_0.29_418  mpc40.res5 ####
            table(sc_persampledata@meta.data$DF_hi_lo_20pc_5res) 
            pc20.res5_DoubletConfidience_df <- as.data.frame(table(sc_persampledata@meta.data[ ,"DF_hi_lo_20pc_5res"]))
            pc20.res5_DoubletConfidience_df
            pc20.res5_Doublet <- as.character(pc20.res5_DoubletConfidience_df[pc20.res5_DoubletConfidience_df$Var1=="Doublet_High_Confidience",2]) 
            pc20.res5_Doublet
            pc20.res5_Singlet <- as.character(pc20.res5_DoubletConfidience_df[pc20.res5_DoubletConfidience_df$Var1=="Doublet_Low_Zero_Confidience",2])
            pc20.res5_Singlet
            pc20.res5_doublet_rate <- round(as.numeric(pc20.res5_Doublet)/ncol(sc_persampledata),2)
            pc20.res5_Singlet_rate <- round(as.numeric(pc20.res5_Singlet)/ncol(sc_persampledata),2)
            rm(pc20.res5_DoubletConfidience_df)

        #### pANN_0.25_0.29_418  Doublet ####
            ncol(sc_persampledata)
            table(sc_persampledata@meta.data$Doublet) 
            Doublet_df <- as.data.frame(table(sc_persampledata@meta.data[ ,"Doublet"]))
            Doublet_df
            Doublet_Doublet <- as.character(Doublet_df[Doublet_df$Var1=="Doublet_High_Confidience",2]) 
            Doublet_Doublet
            Doublet_Singlet <- as.character(Doublet_df[Doublet_df$Var1=="Doublet_Zero_Confidience",2])
            Doublet_Singlet
            Doublet_doublet_rate <- round(as.numeric(Doublet_Doublet)/ncol(sc_persampledata),2)
            Doublet_Singlet_rate <- round(as.numeric(Doublet_Singlet)/ncol(sc_persampledata),2)


        #### fill information  ####
            sampleinfo <- c( samplename, as.character(ncol(sc_persampledata)) , as.character(manpc), as.character(manres), 
                                    as.character(pc20.res1_Doublet), as.character(pc20.res1_Singlet), as.character(pc20.res1_doublet_rate), as.character(pc20.res1_Singlet_rate), 
                                    as.character(pc20.res5_Doublet), as.character(pc20.res5_Singlet), as.character(pc20.res5_doublet_rate), as.character(pc20.res5_Singlet_rate),
                                    as.character(Doublet_Doublet), as.character(Doublet_Singlet), as.character(Doublet_doublet_rate), as.character(Doublet_Singlet_rate) )
            length(sampleinfo)
            sceList1_DF_df <- rbind(sceList1_DF_df, sampleinfo)
    }

    names(sceList1_DF_df) <- colname

    write.csv(sceList1_DF_df, file = paste0(filesDir, "Results/02Clustering/Tables/sceList23samples_DF_count_230426.csv"), quote = FALSE)
    

rm(list = ls())





















# ####  meta.data cleanning ####



#     #load(paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_adddoublet_persample_DF.RData"))
#     #CN1_23filter_pure_persample_DF
#     namelist <- names(CN1_23filter_adddoublet_persample_DF)

# ####  persample meta.data cleansing ####
# CN1_23filter_pure_persample_DF <- lapply(namelist, function(samplename){
#     sc_persampledata <- CN1_23filter_adddoublet_persample_DF[[samplename]]
#     print(samplename)
#     #### remove pANN_columname ####
#     columnname <- colnames(sc_persampledata@meta.data)
#     class(columnname)
#     pANN_columname <- grep ("^pANN", columnname, value = FALSE)
#     pANN_columname
#     sc_persampledata@meta.data <- sc_persampledata@meta.data[,- pANN_columname]

#     ####  remove DF.classifications_columname ####
#     columnname <- colnames(sc_persampledata@meta.data)
#     DF.classifications_columname <- grep("^DF.classification", columnname, value = FALSE)
#     DF.classifications_columname
#     sc_persampledata@meta.data <- sc_persampledata@meta.data[,- DF.classifications_columname]
#     colnames(sc_persampledata@meta.data)

#     ####  remove SCT_snn_res_columname ####
#     columnname <- colnames(sc_persampledata@meta.data)
#     SCT_snn_res_columname <- grep("^SCT_snn_res", columnname, value = FALSE)
#     SCT_snn_res_columname
#     sc_persampledata@meta.data <- sc_persampledata@meta.data[,- SCT_snn_res_columname]
#     colnames(sc_persampledata@meta.data)

#     columnname <- colnames(sc_persampledata@meta.data)
#     print(paste0(samplename, "columnname: ", columnname))
#     columnname 
#     return(sc_persampledata)
#     })

# names(CN1_23filter_pure_persample_DF)
# names(CN1_23filter_pure_persample_DF)<- namelist
# names(CN1_23filter_pure_persample_DF)
# save(CN1_23filter_pure_persample_DF, file = paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_pure_persample_DF.RData"))
# rm(list = ls())





    #rm(list = ls())









 rm(list = ls())    
