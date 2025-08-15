## 02Clustering ####
##Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426.R
#!/bin/bash
#SBATCH -J  Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426
#SBATCH -D /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3 
#SBATCH -p compute
#SBATCH --output=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426.out         
#SBATCH --error=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426.err                  
#SBATCH --nodes=1                          
#SBATCH --ntasks-per-node=56 
#SBATCH --mem=500000    
#SBATCH --nodelist=n01


#**** 这个在MyeloidAnalysis05_CN1UMCN2_annoting_filting3_resfinding 文件后面

# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results


# ###mkdir /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426/

# module load anaconda_global/anaconda_22
# source activate /share/home/sigrid/biosoft_sigrid/anaconda3/envs/R4
# R CMD BATCH --no-restore /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles/02Clustering/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426.R  /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426.Rout
# conda deactivate
# module unload anaconda_global/anaconda_22 

# Error in paste0(filesDir, "Results/02Clustering/Figures/sceList2_MReRNA/",  : 
#   argument is missing, with no default
# Calls: res_persample -> dir.exists -> paste0
# Execution halted

#### prepare core&environment ####
rm(list=ls())
sink("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426.txt", append=TRUE, split=TRUE )
pdf("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426/Clustering04_CN1-23samples_PSCT_puring_rmDoublet_CN1_23samplesUMCN2_Merging_230426.pdf") #
library(future)
plan("multicore", workers = 55) 
options(future.globals.maxSize= 9000000*1024^2)
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



#### prepare folderdirectory ####
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/SaveData/"))){dir.create(paste0(filesDir, "Results/02Clustering/SaveData/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Figures/"))){dir.create(paste0(filesDir, "Results/02Clustering/Figures/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Tables/"))){dir.create(paste0(filesDir, "Results/02Clustering/Tables/"), recursive = TRUE)}

    if(!dir.exists(paste0(filesDir, "Results/03Integration/SaveData/"))){dir.create(paste0(filesDir, "Results/03Integration/SaveData/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/03Integration/Figures/"))){dir.create(paste0(filesDir, "Results/03Integration/Figures/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/03Integration/Tables/"))){dir.create(paste0(filesDir, "Results/03Integration/Tables/"), recursive = TRUE)}

     
#第一遍跑过    CN1_23filter_adddoublet_persample_DF <- readRDS(paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_adddoublet_persample_DF_230426.rds"))


#第一遍跑过####  meta.data cleanning ####
#第一遍跑过    namelist <- names(CN1_23filter_adddoublet_persample_DF)
#第一遍跑过    print(namelist)
#第一遍跑过####  persample meta.data cleansing ####
#第一遍跑过    CN1_23filter_PSCT_pure_persample_DF <- lapply(namelist, function(samplename){
#第一遍跑过    sc_persampledata <- CN1_23filter_adddoublet_persample_DF[[samplename]]
#第一遍跑过    print(samplename)
#第一遍跑过    #### remove pANN_columname ####
#第一遍跑过    columnname <- colnames(sc_persampledata@meta.data)
#第一遍跑过    class(columnname)
#第一遍跑过    pANN_columname <- grep ("^pANN", columnname, value = FALSE)
#第一遍跑过    pANN_columname
#第一遍跑过    sc_persampledata@meta.data <- sc_persampledata@meta.data[,- pANN_columname]
#第一遍跑过
#第一遍跑过    ####  remove DF.classifications_columname ####
#第一遍跑过    columnname <- colnames(sc_persampledata@meta.data)
#第一遍跑过    DF.classifications_columname <- grep("^DF.classification", columnname, value = FALSE)
#第一遍跑过    DF.classifications_columname
#第一遍跑过    sc_persampledata@meta.data <- sc_persampledata@meta.data[,- DF.classifications_columname]
#第一遍跑过    colnames(sc_persampledata@meta.data)
#第一遍跑过
#第一遍跑过    columnname <- colnames(sc_persampledata@meta.data)
#第一遍跑过    DF.classifications_columname <- grep("^DF_hi_lo_20pc", columnname, value = FALSE)
#第一遍跑过    DF.classifications_columname
#第一遍跑过    sc_persampledata@meta.data <- sc_persampledata@meta.data[,- DF.classifications_columname]
#第一遍跑过    colnames(sc_persampledata@meta.data)
#第一遍跑过
#第一遍跑过    ####  remove SCT_snn_res_columname ####
#第一遍跑过    columnname <- colnames(sc_persampledata@meta.data)
#第一遍跑过    SCT_snn_res_columname <- grep("^SCT_snn_res", columnname, value = FALSE)
#第一遍跑过    SCT_snn_res_columname
#第一遍跑过    sc_persampledata@meta.data <- sc_persampledata@meta.data[,- SCT_snn_res_columname]
#第一遍跑过    colnames(sc_persampledata@meta.data)
#第一遍跑过
#第一遍跑过    columnname <- colnames(sc_persampledata@meta.data)
#第一遍跑过    print(paste0(samplename, "columnname: ", columnname))
#第一遍跑过    columnname 
#第一遍跑过    return(sc_persampledata)
#第一遍跑过    })
#第一遍跑过
#第一遍跑过
#第一遍跑过
#第一遍跑过    names(CN1_23filter_PSCT_pure_persample_DF)
#第一遍跑过    names(CN1_23filter_PSCT_pure_persample_DF)<- namelist
#第一遍跑过    names(CN1_23filter_PSCT_pure_persample_DF)
#第一遍跑过    saveRDS(CN1_23filter_PSCT_pure_persample_DF, paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_PSCT_pure_persample_DF_230426.rds"))
#第一遍跑过
#第一遍跑过
#第一遍跑过#### merge sceList1 data ####
#第一遍跑过    # save(sceList1_SCT_preann_pure_persample_DF, file = paste0(filesDir, "Results/02Clustering/SaveData/sceList1_SCT_preann_pure_persample_DF.RData"))
#第一遍跑过    CN1_23filter_PSCT_pure_merge_DF <- merge(CN1_23filter_PSCT_pure_persample_DF[[1]], y = CN1_23filter_PSCT_pure_persample_DF[2:length(CN1_23filter_PSCT_pure_persample_DF)], project = "CN1_23filter_PSCT_merge")
#第一遍跑过
#第一遍跑过    if(!dir.exists(paste0(filesDir, "Results/03Integration/SaveData/merge/"))){dir.create(paste0(filesDir, "Results/03Integration/SaveData/merge/"), recursive = TRUE)}
#第一遍跑过
#第一遍跑过    saveRDS(CN1_23filter_PSCT_pure_merge_DF, file = paste0(filesDir, "Results/03Integration/SaveData/merge/CN1_23filter_PSCT_pure_merge_DF_230426.rds"))
#第一遍跑过 
#第一遍跑过#### remove doublet based on column Doublet  ####
#第一遍跑过    CN1_23filter_PSCT_pure_merge_rmDoublet <- subset(CN1_23filter_PSCT_pure_merge_DF, subset = Doublet %in% c("Doublet_Zero_Confidience", "Doublet_Low_Confidience") )  ##保留的细胞放大#21947*10937
#第一遍跑过    saveRDS(CN1_23filter_PSCT_pure_merge_rmDoublet, paste0(filesDir, "Results/03Integration/SaveData/merge/CN1_23filter_PSCT_pure_merge_rmDoublet_230426.rds"))
#第一遍跑过




    #### 
    CN1_23filter_PSCT_pure_merge_rmDoublet <- readRDS(paste0(filesDir, "Results/03Integration/SaveData/merge/CN1_23filter_PSCT_pure_merge_rmDoublet_230426.rds"))



#### UM ####
#sceList2UM_PSCT_pure_merge_rmDoublet <- readRDS(paste0(filesDir, "/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis2_with_UM/Results/3integration/SaveData/merge/sceList2UM_PSCT_pure_merge_rmDoublet.rds"))

# Error in gzfile(file, "rb") : cannot open the connection
# Calls: readRDS -> gzfile
# In addition: Warning message:
# In gzfile(file, "rb") :
#   cannot open compressed file '/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3//share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis2_with_UM/Results/3integration/SaveData/merge/sceList2UM_PSCT_pure_merge_rmDoublet.rds', probable reason 'No such file or directory'
# Execution halted

sceList2UM_PSCT_pure_merge_rmDoublet <- readRDS(paste0("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis2_with_UM/Results/3integration/SaveData/merge/sceList2UM_PSCT_pure_merge_rmDoublet.rds"))

table(sceList2UM_PSCT_pure_merge_rmDoublet@meta.data$Doublet)
colnames(sceList2UM_PSCT_pure_merge_rmDoublet@meta.data)

#### CN2 ####
sceList3CN2_PSCT_pure_merge_rmDoublet <- readRDS(paste0("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis3_with_CN2/Results/3integration/SaveData/merge/sceList3CN2_PSCT_pure_merge_rmDoublet.rds"))

table(sceList3CN2_PSCT_pure_merge_rmDoublet@meta.data$Doublet)
colnames(sceList3CN2_PSCT_pure_merge_rmDoublet@meta.data)

CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet <- merge(CN1_23filter_PSCT_pure_merge_rmDoublet, y = c(sceList2UM_PSCT_pure_merge_rmDoublet, sceList3CN2_PSCT_pure_merge_rmDoublet), project = "CN1_23samples_PUMCHUMCN2_PSCT_merge")
print(colnames(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet@meta.data))
table(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet$Doublet)
table(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet$orig.ident)
saveRDS(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet, paste0(filesDir, "Results/03Integration/SaveData/merge/CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet_230426.rds"))
    
table(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet@meta.data$Doublet)

rm(list = ls())


# Doublet_Low_Confidience
# Doublet_Zero_Confidience


    #rm(list = ls())



