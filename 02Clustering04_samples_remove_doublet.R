#### prepare core&environment ####
rm(list=ls())
library(future)
plan("multicore", workers = 55) 
options(future.globals.maxSize= 9000000*1024^2)
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

    CN1_23filter_PSCT_pure_merge_rmDoublet <- readRDS(paste0(filesDir, "Results/03Integration/SaveData/merge/CN1_23filter_PSCT_pure_merge_rmDoublet_230426.rds"))



sceList2UM_PSCT_pure_merge_rmDoublet <- readRDS(paste0("sceList2UM_PSCT_pure_merge_rmDoublet.rds"))

table(sceList2UM_PSCT_pure_merge_rmDoublet@meta.data$Doublet)
colnames(sceList2UM_PSCT_pure_merge_rmDoublet@meta.data)

#### CN2 ####
sceList3CN2_PSCT_pure_merge_rmDoublet <- readRDS(paste0("sceList3CN2_PSCT_pure_merge_rmDoublet.rds"))

table(sceList3CN2_PSCT_pure_merge_rmDoublet@meta.data$Doublet)
colnames(sceList3CN2_PSCT_pure_merge_rmDoublet@meta.data)

CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet <- merge(CN1_23filter_PSCT_pure_merge_rmDoublet, y = c(sceList2UM_PSCT_pure_merge_rmDoublet, sceList3CN2_PSCT_pure_merge_rmDoublet), project = "CN1_23samples_PUMCHUMCN2_PSCT_merge")
print(colnames(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet@meta.data))
table(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet$Doublet)
table(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet$orig.ident)
saveRDS(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet, paste0(filesDir, "CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet_230426.rds"))
    
table(CN1_23samples_PUMCHUMCN2_PSCT_pure_merge_rmDoublet@meta.data$Doublet)

rm(list = ls())
