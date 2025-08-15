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
suppressPackageStartupMessages(library(harmony))
options(stringsAsFactors = F)
filesDir<-'/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/'
setwd(filesDir)
getwd()

    if(!dir.exists(paste0(filesDir, "Results/02Clustering/SaveData/"))){dir.create(paste0(filesDir, "Results/02Clustering/SaveData/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Figures/"))){dir.create(paste0(filesDir, "Results/02Clustering/Figures/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Tables/"))){dir.create(paste0(filesDir, "Results/02Clustering/Tables/"), recursive = TRUE)}

 
        CN1_23_filter <- readRDS(paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23_filter_order_230426.rds"))
    # load(paste0(filesDir, "Results/02Clustering/SaveData/CN1_23_filter.RData"))

#### SCTranform for each sample ####
    namelist <- names(CN1_23_filter)
    print(namelist)
    CN1_23filter_SCT = lapply(namelist,function(samplename){
              print(samplename)
              sc_filter_sampledata <- CN1_23_filter[[samplename]]
              sc_filter_sampledata <- SCTransform(sc_filter_sampledata, vars.to.regress = c( "percent.mt"), verbose = FALSE)
              return(sc_filter_sampledata)
        })
    rm(CN1_23_filter)
    names(CN1_23filter_SCT)  
    namelist
    names(CN1_23filter_SCT) =  namelist
    names(CN1_23filter_SCT) 

    saveRDS(CN1_23filter_SCT, paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_SCT_230426.rds"))
    CN1_23filter_SCT <- readRDS(paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_SCT_230426.rds"))

#### explore optimal npcs for RUNPCA  ####
  ####  functions explore optimal npcs for RUNPCA  ####   
    pcafind_persample <- function(sc_persampledata, samplename, method = c("RNA", "SCT"), npc_num = 50 ){
        if(method == "SCT"){
            DefaultAssay(sc_persampledata) <- "SCT"
        }else{
            DefaultAssay(sc_persampledata) <- "RNA"
        }
        npc_str <- as.character(npc_num)
        sc_persampledata <- RunPCA(sc_persampledata, npcs = npc_num, verbose = FALSE)
        ### JackStrawPlot Slow slow slow
            #sc_persampledata  <- JackStraw(object = sc_persampledata, dims = 50)
            #sc_persampledata  <- ScoreJackStraw(sc_persampledata, dims = 1:50)
            #JackStrawPlot(sc_persampledata, dims =1:50)
        #### ElbowPlot for each sample ####
            #The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
            # Determine percent of variation associated with each PC    
        pct <- sc_persampledata[["pca"]]@stdev / sum(sc_persampledata[["pca"]]@stdev) * 100
            # Calculate cumulative percents for each PC
        cumu <- cumsum(pct)
            # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
        col1 <- which(cumu > 90 & pct < 5)[1]
        col1
        # Determine the difference between variation of PC and subsequent PC
        col2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
        # last point where change of % of variation is more than 0.1%.
        col2
        pcs <- min(col1, col2)
        p_elbowplot <- ElbowPlot(object = sc_persampledata, ndims = npc_num) + 
                            labs(title = paste0(samplename,"_elbowplot"), subtitle = paste0("optimal npc:",as.character(col1)," AND",as.character(col2))) ##
        print(p_elbowplot)


            # # Slow slow slow
            # sc_persampledata  <- JackStraw(object = sc_persampledata, dims = npc_num)
            # sc_persampledata  <- ScoreJackStraw(sc_persampledata, dims = 1:npc_num)
            # p_jackJackStrawPlot <- JackStrawPlot(object = sc_persampledata, dims = 1:npc_num)
            # print(p_jackJackStrawPlot)
            # rm(p_jackJackStrawPlot)
            # Error in JackStraw(object = sc_persampledata, dims = npc_num) : 
            #       JackStraw cannot be run on SCTransform-normalized data.
            #              Please supply a non-SCT assay.
            #     Calls: pcafind_persample -> JackStraw

        return(sc_persampledata)
    }
    
    #load(paste0(filesDir, "Results/1preprocess/SaveData/sceList1_SCT.RData"))
    sceList1_SCT <- CN1_23filter_SCT


    head(sceList1_SCT[[1]])
    npc_list <- c(50,100)
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Figures/PSCT_npcfind/"))){dir.create(paste0(filesDir, "Results/02Clustering/Figures/PSCT_npcfind/"), recursive = TRUE)}

for (npc_num in npc_list) {
    if(!is.null(dev.list())){dev.off()}
    pdf(file = paste0(filesDir, "Results/02Clustering/Figures/PSCT_npcfind/", "CN1_23filter_SCT_PCA-",as.character(npc_num),"_test.pdf"), width = 20, height = 12)

    for (i in 1:length(sceList1_SCT)) {
        sc_persampledata <- sceList1_SCT[[i]]
        samplename <- names(sceList1_SCT)[i]
        print(samplename)
        pcafind_persample(sc_persampledata = sc_persampledata,samplename = samplename, method = "SCT", npc_num = npc_num)
        print(paste0("finish ", samplename, " ",npc_num))
    }
   dev.off()  
}

    rm(list=ls())



