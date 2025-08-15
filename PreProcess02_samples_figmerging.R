#### prepare core&environment ####
rm(list=ls())
library(future)
plan("multicore", workers = 55) 
options(future.globals.maxSize= 9000*1024^2)
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
suppressPackageStartupMessages(library(harmony))

options(stringsAsFactors = F)
filesDir<-'/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/'
setwd(filesDir)
getwd()

#### prepare folderdirectory ####


    if(!dir.exists(paste0(filesDir, "Results/01PreProcess/SaveData/"))){dir.create(paste0(filesDir, "Results/01PreProcess/SaveData/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/01PreProcess/Figures/"))){dir.create(paste0(filesDir, "Results/01PreProcess/Figures/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/01PreProcess/Tables/"))){dir.create(paste0(filesDir, "Results/01PreProcess/Tables/"), recursive = TRUE)}




    CN1_23_raw <-  readRDS(paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23sampleraw_230426.rds"))
    CN1_23_filter <- readRDS(paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23samplefilter_230426.rds")) 

    head(CN1_23_raw[[1]])
 
    addingCellType1<-readxl::read_xlsx(paste0(filesDir, "Results/10TotalAnalysis/Tables/CN11UMCN2/sampleinformation/Annotation/CN11UMCN2CN3_AddingGroupInfo_persample_230426.xlsx"))  
    #sceList1_addingCellType1_persample.xlsx
    addingCellType1
    addingCellType1_df<-as.data.frame(addingCellType1)
    rm(addingCellType1)
    head(addingCellType1_df)

    namelist1_raw <- addingCellType1_df$orig.ident[1:23]
 #add percent.mt
    CN1_23_raw_order = lapply(namelist1_raw, function(pro){
                        print(pro)
                        sce_rawdata <- CN1_23_raw[[pro]]
                        sce_rawdata[["percent.mt"]] <- PercentageFeatureSet(sce_rawdata, pattern = "^MT-") 
                  return(sce_rawdata)
        })
    names(CN1_23_raw_order)  
    CN1_23_raw_order
    names(CN1_23_raw_order) =  namelist1_raw
    
    head(CN1_23_raw_order[[1]])
    print(names(CN1_23_raw_order))
    saveRDS(CN1_23_raw_order, paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23_raw_order_230426.rds"))

    CN1_23_raw <- CN1_23_raw_order
  ## rename variables ##
    namelist_filter <- addingCellType1_df$orig.ident[1:23]
    CN1_23_filter <- readRDS(paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23_filter_order_230426.rds"))
    CN1_23_filter_order <- CN1_23_filter

     for (samplename in namelist_filter) {
        
        scedata <- CN1_23_filter_order[[samplename]]

        saveRDS(scedata, paste0(filesDir, "Results/01PreProcess/SaveData/",samplename,"_filter_230426.rds"))
    }    




    preqc_ccounts <- numeric()
    posqc_ccounts <- numeric()
    cellcounts_rownames <- character()
    
    for (i in namelist1_raw) {
      preqc_ccounts <- append(preqc_ccounts, ncol(CN1_23_raw[[i]]))
      posqc_ccounts <- append(posqc_ccounts, ncol(CN1_23_filter[[i]]))
      cellcounts_rownames <- append(cellcounts_rownames, i)
    }
    
    preqc_ccounts
    posqc_ccounts
    cellcounts_rownames

    sc_preprocess_qc1_cellcount_df <- data.frame(preQC_cellcounts = preqc_ccounts , posQC_cellcounts = posqc_ccounts)
    rownames(sc_preprocess_qc1_cellcount_df) <- cellcounts_rownames
    sc_preprocess_qc1_cellcount_df$percent <- 100*posqc_ccounts/preqc_ccounts
    write.csv(sc_preprocess_qc1_cellcount_df, file = paste0(filesDir, "Results/01PreProcess/Tables/", "sc_preprocess_qc1_cellcount.csv"), quote = FALSE)

#### create sc_merge_data ####
    # CN1_23_raw <- merge( CN1_23_raw[[1]],
    #                              y=c(CN1_23_raw[[2]], CN1_23_raw[[3]], 
    #                                  CN1_23_raw[[4]], CN1_23_raw[[5]], CN1_23_raw[[6]], CN1_23_raw[[7]], 
    #                                  CN1_23_raw[[8]], CN1_23_raw[[9]], CN1_23_raw[[10]],
    #                                  CN1_23_raw[[11]], CN1_23_raw[[12]], CN1_23_raw[[13]], CN1_23_raw[[14]])
    #                                 #,add.cell.ids=namelist1_raw
    #                         )

    CN1_23_raw_merge <- merge(CN1_23_raw[[1]], y = CN1_23_raw[2:length(CN1_23_raw)], project = "CN1_23_raw_merge")
    CN1_23_filter_merge <- merge(CN1_23_filter[[1]], y = CN1_23_filter[2:length(CN1_23_filter)], project = "CN1_23_filter_merge")
                                    #,add.cell.ids=namelist1_filter

    save(CN1_23_filter_merge, file = paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23_filter_merge_230426.RData"))

    head(CN1_23_raw_merge)
    tail(CN1_23_raw_merge)
    head(CN1_23_filter_merge)
    tail(CN1_23_filter_merge)




#### plot the pre-post-filtered plots ####
    samplename<- "sc_analysis6_CN1_23samples_merge"
    plot0_1 <- VlnPlot(CN1_23_raw_merge, features = "nFeature_RNA", group.by = "orig.ident", split.by = "orig.ident", pt.size = 0.1) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nFeature_RNA", subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)
    plot0_2 <- VlnPlot(CN1_23_raw_merge, features =  "nCount_RNA",  group.by = "orig.ident", split.by = "orig.ident", pt.size = 0.1) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nCount_RNA", subtitle =paste0(samplename,"_preQC"))  + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)
    plot0_3 <- VlnPlot(CN1_23_raw_merge, features = "percent.mt", group.by = "orig.ident", split.by = "orig.ident", pt.size = 0.1) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "percent.mt", subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)

    plot1 <- FeatureScatter(CN1_23_raw_merge, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1) + 
                            labs(subtitle =paste0(samplename,"_preQC")) + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
    plot2 <- FeatureScatter(CN1_23_raw_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1) + 
                            labs(subtitle =paste0(samplename,"_preQC")) + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
   

    plot3_1 <- VlnPlot(CN1_23_filter_merge, features = "nFeature_RNA", group.by = "orig.ident", split.by = "orig.ident", pt.size = 0) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nFeature_RNA", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5))  + xlab(NULL)
    plot3_2 <- VlnPlot(CN1_23_filter_merge, features =  "nCount_RNA",  group.by = "orig.ident", split.by = "orig.ident", pt.size = 0) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nCount_RNA", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)
    plot3_3 <- VlnPlot(CN1_23_filter_merge, features = "percent.mt", group.by = "orig.ident", split.by = "orig.ident", pt.size = 0) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "percent.mt", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)


    plot4 <- FeatureScatter(CN1_23_filter_merge, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1) + 
    							labs( subtitle = paste0(samplename,"_posQC"))  + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
    plot5 <- FeatureScatter(CN1_23_filter_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1) +
    							labs( subtitle = paste0(samplename,"_posQC"))  + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
    a<- plot_grid(plot1, plot2, plot4, plot5, ncol = 2, align = "vh")
    
    if(!is.null(dev.list())){dev.off()}
    pdf(file = paste0(filesDir, "Results/01PreProcess/Figures/", samplename, "_pre_post_QC_merge_1-1.pdf"), width = 10, height = 6)
    print(plot0_1 / plot3_1) 
    print(plot0_2 / plot3_2)
    print(plot0_3 / plot3_3)
    print(a)
    print(plot4)
    print(plot5)
    dev.off()
    
    if(!is.null(dev.list())){dev.off()}
    pdf(file = paste0(filesDir, "Results/01PreProcess/Figures/", samplename, "_pre_post_QC_merge_1-2.pdf"), width = 20, height = 12)
    print(plot0_1 / plot3_1) 
    print(plot0_2 / plot3_2)
    print(plot0_3 / plot3_3)
    print(a)
    print(plot4)
    print(plot5)
    dev.off()

    print("finish pre-post-filtered 1 plots ")
    
    samplename<- "CN1_23UM_cPDAC_merge"
    plot0_1 <- VlnPlot(CN1_23_raw_merge, features = "nFeature_RNA", group.by = "orig.ident", split.by = "orig.ident", raster=FALSE ) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nFeature_RNA", subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)
    plot0_2 <- VlnPlot(CN1_23_raw_merge, features =  "nCount_RNA",  group.by = "orig.ident", split.by = "orig.ident", raster=FALSE ) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nCount_RNA", subtitle =paste0(samplename,"_preQC"))  + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)
    plot0_3 <- VlnPlot(CN1_23_raw_merge, features = "percent.mt", group.by = "orig.ident", split.by = "orig.ident", raster=FALSE ) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "percent.mt", subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)

    plot1 <- FeatureScatter(CN1_23_raw_merge, feature1 = "nCount_RNA", feature2 = "percent.mt",  raster=FALSE) + 
                            labs(subtitle =paste0(samplename,"_preQC")) + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
    plot2 <- FeatureScatter(CN1_23_raw_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",  raster=FALSE) + 
                            labs(subtitle =paste0(samplename,"_preQC")) + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
   

    plot3_1 <- VlnPlot(CN1_23_filter_merge, features = "nFeature_RNA", group.by = "orig.ident", split.by = "orig.ident",  raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nFeature_RNA", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5))  + xlab(NULL)
    plot3_2 <- VlnPlot(CN1_23_filter_merge, features =  "nCount_RNA",  group.by = "orig.ident", split.by = "orig.ident",  raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nCount_RNA", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)
    plot3_3 <- VlnPlot(CN1_23_filter_merge, features = "percent.mt", group.by = "orig.ident", split.by = "orig.ident",  raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "percent.mt", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 0,hjust = 0.5)) + xlab(NULL)


    plot4 <- FeatureScatter(CN1_23_filter_merge, feature1 = "nCount_RNA", feature2 = "percent.mt",  raster=FALSE) + 
                                labs( subtitle = paste0(samplename,"_posQC"))  + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
    plot5 <- FeatureScatter(CN1_23_filter_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",  raster=FALSE) +
                                labs( subtitle = paste0(samplename,"_posQC"))  + theme(plot.subtitle = element_text(hjust = 0.5)) + theme(legend.title=element_blank(), legend.key.size = unit(10, "pt"), legend.text=element_text(size=5))
    a<- plot_grid(plot1, plot2, plot4, plot5, ncol = 2, align = "vh")
    if(!is.null(dev.list())){dev.off()}
    pdf(file = paste0(filesDir, "Results/01PreProcess/Figures/", samplename, "_pre_post_QC_merge_2-1.pdf"), width = 10, height = 6)
    print(plot0_1 / plot3_1) 
    print(plot0_2 / plot3_2)
    print(plot0_3 / plot3_3)
    print(a)
    print(plot4)
    print(plot5)
       dev.off()

    if(!is.null(dev.list())){dev.off()}
    pdf(file = paste0(filesDir, "Results/01PreProcess/Figures/", samplename, "_pre_post_QC_merge_2-2.pdf"), width = 20, height = 12)
    print(plot0_1 / plot3_1) 
    print(plot0_2 / plot3_2)
    print(plot0_3 / plot3_3)
    print(a)
    print(plot4)
    print(plot5)
    dev.off()

print("finish total ")

rm(list = ls())






