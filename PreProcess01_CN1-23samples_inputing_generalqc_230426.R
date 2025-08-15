### 01PreProcess ####
##PreProcess01_CN1-23samples_inputing_generalqc_230426.R
#!/bin/bash
#SBATCH -J  PreProcess01_CN1-23samples_inputing_generalqc_230426
#SBATCH -D /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3 
#SBATCH -p compute
#SBATCH --output=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/01PreProcess/PreProcess01_CN1-23samples_inputing_generalqc_230426/PreProcess01_CN1-23samples_inputing_generalqc_230426.out         
#SBATCH --error=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/01PreProcess/PreProcess01_CN1-23samples_inputing_generalqc_230426/PreProcess01_CN1-23samples_inputing_generalqc_230426.err                  
#SBATCH --nodes=1                          
#SBATCH --ntasks-per-node=56 
#SBATCH --mem=500000    
#SBATCH --nodelist=n01


#**** 这个在MyeloidAnalysis05_CN1UMCN2_annoting_filting3_resfinding 文件后面

# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results


# ###mkdir /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/01PreProcess/PreProcess01_CN1-23samples_inputing_generalqc_230426/

# module load anaconda_global/anaconda_22
# source activate /share/home/sigrid/biosoft_sigrid/anaconda3/envs/R4
# R CMD BATCH --no-restore /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles/01PreProcess/PreProcess01_CN1-23samples_inputing_generalqc_230426.R  /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/01PreProcess/PreProcess01_CN1-23samples_inputing_generalqc_230426/PreProcess01_CN1-23samples_inputing_generalqc_230426.Rout
# conda deactivate
# module unload anaconda_global/anaconda_22 

# Error in paste0(filesDir, "Results/01PreProcess/Figures/CN1_23UM_MReRNA/",  : 
#   argument is missing, with no default
# Calls: res_persample -> dir.exists -> paste0
# Execution halted

#### prepare core&environment ####
rm(list=ls())
sink("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/01PreProcess/PreProcess01_CN1-23samples_inputing_generalqc_230426/PreProcess01_CN1-23samples_inputing_generalqc_230426.txt", append=TRUE, split=TRUE )
pdf("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/01PreProcess/PreProcess01_CN1-23samples_inputing_generalqc_230426/PreProcess01_CN1-23samples_inputing_generalqc_230426.pdf") #
library(future)
plan("multicore", workers = 55) 
options(future.globals.maxSize= 9000*1024^2)
future.seed = NULL
nbrOfWorkers()



# save(CN1_23UM_SCT1_1, file = paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23UM_SCT1_1.RData"))
# save(CN1_23UM_SCT1_3, file = paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23UM_SCT1_3.RData"))
# save(CN1_23UM_SCT1_3, file = paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23UM_SCT1_3.RData"))

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
#    suppressPackageStartupMessages(library(sctransform))
#    library(viridis)
# library(maps)
# library(fields)
# library(spam)
#    suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(harmony))
#suppressPackageStartupMessages(library(clusterProfiler))
#suppressPackageStartupMessages(library(org.Hs.eg.db))
options(stringsAsFactors = F)
filesDir<-'/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/'
setwd(filesDir)
getwd()




#### prepare folderdirectory ####
	if(!dir.exists(paste0(filesDir, "Results/01PreProcess/SaveData/"))){dir.create(paste0(filesDir, "Results/01PreProcess/SaveData/"), recursive = TRUE)}
	if(!dir.exists(paste0(filesDir, "Results/01PreProcess/Figures/"))){dir.create(paste0(filesDir, "Results/01PreProcess/Figures/"), recursive = TRUE)}
	if(!dir.exists(paste0(filesDir, "Results/01PreProcess/Tables/"))){dir.create(paste0(filesDir, "Results/01PreProcess/Tables/"), recursive = TRUE)}

#### establish rawdata_list ####
	dir="/share/home/sigrid/IPMN_PROGRAM/Rawdata/matrixdata_nointrons_c/matrixdata_nointrons"
	samples=list.files(dir)
	samples
	CN1_23 = lapply(samples,function(pro){
		  print(pro)
		  folder=file.path(dir,pro)	
		  sce.data <- Read10X(folder)
		  colnames(sce.data) <- paste0(pro, "_",colnames(sce.data))
		  head(colnames(sce.data))
		  sce=CreateSeuratObject(	counts = sce.data,
		                     		project = pro,
		                     		min.cells = 3,
		                    		min.features = 200)
		  print(head(sce))
		  return(sce)
	})


	##nomenclature
	
	names(CN1_23)  
	samples
	names(CN1_23) =  samples
	saveRDS(CN1_23, paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23sampleraw_230426.rds"))
	head(CN1_23[[1]]@meta.data)
	head(CN1_23[[2]]@meta.data)




filesDir<-'/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/'
setwd(filesDir)
getwd()

    if(!dir.exists(paste0(filesDir, "Results/01PreProcess/SaveData/"))){dir.create(paste0(filesDir, "Results/01PreProcess/SaveData/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/01PreProcess/Figures/"))){dir.create(paste0(filesDir, "Results/01PreProcess/Figures/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/01PreProcess/Tables/"))){dir.create(paste0(filesDir, "Results/01PreProcess/Tables/"), recursive = TRUE)}



#### filter rawdata 函数 ####
  preprocess_filter_sc_func <- function( sce_rawdata, samplename, featureMin = 200, featureMax = 8000, 
                                         CountRNAMin = 1000, percentMT = 15) {  
    #### filter rawdata 函数 ####
    #### library调用 ####
    #### filter rawdata ####
    sce_rawdata[["percent.mt"]] <- PercentageFeatureSet(sce_rawdata, pattern = "^MT-")
    fivenum(sce_rawdata$percent.mt)
    sce_rawdata_filter <- subset(sce_rawdata, subset = nFeature_RNA >= featureMin & nFeature_RNA <= featureMax & nCount_RNA >= CountRNAMin & percent.mt <= percentMT ) #其他文献中都没有添加≥
    # selected_f <- rownames(sce_rawdata_filter)[Matrix::rowSums(sce_rawdata_filter@assays$RNA@counts > 0 ) > 3] # 至少在3个细胞中表达 
    # sce_rawdata_filter<-subset(sce_rawdata_filter,features = selected_f)
    
    #### plot the pre-post-filtered plots ####
    plot0_1 <- VlnPlot(sce_rawdata, features = "nFeature_RNA", group.by = "orig.ident", raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nFeature_RNA", subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot0_2 <- VlnPlot(sce_rawdata, features =  "nCount_RNA",  group.by = "orig.ident", raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nCount_RNA", subtitle =paste0(samplename,"_preQC"))  + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot0_3 <- VlnPlot(sce_rawdata, features = "percent.mt", group.by = "orig.ident", raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "percent.mt", subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot0 = plot0_1|plot0_2|plot0_3 

    plot1 <- FeatureScatter(sce_rawdata, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
                            labs(subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot2 <- FeatureScatter(sce_rawdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
                            labs(subtitle =paste0(samplename,"_preQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    
    plot3_1 <- VlnPlot(sce_rawdata_filter, features = "nFeature_RNA", group.by = "orig.ident", raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nFeature_RNA", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot3_2 <- VlnPlot(sce_rawdata_filter, features =  "nCount_RNA",  group.by = "orig.ident", raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "nCount_RNA", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot3_3 <- VlnPlot(sce_rawdata_filter, features = "percent.mt", group.by = "orig.ident", raster=FALSE) + geom_boxplot(width=.2,col="black",fill="white") + 
                            labs(title = "percent.mt", subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot3 = plot3_1|plot3_2|plot3_3 

    plot4 <- FeatureScatter(sce_rawdata_filter, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
                                labs( subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    plot5 <- FeatureScatter(sce_rawdata_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
                                labs( subtitle = paste0(samplename,"_posQC")) + NoLegend() + theme(plot.subtitle = element_text(hjust = 0.5))
    
    if(!is.null(dev.list())){dev.off()}
    pdf(file = paste0(filesDir, "Results/01PreProcess/Figures/", samplename, "_1_pre_post_QC.pdf"), width = 10, height = 6)
    print(plot0) #待测试)
    print(plot1 + plot2)
    print(plot3)
    print(plot4 + plot5)
    dev.off()
    
    
    #### statistic the QC results  ####
    qc <- c(ncol(sce_rawdata), nrow(sce_rawdata), median(sce_rawdata$nFeature_RNA), median(sce_rawdata$nCount_RNA), median(sce_rawdata$percent.mt),
            mean(sce_rawdata$nFeature_RNA), mean(sce_rawdata$nCount_RNA), mean(sce_rawdata$percent.mt))
    qcPost <- c(ncol(sce_rawdata_filter), nrow(sce_rawdata_filter), median(sce_rawdata_filter$nFeature_RNA), median(sce_rawdata_filter$nCount_RNA), 
                median(sce_rawdata_filter$percent.mt),mean(sce_rawdata_filter$nFeature_RNA), mean(sce_rawdata_filter$nCount_RNA), mean(sce_rawdata_filter$percent.mt))
    stat <- data.frame(qcBefore = qc, qcPost = qcPost)
    rownames(stat) <- c("cells", "genes", "medianFeatures", "medianCounts", "medianmt", "meanFeatures", "meanCounts", "meanmt")
    stat$percent <- 100*qcPost/qc
    write.csv(stat, file = paste0(filesDir, "Results/01PreProcess/Tables/", samplename, "_QCstat.csv"), quote = FALSE)
    
    return(sce_rawdata_filter)
  }

#### filter rawdata QC ####

    CN1_23 <-  readRDS(paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23sampleraw_230426.rds"))
    namelist <- names(CN1_23)
    namelist
    CN1_23_filter = lapply(namelist,function(pro){
              print(pro)
              scrawdata <- CN1_23[[pro]]
              sce_rawdata_filter <- preprocess_filter_sc_func(scrawdata, pro, featureMin = 200, featureMax = 8000, CountRNAMin = 1000, percentMT = 15)
              return(sce_rawdata_filter)
        })
    names(CN1_23_filter)  
    namelist
    names(CN1_23_filter) =  namelist
    for (i in names(CN1_23_filter)) {
        sce <- CN1_23_filter[[i]]
        print(paste0(i, " cell_filter count is ", ncol(sce)))
    }
    head(CN1_23_filter[[1]]@meta.data)
    colnames(CN1_23_filter[[1]]@meta.data)


#### 调整顺序 ####
 
    addingCellType1<-readxl::read_xlsx(paste0(filesDir, "Results/10TotalAnalysis/Tables/CN11UMCN2/sampleinformation/Annotation/CN11UMCN2CN3_AddingGroupInfo_persample_230426.xlsx"))  
    #sceList1_addingCellType1_persample.xlsx
    addingCellType1
    addingCellType1_df<-as.data.frame(addingCellType1)
    rm(addingCellType1)
    head(addingCellType1_df)

  ## rename variables ##
    namelist_filter <- addingCellType1_df$orig.ident[1:23]
    namelist_filter
    CN1_23_filter_order = lapply(namelist_filter, function(pro){
                                        print(pro)
                                        sce_rawdata <- CN1_23_filter[[pro]]
                                  return(sce_rawdata)
        })
    names(CN1_23_filter_order)  
    CN1_23_filter_order
    names(CN1_23_filter_order) =  namelist_filter

    names(CN1_23_filter_order)  
    colnames(CN1_23_filter_order[[1]]@meta.data)  

    saveRDS(CN1_23_filter_order, paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23_filter_order_230426.rds"))



rm(list = ls())