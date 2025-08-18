### 07MyeloidAnalysis ####
##MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin.R
#!/bin/bash
#SBATCH -J  MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin
#SBATCH -D /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3 
#SBATCH -p compute
#SBATCH --output=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/07MyeloidAnalysis/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin.out         
#SBATCH --error=/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/07MyeloidAnalysis/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin.err                  
#SBATCH --nodes=1                          
#SBATCH --ntasks-per-node=56 
#SBATCH --mem=500000    
#SBATCH --nodelist=n02



# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles
# /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results


# ###mkdir /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/07MyeloidAnalysis/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin/

# module load anaconda_global/anaconda_22
# source activate /share/home/sigrid/biosoft_sigrid/anaconda3/envs/R1
# R CMD BATCH --no-restore /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/sc_scriptfiles/07MyeloidAnalysis/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin.R  /share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/07MyeloidAnalysis/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin.Rout
# conda deactivate
# module unload anaconda_global/anaconda_22 

# Error in paste0(filesDir, "Results/07MyeloidAnalysis/Figures/sceList2UM_MReRNA/",  : 
#   argument is missing, with no default
# Calls: res_persample -> dir.exists -> paste0
# Execution halted

#### prepare core&environment ####
    rm(list=ls())
    # sink("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/07MyeloidAnalysis/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin.txt", append=TRUE, split=TRUE )
    # pdf("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/07MyeloidAnalysis/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin/MyeloidAnalysis15-01_CN11UMCN2_Macrophage_exp_score-samplestage2_240527_inwin.pdf") #
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
      library(ggalluvial)
    library(RColorBrewer)
    library(ggpubr)
    library(rstatix)
    library(ggsignif)
    library(ggalluvial)
 
 #    suppressPackageStartupMessages(library(sctransform))
 #    library(viridis)
    # library(maps)
    # library(fields)
    # library(spam)
 #    suppressPackageStartupMessages(library(DoubletFinder))
    #suppressPackageStartupMessages(library(harmony))
    #suppressPackageStartupMessages(library(clusterProfiler))
    #suppressPackageStartupMessages(library(org.Hs.eg.db))
    options(stringsAsFactors = F)
    filesDir<-'I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/07MyeloidAnalysis/'
    setwd(filesDir)
    getwd()
    identlevel <- "Macrophage"
    samplename <- "IPMClineage20sample_240527"
    method = "RNA"
    
    FigDir <- paste0(filesDir, "Figures/", samplename,"/",identlevel,"/")
    TableDir <- paste0(filesDir, "Tables/", samplename,"/",identlevel,"/")
    DataDir <- paste0(filesDir, "SaveData/", samplename,"/",identlevel,"/")
    
    
# Th1 = c("IL2", )
     if(!dir.exists(paste0(FigDir))){dir.create(paste0(FigDir), recursive = TRUE)}
     if(!dir.exists(paste0(TableDir))){dir.create(paste0(TableDir), recursive = TRUE)}
     if(!dir.exists(paste0(DataDir))){dir.create(paste0(DataDir), recursive = TRUE)}

#### INPUT DATA ####


totalcelldata <- readRDS( paste0("I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/07MyeloidAnalysis/SaveData/Macrophage/IPMClineage_231116/IPMClineage_macrophage_seurat.rds"))
    ncol(totalcelldata)
    table(totalcelldata$SampleStage) #5777
    totalcelldata <- subset(totalcelldata, subset = SampleStage %in% c("CTRL","LGD-IPMN","HGD-IPMN","IPMC"))
    ncol(totalcelldata)
  totalcelldata$SampleStage2 <- as.character(totalcelldata$SampleStage2)
  totalcelldata$SampleStage2[totalcelldata$SampleStage2== "CTRL"] <- "UNIN"
  totalcelldata$SampleStage2[totalcelldata$SampleStage2== "IPMC"] <- "PDAC"
  table(totalcelldata$SampleStage2)
  table(totalcelldata$SampleStage)

    totalcelldata@meta.data$SampleStage2 <- factor(totalcelldata@meta.data$SampleStage2, levels = c("UNIN","IPMN","PDAC"))
    table(totalcelldata$SubCellType0_CN11UMCN2)
    table(totalcelldata$SampleStage2)
    table(totalcelldata$CellType_M1CN11UMCN2)
    table(totalcelldata$SubCellType0_CN11UMCN2)
    table(totalcelldata$SubCellType1_CN11UMCN2)
    table(totalcelldata$SubCellTypeup_CN11UMCN2)

    scorelist = list(
    M1 = c("NOS2","IL12A" ,"IL12B" ,"IL6"  ,"TNF" ,"CD86","IL1B" ,"CXCL9","CXCL10" ,"CXCL11","CXCL12","STAT1","AIM2","IL23"),
    M2 = c("MRC1", "ARG1" , "CCL22" , "CCL17" , "CCL24","CCL18","CXCR4","NRP1","VEGFA","VEGFB",
           "PDGFB","LLGL1","MMP12","CD209","IL10","WNT5A","ADORA3","STAT6","SOCS3","IRF4",
           "LTA","TFRC","IL4"),
    
    TAM = c("IL1B","IL6","TNF","CCL2","CXCL8","CXCL10","VEGFA",
            "VEGFB","PDGFB","PDGFA","TGFB1","TGFB2","TGFB3","ARG1","IDO1",
            "IL10","CD274") ,
    
    embryonically_derived_TAM = c("CXCR4","COL1A2","COL3A1","COL6A1","COL5A1","COL5A1","COL10A1","COL17A1","COL18A1",
                                  "TNC","ELN","SPARC","HAS2","HAS3","LOX","LOXL1","ADAMTS12","MMP2","YAP1"), #exhibit a pro-fibrotic transcriptional profile.
    HSC_D_TAM = c("FLT3","IL12A","IL4","CCL17","IFNB1",
                  "CIITA","ERAP1","PSME1",
                  "HLA-DRB5", "HLA-DPA1","HLA-DMA",
                  "HLA-DRB1","HLA-DQB1", "HLA-DRA", "HLA-DQA1", "HLA-DRB3", "HLA-DPB1","HLA-DRB4"), ### 再添加一点 HLA-DR
    antigen_presentation = c("CIITA","ERAP1","PSME1",  "HLA-DRB5", "HLA-DPA1","HLA-DMA",
                             "HLA-DRB1","HLA-DQB1", "HLA-DRA", "HLA-DQA1", "HLA-DRB3", "HLA-DPB1","HLA-DRB4"),
    Fibrosis = c("TGFB2","ICAM1","IL6R","TGFA","TGFBR1"),
    ECM = c("COL1A2","COL3A1","COL6A1","COL5A1","COL5A1","COL10A1","COL17A1","COL18A1",
            "TNC","ELN","SPARC","HAS2","HAS3","LOX","LOXL1","ADAMTS12","MMP2","YAP1"),
    Metastasis = c("CCND1","WNT5A","PIK3R6","IL6R",'PIK3R1'),
    MMP = c("MMP12","MMP9","MMP2","MMP14","MMP25"),
    Angiogenesis = c("ADM","HIF1A","VEGFA","FLT1","PLAU"),
    LYVEhigh_TRM = c("LYVE1","SIGLEC1","FOLR2","MRC1","MSR1","CD209","PDGFB","PDGFC","CCL8","CCL2","F13A1"), ###PDFGC 写错一个 #37563309
    # profibrosis= c("COL1A2","COL3A1","COL6A1","COL5A1","COL5A1","COL10A1","COL17A1","COL18A1",
    #                "TNC","ELN","SPARC","HAS2","HAS3","LOX","LOXL1","ADAMTS12","MMP2","YAP1"), #exhibit a pro-fibrotic transcriptional profile. "CXCR4",
    ECM_and_ECM_remodeling= c("COL1A2","COL3A1","COL6A1","COL5A1","COL4A4","COL10A1","COL17A1","COL18A1", # 28813661 
                              "TNC","ELN","SPARC","HAS2","HAS3","LOX","LOXL1","ADAMTS12","MMP2","YAP1"),
    collogen_gene = c("COL1A1","COL4A3","COL5A1","COL11A1","COL8A1","COL8A2","COL22A1"),
    M1_M2 = c("ARG1", "ARG2", "CCL13", "CCL17", "CCL18", "CCL20", "CCL22", "CCL24", "CCL4", 
				"CCL5", "CD163", "CD200R1", "CD274", "CD276", "CD32", "CD40", "CD80", "CD86", 
				"CLEC7A", "CTSA", "CTSB", "CTSC", "CTSD", "CXCL10", "CXCL11", "CXCL9", "EGF", 
				"FASL", "FCER1G", "FCGR1", "FCGR2B", "FN1", "H2-AB1", "IDO1", "IL10", "IL12A", 
				"IL12B", "IL1A", "IL1B", "IL1R2", "IL1RN", "IL23A", "IL4RA", "IL6", 
				"IRF1", "IRF4", "IRF5", "KYNU", "LYVE1", "MARCO", "MMP14", "MMP19", "MMP9", 
				"MRC1", "MSR1", "NOS2", "PDCD1LG2", "TGFB1", "TGFB2", "TGFB3", "TNF", "TNFSF12", 
				"TNFSF8", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "VTCN1", "WNT7B"),  #35613577

Angiogenesis_cancercell= c("CCND2", "CCNE1", "CD44", "CXCR4", "E2F3", "EDN1", 
"EZH2", "FGF18", "FGFR1", "FYN", "HEY1", "ITGAV",
 "JAG1", "JAG2", "MMP9", "NOTCH1", "PDGFA", "PTK2", 
 "SPP1", "STC1", "TNFAIP6", "TYMP", "VAV2", "VCAN", "VEGFA"), #33545035


Phagocytosis_cancercell=c("MRC1", "CD163", "MERTK", "C1QB"), #33545035

M1_cancercell=c("IL23", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD86",
 "IL1A", "IL1B", "IL6", "CCL5", "IRF5", "IRF1",
  "CD40", "IDO1", "KYNU", "CCR7"), #33545035

M2_cancercell= c("IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", "CCL24",
"LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", "CTSA", "CTSB",
"CTSC", "CTSD", "TGFB1", "TGFB2", "TGFB3", "MMP14", "MMP19", "MMP9",
"CLEC7A", "WNT7B", "FASL", "TNFSF12", "TNFSF8", "CD276", "VTCN1",
"MSR1", "FN1", "IRF4")  #33545035
    
    
  )     
  


#### SCORE ####
    table(totalcelldata$SubCellTypeup_CN11UMCN2)
    table(totalcelldata$SubCellTypeup_CN11UMCN2)
    ncol(totalcelldata)
    Idents(totalcelldata) <- totalcelldata$SubCellTypeup_CN11UMCN2
    themefamily <- 
      theme(axis.text.x = element_text( color = "#000000"),
            axis.text.y = element_text( color = "#000000"))
    
    
    scoregenenamelist <- names(scorelist)
    genelist_filt0 = lapply(scoregenenamelist,function(genesetname){
      geneset <- scorelist[[genesetname]]
      geneset <- geneset[geneset%in%rownames(totalcelldata)]
      geneset <- unique(geneset)
      return(geneset)
    })
    scoregenenamelist
    names(genelist_filt0) <- scoregenenamelist
    names(genelist_filt0)
    
    scoregenenamelist <- names(genelist_filt0)
    genelistkegg_filt = lapply(scoregenenamelist,function(genesetname){
      geneset <- genelist_filt0[[genesetname]]
      geneexpr <-GetAssayData(object = totalcelldata[["RNA"]], slot = "data")[geneset,]
      class(geneexpr)
      geneexpr <- t(geneexpr)
      class(geneexpr)
      geneexpr <- as.data.frame(geneexpr)
      class(geneexpr)
      head(geneexpr)
      max_values <- apply(geneexpr, 2, max)
      max_values[1]
      max_values_new <- max_values[which(max_values > 0)] 
      names(max_values_new)
      geneset <- geneset[geneset%in%names(max_values_new)]
      geneset <- unique(geneset)
      return(geneset)
    })
    scoregenenamelist
    names(genelistkegg_filt) <- scoregenenamelist
    names(genelistkegg_filt)
    
    for (genesetname in scoregenenamelist) {
      featureset <- list(genelistkegg_filt[[genesetname]])
      print(genesetname)
      print(featureset)
      totalcelldata <- AddModuleScore(
        totalcelldata,
        featureset,
        name = genesetname
      )
    }
    #### metadata export ####
    metadata  <-  totalcelldata@meta.data 
    colnames(metadata) 
    saveRDS(metadata,paste0(DataDir,samplename,"_score_metadata_240527.rds"))
    table(totalcelldata$SampleStage2)
    reduction <- "umap"
    for (genesetname in scoregenenamelist) {
      #if(!dir.exists(paste0(FigDir,""))){dir.create(paste0(FigDir), recursive = TRUE)}
      if(!dir.exists(paste0(FigDir,"score_240527/"))){dir.create(paste0(FigDir,"score_240527/"), recursive = TRUE)}
      
      print(colnames(totalcelldata@meta.data))
      length(colnames(totalcelldata@meta.data))
      Idents(totalcelldata) <- paste0("SubCellTypeup_CN11UMCN2") 
      
      locatedir <- paste0(FigDir,"score_240527/")
      if(!is.null(dev.list())){dev.off()}
      pdf(file = paste0(locatedir, samplename,"_",genesetname,"_totalscore-1.pdf"), width = 5, height = 8)
      
      plot <- VlnPlot(totalcelldata, features = paste0(genesetname,"1"), group.by = "SubCellType0_CN11UMCN2") +
        themefamily+ NoLegend() +
        labs(title = genesetname) 
      print(plot)
      plot <- VlnPlot(totalcelldata, features = paste0(genesetname,"1"), group.by = "SubCellTypeup_CN11UMCN2") +
        themefamily+ NoLegend() +
        labs(title = genesetname) 
      print(plot)
      # ggsave(paste0(FigDir, genesetname,"_VlnPlot.pdf"),width = 10, height = 6)
      
      
      plot <- VlnPlot(totalcelldata, features = paste0(genesetname,"1"), group.by = 'SampleStage2' ) +
        themefamily+ NoLegend() +
        labs(title = genesetname) 
      print(plot)
      dev.off()
      
      if(!is.null(dev.list())){dev.off()}
      pdf(file = paste0(locatedir, samplename,"_",genesetname,"_totalscore-2.pdf"), width = 12, height = 8)
      genesublistVlnplot<-VlnPlot(totalcelldata, 
                                  features =  paste0(genesetname,"1"), 
                                  group.by = "SubCellTypeup_CN11UMCN2", 
                                  pt.size = 0.1, raster=FALSE
      ) + 
        geom_boxplot(width=.2,col="black",fill="white") + NoLegend() +
        labs(title = genesetname) +themefamily+
        theme(plot.title = element_text(hjust = 0.5), 
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x=element_text(angle = 45,hjust = 1)) 
      
      print(genesublistVlnplot)
      
      genesublistVlnplot<-VlnPlot(totalcelldata, 
                                  features =  paste0(genesetname,"1"), 
                                  group.by = "SubCellTypeup_CN11UMCN2", 
                                  pt.size = 0.1, raster=FALSE
      ) + 
        geom_boxplot(width=.2,col="black",fill="white") + NoLegend() +
        labs(title = genesetname) +themefamily+
        theme(plot.title = element_text(hjust = 0.5), 
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.x=element_text(angle = 45,hjust = 1)) 
      
      print(genesublistVlnplot)
      
      
      
      genesublistFeatureplot <- FeaturePlot(totalcelldata, features = paste0(genesetname,"1") , pt.size = 0.1, raster=FALSE)  +
        coord_fixed(ratio=1/1) +themefamily+
        theme(plot.subtitle = element_text(hjust = 0.5))  + labs(subtitle = paste0(samplename,"_",method, "_",  genesetname, "_genes"))                         
      print(genesublistFeatureplot)
      
      dev.off()
      if(!is.null(dev.list())){dev.off()}
      mysubtitle <-  labs(subtitle = paste0(samplename))  #+ theme(plot.subtitle = element_text(hjust = 0.5))
      pdf(file = paste0(locatedir,samplename, "_",  genesetname, "_1.pdf"), width = 20, height = 12)
      total_dimplot <- DimPlot(object = totalcelldata, reduction = reduction, label = TRUE, group.by = "SubCellTypeup_CN11UMCN2", pt.size = 0.2, raster=FALSE) + 
        coord_fixed(ratio=1/1) + mysubtitle + 
        theme(plot.subtitle = element_text(hjust = 0.5)) + labs(title = NULL) +
        theme(legend.title=element_blank(), legend.key.size = unit(40, "pt"), legend.text=element_text(size=20))
      print(total_dimplot)
      rm(total_dimplot)
      if(length(genelistkegg_filt[[genesetname]]) > 1 ){
        genesublistVlnplot<-VlnPlot(totalcelldata, features = genelistkegg_filt[[genesetname]] , group.by = "SubCellTypeup_CN11UMCN2", pt.size = 0.2, raster=FALSE,
                                    stack = TRUE, flip = T) + geom_boxplot(width=.2,col="black",fill="white")+ NoLegend() +
          xlab(NULL) + mysubtitle + labs(title = paste0(genesetname, "_genes")) + 
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 
        print(genesublistVlnplot)
        rm(genesublistVlnplot)
      }
      dev.off()
      
      if(!is.null(dev.list())){dev.off()}
      mysubtitle <-  labs(subtitle = paste0(samplename))  #+ theme(plot.subtitle = element_text(hjust = 0.5))
      pdf(file = paste0(locatedir,samplename, "_",  genesetname, "_2.pdf"), width = 20, height = 12)
      for (i in 1:length(genelistkegg_filt[[genesetname]])){
        genesublistFeatureplot <- FeaturePlot(totalcelldata, features = genelistkegg_filt[[genesetname]][i] , pt.size = 0.2, raster=FALSE)  +
          coord_fixed(ratio=1/1) +
          theme(plot.subtitle = element_text(hjust = 0.5))  + labs(subtitle = paste0(samplename,"_", genesetname, "_genes"))                         
        print(genesublistFeatureplot)
        rm(genesublistFeatureplot)
      }
      dev.off()
      
    } 
    
    table(totalcelldata$SubCellTypeup_CN11UMCN2)
    

    
    #### score vlnplot SubCellTypeup_CN11UMCN2 ####
 
    
        
        medata  <-  totalcelldata@meta.data 
        
    medata <- medata %>% mutate(M2_M11 = M21 - M11)
    colnames(medata)
    colnames(medata)[which(colnames(medata) == "antigen_presentation1")] <- "Antigen_presentation_score"
    colnames(medata)[which(colnames(medata) == "M11")] <- "M1_score"
    colnames(medata)[which(colnames(medata) == "M21")] <- "M2_score"
    colnames(medata)[which(colnames(medata) == "M2_M11")] <- "M2_M1_score"
    #colnames(medata)[which(colnames(medata) == "Angiogenesis1")] <- "Angiogenesis_score"
    colnames(medata)[which(colnames(medata) == "LYVEhigh_TRM1")] <- "LYVEhigh_TRM_score"
    colnames(medata)[which(colnames(medata) == "ECM_and_ECM_remodeling1")] <- "ECM_and_ECM_remodeling_score"
    colnames(medata)[which(colnames(medata) == "Angiogenesis_cancercell1")] <- "Angiogenesis_score"
    #colnames(medata)[which(colnames(medata) == "Phagocytosis_cancercell1")] <- "Phagocytosis_cancercell_score"
    colnames(medata)[which(colnames(medata) == "M1_cancercell1")] <- "M1_cancercell_score"
    colnames(medata)[which(colnames(medata) == "M2_cancercell1")] <- "M2_cancercell_score"
    #colnames(medata)[which(colnames(medata) == "M1_M21")] <- "M1_M2_score"
    
    saveRDS(medata, paste0(DataDir,"Mac_IPMClineage20sample_score_metadata.rds"))
    medata <- readRDS(paste0(DataDir,"Mac_IPMClineage20sample_score_metadata.rds"))
    scorelevel <- c("Antigen_presentation_score","M1_score","M2_score","M2_M1_score",
                    "Angiogenesis_score","LYVEhigh_TRM_score",
                    "ECM_and_ECM_remodeling_score"
                   #"Phagocytosis_cancercell_score",
                   # "M1_cancercell_score", "M2_cancercell_score", "M1_M2_score"
                   )
    medata$SubCellTypeup_CN11UMCN2 <- factor(medata$SubCellTypeup_CN11UMCN2, levels = names(mycolor))
    
    #### source data ####
    colnames(medata)
    #view(df)
    sourcedata <- medata
    sourcedata <- sourcedata %>% dplyr::select(orig.ident, SubCellTypeup_CN11UMCN2, LYVEhigh_TRM_score, ECM_and_ECM_remodeling_score)
    view(sourcedata)
    saveRDS(sourcedata, paste0(DataDir,"v1_figS4f_sourcedata.rds"))
   write.csv(sourcedata, paste0("I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results/V1-240423/Fibroblast/plot-V2-241106/sourcedata/v1_figS4f_sourcedata.csv"))
    
    
    comparelevel <- list(
      c("Mac_FOLR2", "Mac_SPP1"))
    
#### 画 macrophage subpopulation ####    
    themefamily <- theme(panel.grid = element_blank(),
                         axis.text.x = element_text(angle = 45, hjust = 1, color = "#000000"),
                         #axis.text.x = element_text(angle = 90, hjust = 1, color = "#000000"),
                         axis.text.y =  element_text(angle = 0, hjust = 0.5, color = "#000000"),
                         plot.title = element_text(hjust = 0.5,color = "#000000"),
                         axis.title.y =element_text( color="#000000",hjust=0.5),
                         axis.title.x =element_text( color="#000000",hjust=0.5),
                         panel.border = element_rect(color = "#000000", size = 0.5, linetype="solid", fill = NA),
                         legend.position = "none",
                         strip.background = element_rect(
                           fill="#efebbe", size= NULL, linetype= NULL, color =  NULL),
                         plot.margin = margin(t = 5, r = 10, b = 0, l = 10, unit = 'pt'),
                         strip.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5),
                                                     size = 10,  hjust = 0.5, color = "#000000")
    )
    
    
    plistscore_folr <- list()
    
    for (scorename in scorelevel) {
      if(!dir.exists(paste0(FigDir,"score_240527/"))){dir.create(paste0(FigDir,"score_240527/"), recursive = TRUE)}
      #scorename <- "Angiogenesis_score"
      #scorename <- "normal_fibroblast_pancreas"
      p4 <- ggplot(data = medata, aes_string(x = "SubCellTypeup_CN11UMCN2", y = scorename)) +
        geom_point(position = 'jitter', color = 'grey',
                   size = 0.1, alpha = 0.8) +
        geom_violin(aes(fill = SubCellTypeup_CN11UMCN2),
                    color = 'grey', alpha = 1,
                    scale = 'width',
                    linewidth = 0.5, #外轮廓粗细
                    trim = TRUE) +
        geom_boxplot(color = 'white',
                     outlier.color = 'grey',
                     outlier.size = 0.1,
                     width = 0.35, #箱子宽度
                     size = 0.5, #外轮廓描边粗细
                     fill = NA) +
        theme_bw()+
        themefamily + theme(legend.position = "none") +
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_blank()) +
        ggtitle(paste(scorename)) +
        scale_fill_manual(values = mycolor) +
        geom_hline(aes(yintercept=0), colour="#000000", linetype="dashed",size = 0.2) +
        
        scale_color_manual(values = mycolor) + #+RotatedAxis() +
        geom_signif(aes_string(x = "SubCellTypeup_CN11UMCN2", y = scorename),
                    comparisons = comparelevel,
                    map_signif_level=T, # T显示显著性，F显示p值
                    #tip_length=rep(0.01, each = length(comparelevel)), # 修改显著性线两端的长短
                    #y_position = y_positionlist, # 设置显著性线的位置高度
                    size = 0.3, # 修改线的粗细
                    textsize =5, # 修改显著性标记的大小
                    extend_line = 0, # 线的长短：默认为0；
                    vjust=0.5,
                    color="#000000",
                    test = "wilcox.test"#,
        ) # 检验的类型
      #theme_classic() + 
      #coord_flip() 
      p4
      plistscore_folr[[scorename]] <- p4
      ggsave(paste0(FigDir,"score_240527/",scorename,"_vlnplot2.pdf"),width = 2, height = 3)
      ggsave(paste0(FigDir,"score_240527/",scorename,"_vlnplot2.tiff"),width = 2, height = 3)
      
      
    }
    
    plot_grid(plistscore_folr$Antigen_presentation_score,plistscore_folr$M1_score,
              plistscore_folr$M2_score,plistscore_folr$M2_M1_score,
              plistscore_folr$Angiogenesis_score, plistscore_folr$LYVEhigh_TRM_score, 
              plistscore_folr$ECM_and_ECM_remodeling_score,
               align = "vh",  ncol = 4)
    
    FigDir
    
    ggsave(paste0(FigDir,"score_240527/totalscore_vlnplot2_sig.pdf"),width = 8, height = 6)
    ggsave(paste0(FigDir,"score_240527/totalscore_vlnplot2_sig.tiff"),width = 8, height = 6)
    
    
    plot_grid(plistscore_folr$Antigen_presentation_score,plistscore_folr$M1_M2_score,
              plistscore_folr$M1_score,plistscore_folr$M1_cancercell_score,
              plistscore_folr$M2_score,plistscore_folr$M2_cancercell_score,
              plistscore_folr$Angiogenesis_score,  plistscore_folr$Angiogenesis_cancercell_score, 
              plistscore_folr$LYVEhigh_TRM_score, plistscore_folr$ECM_and_ECM_remodeling_score,
               align = "vh",  ncol = 2)
    
    FigDir
    
    ggsave(paste0(FigDir,"score_240527/totalscore_vlnplot2_sig_2.pdf"),width = 4, height = 12)
    ggsave(paste0(FigDir,"score_240527/totalscore_vlnplot2_sig_2.tiff"),width = 4, height = 12)
    

    
    
####-----------------------------------------    
    scorelevel <- c("Antigen_presentation_score","M1_score","M2_score","M2_M1_score",
                    "Angiogenesis_score","LYVEhigh_TRM_score",
                    "ECM_and_ECM_remodeling_score"
                    #"Phagocytosis_cancercell_score",
                    # "M1_cancercell_score", "M2_cancercell_score", "M1_M2_score"
    )
    
  
    
    
    ####宽数据转成长数据 #####
    colnames(medata)
    medata$cellbarcode <- rownames(medata)
    saveRDS(medata, paste0(DataDir,"macrophage_score_data_250704.rds"))
    medata <- readRDS(paste0(DataDir,"macrophage_score_data_250704.rds"))
    medatalong <- medata %>% 
      pivot_longer(cols = scorelevel,
                   names_to = "Mac_score_type",
                   values_to = "score",
                   values_drop_na = T
      )
    
    colnames(medatalong)
    unique(medatalong$Mac_score_type)
    
    #### source data ####
    colnames(medatalong)
    #view(df)
    sourcedata <- medatalong
    sourcedata <- sourcedata %>% dplyr::select(cellbarcode,orig.ident, SampleStage2,  Mac_score_type,score)
    view(sourcedata)
    saveRDS(sourcedata, paste0(DataDir,"v1_fig1k_s4g_sourcedata.rds"))
    write.csv(sourcedata, paste0("I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results/V1-240423/Fibroblast/plot-V2-241106/sourcedata/v1_fig1k_s4g_sourcedata.csv"))
    
    scorename <- "score"
    

    
    themefamily <- theme(panel.grid = element_blank(),
                         axis.text.x = element_text(angle = 0, hjust = 0.5,  color = "#000000",size = 8),
                         #axis.text.x = element_text(angle = 90, hjust = 1, color = "#000000"),
                         axis.text.y =  element_text(angle = 0, hjust = 0.5, color = "#000000",size = 8),
                         plot.title = element_text(hjust = 0.5,color = "#000000"),
                         axis.title.y =element_text( color="#000000",hjust=0.5,size = 10),
                         axis.title.x =element_text( color="#000000",hjust=0.5,size = 10),
                         panel.border = element_rect(color = "#000000", size = 0.5, linetype="solid", fill = NA),
                         legend.position = "none",
                         strip.background = element_rect(
                           fill="#efebbe", size= NULL, linetype= NULL, color =  NULL),
                         plot.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = 'pt'),
                         strip.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5),
                                                     size = 10,  hjust = 0.5, color = "#000000")
    )
    
    
#### Mac_subset vlnplot ####
    
    mycolor <- c( 'Mac_FOLR2' = "#ffadaa",
                  'Mac_SPP1' = "#ffc995"
    )
    
    comparelevel <- list(c("Mac_FOLR2",'Mac_SPP1'))
    
    positionfunc <- function(maxnum){
      positionlist <- c(maxnum
      )
      return(positionlist)
    }
    medatalong$SubCellTypeup_CN11UMCN2 <- factor(medatalong$SubCellTypeup_CN11UMCN2,
                                                 levels = c('Mac_SPP1', "Mac_FOLR2"))
    plistscore_folr <- list()
    for (scorelevelname in scorelevel) {
      #scorelevelname <- "naive_score"
      medatalong1 <- medatalong %>% filter(Mac_score_type %in% scorelevelname)
      maxnum <-  max(medatalong1[,scorename])
      y_positionlist <- positionfunc(maxnum)
      p6 <- ggplot(data = medatalong1, aes_string(x = "SubCellTypeup_CN11UMCN2", y = scorename)) +
        facet_wrap(vars(Mac_score_type), ncol = 1, scales = "free_y") + 
        geom_point(position = 'jitter', color = 'grey',
                   size = 0.1, alpha = 0.8) +
        geom_violin(aes(fill = SubCellTypeup_CN11UMCN2),
                    color = 'grey', alpha = 0.8,
                    scale = 'width',
                    linewidth = 0.6, #外轮廓粗细
                    trim = TRUE) +
        geom_boxplot(color = 'white',
                     outlier.color = NA,
                     width = 0.4, #箱子宽度
                     size = 0.2, #外轮廓描边粗细
                     fill = NA) +
        theme_bw()+
        geom_hline(aes(yintercept=0), colour="#000000", linetype="dashed", size = 0.2) +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 0, hjust = 0.5,  color = "#000000",size = 6),
              #axis.text.x = element_text(angle = 90, hjust = 1, color = "#000000"),
              axis.text.y =  element_text(angle = 0, hjust = 0.5, color = "#000000",size = 6),
              plot.title = element_text(hjust = 0.5,color = "#000000"),
              axis.title.y =element_text( color="#000000",hjust=0.5,size = 6),
              axis.title.x =element_text( color="#000000",hjust=0.5,size = 6),
              panel.border = element_rect(color = "#000000", size = 0.5, linetype="solid", fill = NA),
              legend.position = "none",
              strip.background = element_rect(
                fill="#efebbe"),
              plot.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = 'pt'),
              strip.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5),
                                          size = 8,  hjust = 0.5, color = "#000000")) + 
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_blank()) +
        #ylab("CD4+ T cells")+ 
        theme(  legend.position = "none")+
        #strip.background = element_rect(
        # fill="#efebbe", size=1, linetype="solid", color = "#")) +
        #ggtitle(paste("CD4Ts")) +
        scale_fill_manual(values = mycolor) +
        scale_color_manual(values = mycolor) +# + #+RotatedAxis() +
        #theme_classic() + 
        #差异显著性检验：
        geom_signif(aes_string(x = "SubCellTypeup_CN11UMCN2", y = scorename),
                    comparisons = comparelevel,
                    map_signif_level=F, # T显示显著性，F显示p值
                    #tip_length=rep(0.01, each = length(comparelevel)), # 修改显著性线两端的长短
                    y_position = y_positionlist, # 设置显著性线的位置高度
                    size = 0.3, # 修改线的粗细
                    textsize = 3, # 修改显著性标记的大小
                    extend_line = 0, # 线的长短：默认为0；
                    vjust=0.5,
                    color="#000000",
                    test = "wilcox.test"#,
        ) + # 检验的类型 
        coord_flip() 
      p6
      
      plistscore_folr[[scorelevelname]] <- p6
      
      #plistscore_folr[[scorelevelname]] <- p6
      ggsave(paste0(FigDir,"Mac_score_",scorelevelname,"_vlnplot_SubCellTypeup_CN11UMCN2_sig_pvalue.pdf"),width = 2, height = 2.5)
      ggsave(paste0(FigDir,"Mac_score_",scorelevelname,"_vlnplot_SubCellTypeup_CN11UMCN2_sig_pvalue.tiff"),width = 2, height = 1.5)
      
    }    
    
    
    plot_grid(plistscore_folr$Antigen_presentation_score,plistscore_folr$M1_score,
              plistscore_folr$M2_score,plistscore_folr$M2_M1_score,
              plistscore_folr$LYVEhigh_TRM_score, 
              plistscore_folr$ECM_and_ECM_remodeling_score,
              align = "vh",  ncol = 1,scale = 0.99)
    ggsave(paste0(FigDir,"score_240527/totalscore_vlnplot_SubCellTypeup_CN11UMCN2_sig_horizon.pdf"),width = 2, height = 8*0.8)
    ggsave(paste0(FigDir,"score_240527/totalscore_vlnplot_SubCellTypeup_CN11UMCN2_sig_horizon.tiff"),width = 2, height = 8*0.8)
    
#### samplestage2 vlnplot ####  
    mycolor <- c("UNIN" = "#5ad9a6", 
                 # "Gastric" ="#bebada" ,
                 #"Pancreatobiliary" ="#61a0d5" ,
                 "IPMN" ="#2f6490" ,
                 "JUNC" = "#5d7091", 
                 "PDAC" = "#f6bd16", 
                 "cPDAC" ="#e86949", 
                 "PASC" = "#ff98c3")
    
    comparelevel <- list(
      c("UNIN", "IPMN"), c("IPMN", "PDAC"),
      c("UNIN", "PDAC"))
    
    
    color <- c("UNIN" = "#5ad9a6", 
               # "Gastric" ="#bebada" ,
               #"Pancreatobiliary" ="#61a0d5" ,
               "IPMN" ="#2f6490" ,
               "JUNC" = "#5d7091", 
               "PDAC" = "#f6bd16", 
               "cPDAC" ="#e86949", 
               "PASC" = "#ff98c3")
    
    positionfunc <- function(maxnum){
      positionlist <- c(maxnum,maxnum,
                        maxnum*1.08
      )
      return(positionlist)
    }
    
    
    
    library(rstatix)
    library(ggsignif)
    library(ggpubr)
    unique(medatalong$SampleStage2)
    #medatalong$SampleStage2 <- factor(medatalong$SampleStage2,levels = c("UNIN","IPMN","PDAC"))
    medatalong$SampleStage2 <- factor(medatalong$SampleStage2,levels = c("PDAC","IPMN","UNIN"))
    table(medatalong$SampleStage2)
    table(medatalong$SubCellTypeup_CN11UMCN2)
    
    table(medatalong$SampleStage2)
    plist <- list()
    for (scorelevelname in scorelevel) {
      #scorelevelname <- "M2_M1_score"
      medatalong1 <- medatalong %>% filter(Mac_score_type %in% scorelevelname)
      table(medatalong1$SampleStage2)
      table(medatalong1$SubCellTypeup_CN11UMCN2)
      maxnum <-  max(medatalong1[,scorename])
      y_positionlist <- positionfunc(maxnum)
      p6 <- ggplot(data = medatalong1, aes_string(x = "SampleStage2", y = scorename)) +
        facet_wrap(vars(Mac_score_type), ncol = 1, scales = "free_y") + 
        geom_point(position = 'jitter', color = 'grey',
                   size = 0.1, alpha = 0.8) +
        geom_violin(aes(fill = SampleStage2),
                    color = 'grey', alpha = 0.8,
                    scale = 'width',
                    linewidth = 0.6, #外轮廓粗细
                    trim = TRUE) +
        geom_boxplot(color = 'white',
                     outlier.color = NA,
                     width = 0.4, #箱子宽度
                     size = 0.2, #外轮廓描边粗细
                     fill = NA) +
        theme_bw()+
        geom_hline(aes(yintercept=0), colour="#000000", linetype="dashed", size = 0.2) +
        themefamily + 
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_blank()) +
        #ylab("CD4+ T cells")+ 
        theme(  legend.position = "none")+
        #strip.background = element_rect(
        # fill="#efebbe", size=1, linetype="solid", color = "#")) +
        #ggtitle(paste("CD4Ts")) +
        scale_fill_manual(values = mycolor) +
        scale_color_manual(values = mycolor) +# + #+RotatedAxis() +
        #theme_classic() + 
               #差异显著性检验：
        geom_signif(aes_string(x = "SampleStage2", y = scorename),
                    comparisons = comparelevel,
                    map_signif_level=F, # T显示显著性，F显示p值
                    #tip_length=rep(0.01, each = length(comparelevel)), # 修改显著性线两端的长短
                    y_position = y_positionlist, # 设置显著性线的位置高度
                    size = 0.3, # 修改线的粗细
                    textsize = 3, # 修改显著性标记的大小
                    extend_line = 0, # 线的长短：默认为0；
                    vjust=0.5,
                    color="#000000",
                    test = "wilcox.test"#,
        ) + # 检验的类型 
      coord_flip() 
      p6
      
      plist[[scorelevelname]] <- p6
      
      #plist[[scorelevelname]] <- p6
      ggsave(paste0(FigDir,"Mac_score_",scorelevelname,"_vlnplot_Samplestage2_sig_pvalue.pdf"),width = 4, height = 2)
      ggsave(paste0(FigDir,"Mac_score_",scorelevelname,"_vlnplot_Samplestage2_sig_pvalue.tiff"),width = 4, height = 2)
      
    }    
    
    plot_grid(plist$Antigen_presentation_score,plist$M1_score,
              plist$M2_score,plist$M2_M1_score,
               plist$LYVEhigh_TRM_score, 
              plist$ECM_and_ECM_remodeling_score,plist$Angiogenesis_score,
              align = "vh",  ncol = 2,scale = 0.99)
    
    FigDir
    
    ggsave(paste0(FigDir,"score_240527/totalscore__vlnplot_Samplestage2_sig_horizon.pdf"),width = 2, height = 8)
    ggsave(paste0(FigDir,"score_240527/totalscore_vlnplot_Samplestage2_sig_horizon.tiff"),width = 2, height = 8)
    
    
    
  P1 <-   plot_grid(plist$M1_score,
              plist$M2_score,plist$M2_M1_score,
              align = "vh",  ncol = 1,scale = 0.99)
  P1
    ggsave(paste0(FigDir,"score_240527/M1-2_vlnplot_Samplestage2_sig_horizon.pdf"),width = 2, height = 3)
    ggsave(paste0(FigDir,"score_240527/M1-2_vlnplot_Samplestage2_sig_horizon.tiff"),width = 2, height =3)
    
    P2 <-   plot_grid(plist$LYVEhigh_TRM_score,
              plist$ECM_and_ECM_remodeling_score,
              align = "vh",  ncol = 1,scale = 0.99)
    P2
    ggsave(paste0(FigDir,"score_240527/LYVEhigh_TRM_vlnplot_Samplestage2_sig_horizon.pdf"),width = 2, height = 2)
    ggsave(paste0(FigDir,"score_240527/LYVEhigh_TRM_vlnplot_Samplestage2_sig_horizon.tiff"),width = 2, height =2)
    
    
    
    
#### dotplot ####
    #### dotplot ####
    
    
    mycolor = c( 'Mac_FOLR2' = "#ffadaa",
                 'Mac_SPP1' = "#ffc995"
    )
    genelist <- 
      c("LYVE1","SIGLEC1","FOLR2","MRC1","MSR1","CD209","PDGFB","PDGFC","CCL8","CCL2","F13A1"
      )
    
      
    
    
    #### dimplot ####

    
    
    #### IPMC lineage ####
    dotplot <- DotPlot(totalcelldata, features = unique(genelist),group.by = "SubCellTypeup_CN11UMCN2")+ #coord_flip()+
      theme_bw()+
      theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1), axis.text.y = element_text(size = 10))+
      labs(x=NULL,y=NULL)+
      guides(size=guide_legend(order=3))+
      scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
    print(dotplot)
    
    if(!is.null(dev.list())){dev.off()}
    pdf(file = paste0(FigDir,samplename,  "_",  "deg_mannual_LYVE.pdf"), width = 20, height = 14)
    print(dotplot)
    dev.off()
    
    df<- dotplot$data
    df
    saveRDS(df, paste0(DataDir, "_Mac_deg_dotplot_LYVE_df.rds"))
    write.csv(df, paste0(DataDir, "_figs2h_Mac_deg_dotplot_LYVE_df.csv"))

    #### source data ####
    sourcedata <- df
    write.csv(sourcedata, paste0("I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results/V1-240423/Fibroblast/plot-V2-241106/sourcedata/v1_figs2h_sourcedata.csv"))
    
    
    library("RColorBrewer")
    col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
    levels(df$id)
    head(df)
    df$id <- factor(df$id,levels = rev(levels(df$id)))
    p2 <-  ggplot(df, aes(x = id, y = features.plot)) +
      geom_point(aes(color = avg.exp.scaled, size = pct.exp)) +
      
      theme_bw()+
      theme( #panel.grid = element_blank(),
        #axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5, color = "#000000"),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1, color = "#000000",face="italic"),
        axis.text.y = element_text(color = "#000000"),
        aspect.ratio = length(unique(df$id))/length(unique(df$features.plot)))+
      # 坐标轴label：
      xlab(NULL)+
      ylab(NULL)+
      #scale_color_gradientn(colors = c("#BBBBBB","#CC0000"), limits = c(-1.3, 2.5)) +
      guides(size=guide_legend(order=3))+
      scale_color_gradientn(values = seq(0,1,0.2),colours = col) +
      coord_flip()
    p2
    
    
    ggsave(paste0(FigDir, "Mac_deg_maunal_dotplot_LYVE_240527.pdf"), width = length(unique(df$features.plot))/2.5*0.8, height = length(unique(df$id))/1.5*0.8, limitsize = F)
    ggsave(paste0(FigDir, "Mac_deg_maunal_dotplot_LYVE_240527.tiff"), width = length(unique(df$features.plot))/2.5*0.8, height = length(unique(df$id))/1.5*0.8, limitsize = F)
    
    ggsave(paste0(FigDir, "Mac_deg_maunal_dotplot_LYVE_240527_2.pdf"), limitsize = F,  width = length(unique(df$features.plot))/2.5, height = length(unique(df$id)))
    ggsave(paste0(FigDir, "Mac_deg_maunal_dotplot_LYVE_240527_2.tiff"), limitsize = F,  width = length(unique(df$features.plot))/2.5, height = length(unique(df$id)))
    
    ggsave(paste0(FigDir, "Mac_deg_maunal_dotplot_LYVE_240527_3.pdf"), limitsize = F,  width = length(unique(df$features.plot))/2.5+0.4, height = length(unique(df$id)))
    ggsave(paste0(FigDir, "Mac_deg_maunal_dotplot_LYVE_240527_3.tiff"), limitsize = F,  width = length(unique(df$features.plot))/2.5+0.4, height = length(unique(df$id)))
    
    
#### score gene samplestage2 ####
    names(scorelist)
    
    genesetnamelist <- c( "antigen_presentation","M1" ,"M2",
                          "Angiogenesis_cancercell" ,"LYVEhigh_TRM", 
                          "ECM_and_ECM_remodeling","M1_cancercell","M2_cancercell")
    for (genesetname in genesetnamelist) {
      
      #genesetname <- "cytotoxic_score_gut"
      VlnPlot(totalcelldata, features = genelistkegg_filt[[genesetname]] , group.by = "SampleStage2", pt.size = 0.2, raster=FALSE,
              stack = TRUE, flip = T) + geom_boxplot(width=.2,col="black",fill="white")+ NoLegend() +
        xlab(NULL) +  labs(title = paste0(genesetname, "_genes")) + 
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) 
      x=length(genelistkegg_filt[[genesetname]])
      ggsave(paste0(FigDir,"Mac_gene_",genesetname,"_score_SampleStage2_vlnplot.pdf"),width = 4, height = 2+x*0.5)
      ggsave(paste0(FigDir,"Mac_gene_",genesetname,"_score_SampleStage2_vlnplot.tiff"),width = 4, height = 2+x*0.5)
      
    }
    
#### Samplestage2 radar plot ####
    library(devtools)
    library(fmsb) #雷达图绘制包1
    library(ggradar) #雷达图绘制包2
    
    scorelevel <- c("Antigen_presentation_score","M1_score","M2_score",
                    "Angiogenesis_score","LYVEhigh_TRM_score",
                    "ECM_and_ECM_remodeling_score"
                    #"Phagocytosis_cancercell_score",
                    # "M1_cancercell_score", "M2_cancercell_score", "M1_M2_score"
    )
    
    
    
    
    #### SAMPLESTAGE2 ####
    
    table(medata$orig.ident)
    unique(medata$SampleStage2)
    
    summarydf <-  medata %>% dplyr::select(SampleStage2,M1_score, M2_score, Angiogenesis_score,
                                    LYVEhigh_TRM_score,ECM_and_ECM_remodeling_score ) %>% 
      group_by(SampleStage2) %>% 
      summarise_if(is.numeric, mean, na.rm = TRUE) 
    
    
    summarydf <- as.data.frame(summarydf)
    rownames(summarydf) <- summarydf$SampleStage2
    
    summarydf$SampleStage2 <- NULL
    
    library(scales)
    
    #scale(summarydf) #, center = TRUE, scale = TRUE 标准差
    summarymatrix <- as.matrix(summarydf)
    summarymatrix1 <- t(summarymatrix)
    #rownames(summarymatrix) <- rownames(summarydf)
    #scorelevel <-   c("naive_score","immunosuppressive_CAF", "myCAF","iCAF","apCAF")
    aa <- rescale(as.matrix(summarymatrix1), to = c(0, 1)) #### 归一化不行 
    
    
    scorelevel <- c("Antigen_presentation_score","M1_score","M2_score","M2_M1_score",
                    "Angiogenesis_score","LYVEhigh_TRM_score",
                    "ECM_and_ECM_remodeling_score"
                    #"Phagocytosis_cancercell_score",
                    # "M1_cancercell_score", "M2_cancercell_score", "M1_M2_score"
    )
    
    
    #### 归一化到0到1之间 #### 
    summaryscale <-   summarydf %>% 
      dplyr::select(M1_score, M2_score, Angiogenesis_score,
             LYVEhigh_TRM_score,ECM_and_ECM_remodeling_score) %>%
      scale(center=apply(., 2, min),
            scale=apply(., 2, max) - apply(., 2, min))
    
    summaryscaledf <- as.data.frame(summaryscale)
    
    #### source data ####
    colnames(summaryscaledf)
    #view(df)
    sourcedata <- summaryscaledf
    #sourcedata <- sourcedata %>% dplyr::select(orig.ident, SubCellTypeup_CN11UMCN2, LYVEhigh_TRM_score, ECM_and_ECM_remodeling_score)
    view(sourcedata)
    saveRDS(sourcedata, paste0(DataDir,"v1_fig1j_sourcedata.rds"))
    write.csv(sourcedata, paste0("I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results/V1-240423/Fibroblast/plot-V2-241106/sourcedata/v1_fig1j_sourcedata.csv"))
    
    
    #参数优化：
    mycol <- c("UNIN" = "#5ad9a6", 
               #"Gastric" ="#bebada" ,
               #"Pancreatobiliary" ="#61a0d5" ,
               "IPMN" ="#2f6490" ,
               #"proIPMN" = "#5d7091", 
               "PDAC" = "#f6bd16"#, 
               #"cPDAC" ="#e86949", 
               #"PASC" = "#ff98c3"
    )
    
    #mycol <- c4a('superfishel_stone', 5) #自定义配色
    class(mycol) 
    #### ggradar ####
    ggradar(summaryscaledf[1,]#,
            # grid.min = 0,
            # grid.mid = 2.1,
            # grid.max = 4.3 #需要根据数据限制网格线数值范围，若数据不在区间内(默认0-1)会报错
    )
    
    summaryscaledf$SampleStage2 <- rownames(summaryscaledf)
    scorelevel
    summaryscaledf <- summaryscaledf %>% select(SampleStage2,M1_score, M2_score, Angiogenesis_score,
                                                LYVEhigh_TRM_score,ECM_and_ECM_remodeling_score)
    summaryscaledf$SampleStage2 <- factor(summaryscaledf$SampleStage2, levels = c("UNIN","IPMN","PDAC"))
    colnames(summaryscaledf) <- c("SampleStage2","M1_score", "M2_score", "Angiogenesis_score",
                                  "LYVEhigh_TRM_score","ECM_and_ECM_remodeling_score")
    colnames(summaryscaledf)
    ggradar(summaryscaledf#,
            #$grid.min = 0,
            #grid.mid = 2.1,
            #grid.max = 4.3
    )
    mycol <-  rev(mycol)
    radar1 <-  ggradar(summaryscaledf,
                       base.size=12,
                       grid.min = 0, #网格线最小值
                       #grid.mid = 2.1, #网格线均值
                       grid.max = 1, #网格线最大值
                       values.radar = c(0,0.5,1), #轴标签显示
                       group.colours = mycol,
                       group.point.size = 4,#分组点大小
                       group.line.width = 1, #线宽
                       background.circle.colour = 'grey', #背景填充色
                       #axis.label.size = 2,
                       #grid.label.size = 20,
                       background.circle.transparency = 0.1, #背景填充不透明度(这里改为0可去掉背景填充)
                       plot.legend=T,#
                       legend.position = 'right', #图例位置
                       legend.text.size = 12, #图例标签大小
                       fill = TRUE, #各分组是否填充色
                       fill.alpha = 0.3 #分组填充色不透明度
    )+# +geom_richtext()
      theme(plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = 'pt')) 
    radar1
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot.pdf"),width = 8, height = 4)
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot.tiff"),width = 8, height = 4)
    
    
    
    radar2_inside <-  ggradar(summaryscaledf,
                              #base.size=12,
                              axis.labels = rep(NA, 5),
                              grid.min = 0, #网格线最小值
                              #grid.mid = 2.1, #网格线均值
                              grid.max = 1, #网格线最大值
                              values.radar = c(0,0.5,1), #轴标签显示
                              group.colours = mycol,
                              group.point.size = 3,#分组点大小
                              group.line.width = 1, #线宽
                              background.circle.colour = 'grey', #背景填充色
                              #axis.label.size = 2,
                              #grid.label.size = 20,
                              background.circle.transparency = 0.1, #背景填充不透明度(这里改为0可去掉背景填充)
                              plot.legend=T,#
                              legend.position = "none",
                              #legend.position = 'right', #图例位置
                              legend.text.size = 10, #图例标签大小
                              fill = TRUE, #各分组是否填充色
                              fill.alpha = 0.3 #分组填充色不透明度
    )+# +geom_richtext()+
      theme(plot.background = element_blank(),
            panel.background = element_blank())
    # theme(plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = 'pt'))  
    radar2_inside
    
    #### update raderplot ####
    tmp <- data.frame(x = rep(1, 5),
                      y = rep(6, 5),
                      group = c("M1_score", "M2_score", "Angiogenesis_score",
                                "LYVEhigh_TRM_score","ECM_and_ECM_remodeling_score"))
    genename <- c("IL1B","TNF","CD86",
                  "NRP1","WNT5A","ARG1",
                  "VEGFA","VCAN","PDGFA",
                  "LYVE1","SIGLEC1","FOLR2",
                  "COL1A2","COL3A1","SPARC"
    )
    
    p1 <- ggplot()+
      # 圆环：
      geom_bar(data = tmp, aes(x, y, fill = group), stat = "identity", position = "dodge")+
      # 文本注释：
      geom_text(aes(x = rep(1,15), y = rep(2, 15),
                    label = genename, group = 1:15),
                color = "black", size =3,
                position = position_dodge(width = 0.9),
                hjust = 0.5)+
      # geom_text(aes(x, y, label = gsub("[.]", " ", df$group), group = group),
      #           color = "white",
      #           position = position_dodge(width = 0.9))+
      scale_fill_manual(values = c("#999999", "#ebbfbf", "#c0cad9", "#f6dcc2", "#bfdbdc"))+
      #ylim(-5.5,3)+
      #ylim(-8,4)+
      ylim(-15,6)+
      # 0.63是计算得来，5个色块，第一个色块的正中心要对准0的位置，
      # 所以2pi/10=0.628即为第一个色块左边界的位置
      # 所以pi/3=0.628即为第一个色块左边界的位置
      coord_polar(start = -0.628)+
      theme_void()+
      theme(legend.position = "none")
    
    p1
    
    p3 <-   p1 + inset_element(radar2_inside, left = 0, bottom = 0, right = 0.99, top = 0.99)
    p3 
    
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot_2.pdf"),width = 4, height = 4)
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot_2.tiff"),width = 4, height = 4)
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot_3.pdf"),width = 4, height = 4)
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot_3.tiff"),width = 4, height = 4)
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot_4.pdf"),width = 4, height = 4)
    ggsave(paste0(FigDir,"Mac","_score_SampleStage2_razerplot_4.tiff"),width = 4, height = 4)
    
#### heatmap SampleStage2 ####
    library(dplyr)
    library(pheatmap)
    library(cols4all) #no cols4all
    library(ComplexHeatmap)
    library(circlize)
    
   
    scorelist_heatmap = list(
      M1 = c("NOS2","IL12A" ,"IL12B" ,"IL6"  ,"TNF" ,"CD86","IL1B" ,"CXCL9","CXCL10" ,"CXCL11","CXCL12","STAT1","AIM2"), #IL23, VEGFA, MRC1, CD209, PDGFC, 
      M2 = c("MRC1", "ARG1" , "CCL22" , "CCL17" , "CCL24","CCL18","NRP1","VEGFB",
                   "LLGL1","MMP12","CD209","IL10","WNT5A","ADORA3","STAT6","SOCS3","IRF4",
                   "LTA","TFRC","IL4"),
      
      Angiogenesis= c("CCND2", "CCNE1", "CD44", "CXCR4", "E2F3", "EDN1", 
                            "EZH2", "FGF18", "FGFR1", "FYN", "HEY1", "ITGAV",
                            "JAG1", "JAG2", "MMP9", "NOTCH1", "PDGFA", "PTK2", 
                            "SPP1", "STC1", "TNFAIP6", "TYMP", "VAV2", "VCAN", "VEGFA"), #33545035
      
      LYVEhigh_TRM = c("LYVE1","SIGLEC1","FOLR2","MSR1","PDGFB","PDGFC","CCL8","CCL2","F13A1"), ###PDFGC 写错一个 #37563309
      # profibrosis= c("COL1A2","COL3A1","COL6A1","COL5A1","COL5A1","COL10A1","COL17A1","COL18A1",
      #                "TNC","ELN","SPARC","HAS2","HAS3","LOX","LOXL1","ADAMTS12","MMP2","YAP1"), #exhibit a pro-fibrotic transcriptional profile. "CXCR4",
      ECM_remodeling= c("COL1A2","COL3A1","COL6A1","COL5A1","COL4A4","COL10A1","COL17A1","COL18A1", # 28813661 
                              "TNC","ELN","SPARC","HAS2","HAS3","LOX","LOXL1","ADAMTS12","MMP2","YAP1")
      
      
      
      
    )    
  
    genescore_df <-  as.data.frame(unlist(scorelist_heatmap))
   
   soretypen <- c(rep("M1",13),
                  rep("M2",20),
                  rep("Angiogenesis",25),
                  rep("LYVEhigh_TRM",9),
                  rep("ECM_remodeling",18))
   length(soretypen)
   genename <- genescore_df$`unlist(scorelist_heatmap)`
   length(genename)
    genescoretype1 <- as.matrix(soretypen)
    
    # genename <- unlist(c(scorelist$M1, scorelist$M2,scorelist$Angiogenesis_cancercell,
    #               scorelist$LYVEhigh_TRM,scorelist$ECM_and_ECM_remodeling))
    length(genename)

    length(soretypen)
      length(scorelist$ECM_and_ECM_remodeling)

    
    
class(genescoretype)
    aver_dt <- AverageExpression(totalcelldata,
                                 features = genename,
                                 group.by = 'SampleStage2',
                                 slot = 'data')
    
    aver_dt$RNA
    aver_dtt <-scale(t(aver_dt$RNA))
    view(aver_dtt)
    #### source data ####
    colnames(aver_dtt)
    #view(df)
    sourcedata <- aver_dtt
    #sourcedata <- sourcedata %>% dplyr::select(orig.ident, SubCellTypeup_CN11UMCN2, LYVEhigh_TRM_score, ECM_and_ECM_remodeling_score)
    view(sourcedata)
    saveRDS(sourcedata, paste0(DataDir,"v1_figs4k_sourcedata.rds"))
    write.csv(sourcedata, paste0("I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/Results/V1-240423/Fibroblast/plot-V2-241106/sourcedata/v1_figs4k_sourcedata.csv"))
    
    
####     
    aver_dtt2 <- t(scale(t(aver_dt$RNA)))
    
    splitlabel <- factor(as.character(genescoretype1), levels = c("M1","M2","Angiogenesis",
                                                          "LYVEhigh_TRM","ECM_remodeling"))
  
    library(ComplexHeatmap)
  pheat <-   Heatmap(
      aver_dtt, name = 'Expression', 
      col = colorRampPalette(c('#240EFF', '#B491F9', '#ECEAEF', '#FFA48C', '#FF1A09'))(100),  #基因表达的渐变色设置
      row_names_side = "left",
      cluster_rows = F, cluster_columns = T,  #行列根据基因表达聚类
      row_names_gp = gpar(fontsize = 10), 
      column_names_gp = gpar(fontsize = 10, fontface = "italic"),  #行（基因名）列（样本名）字体大小
      column_names_rot = 45,
      column_split =splitlabel,
      rect_gp = gpar(col = "white", lwd = 0.5),    ##热图矩阵线设置
      border_gp = gpar(col = "black", lty = 5), ##整个边界线设置
      #根据基因类型定义行的颜色
      # right_annotation = HeatmapAnnotation(
      #   Class = gene_anno, which = 'row', show_annotation_name = FALSE, 
      #   col = list(Class= c('Carbohydrate metabolism' = '#FFAD30', 'Lipid metabolism' = '#634FB8', 'Amino acid metabolism' = '#68AD30'))
      # ),
      #根据样本分组定义列的颜色
      top_annotation = HeatmapAnnotation(
        Score = genescoretype1, which = 'column', show_annotation_name = F,
        col = list(Score = c("M1" ="#999999",
                   "M2" = "#ebbfbf",
                   "Angiogenesis" = "#c0cad9",
                   "LYVEhigh_TRM" = "#f6dcc2",
                   "ECM_remodeling" ="#bfdbdc") )
      )
    )
  if(!dir.exists(paste0(FigDir,"Heatmap_gene/"))){dir.create(paste0(FigDir,"Heatmap_gene/"), recursive = TRUE)}
    pdf(file = paste0(FigDir,"Heatmap_gene/","Heatmap_Mac_gene_2.pdf"), width = 17, height =2.2)
    print(pheat)

dev.off()


#### kegg ####

if(!dir.exists(paste0(FigDir,"KEGG/"))){dir.create(paste0(FigDir,"KEGG/"), recursive = TRUE)}

#https://mp.weixin.qq.com/s/8XoLkXLbDEDpyof1c8QrPA
# #### hallmark 基因集 ####
library(clusterProfiler)
library('msigdbr')
msigdbr_species() #列出有的物种
#KEGG 选择基因集合
human_KEGG = msigdbr(species = "Homo sapiens", #物种
                     category = "C2",
                     subcategory = "KEGG") %>% 
  dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol或者ID
human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)#list


#### 方法1 基于细胞类型 不可行 ####
#计算各celltype的表达量均值
Idents(totalcelldata) <- "SubCellTypeup_CN11UMCN2" 
expr <- AverageExpression(totalcelldata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)
head(expr)


library(GSVA)
gsva.kegg <- gsva(expr, gset.idx.list = human_KEGG_Set, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=40)
head(gsva.kegg)



library(limma)
# 构建分组文件
group_list <- data.frame(celltype = colnames(gsva.kegg), group = c(rep("Mac_SPP1", 1), rep("Mac_FOLR2", 1)))

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva.kegg)

# 构建差异比较矩阵
contrast.matrix <- makeContrasts('Mac_SPP1-Mac_FOLR2', levels = design)

summary(gsva.kegg)
# 差异分析，Mac_FOLR2 vs. Mac_SPP1
fit <- lmFit(gsva.kegg, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff)


#设置分组
diff$group <- ifelse( diff$logFC > 0 & diff$adj.P.Val < 0.05 ,"up" ,
                      ifelse(diff$logFC < 0 & diff$adj.P.Val < 0.05 ,"down","noSig")
)


#设置label的位置
diff2 <- diff %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  filter(group != "noSig") %>% #注释掉该行即可
  arrange(t) %>% 
  rownames_to_column("ID")

#通过factor设置展示顺序

diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))

#### 方法2 ####
Mac.data <- totalcelldata[["RNA"]]@data
Mac.data <- as.matrix(Mac.data)



# go C5


gsva <- gsva(Mac.data, human_KEGG_Set,
             kcdf="Gaussian",
             method = "gsva",
             parallel.sz=40
) #BiocManager::install("GSVA") ,method = "ssgsea"

saveRDS(gsva, paste0(DataDir, "Mac_kegg_gsva.rds"))
gsva[1:4,1:4] 


#### 处理gsva结果 ####

#之前定义过分组信息

group_list2 <- totalcelldata@meta.data[,c("SubCellTypeup_CN11UMCN2")]

#标准limma分析
design <- model.matrix(~0+factor(group_list2))
colnames(design)=levels(factor(group_list2))
rownames(design)=colnames(expr)

contrast.matrix<-makeContrasts(contrasts = Mac_FOLR2 - Mac_SPP1,levels = design)

fit <- lmFit(gsva.kegg2,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff2 <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

#### heatmap ####
library(ComplexHeatmap)
library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)



ggplot(diff2, aes(ID, t,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","up"), #设置颜色
                    values = c("#008020","#08519C"))+
  geom_text(data = diff2, aes(label = diff2$ID, #添加通路标签
                              y = diff2$nudge_y),
            nudge_x =0,nudge_y =0,hjust =diff2$hjust,
            size = 3)+ #设置字体大小
  labs(x = "KEGG pathways", #设置标题 和 坐标
       y=paste0("t value of GSVA score\n","celltype-unknown"),
       title = "GSVA")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(), #主题微调
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none" #去掉legend
  )


human_GO = msigdbr(species = "Homo sapiens", #物种
                   category = "C5",
                   subcategory = "GO:BP")

help(msigdbr)
aa <- msigdbr_collections()
aa$gs_cat
aa$gs_subcat
table(totalcelldata$SubCellTypeup_CN11UMCN2)
Idents(totalcelldata) <- "SubCellTypeup_CN11UMCN2"
filtercelldata_markers <- FindAllMarkers(totalcelldata, 
                                         only.pos = FALSE, 
                                         min.pct = 0.25, 
                                         logfc.threshold = 0.25)
filtercelldata_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) -> top10
filtercelldata_markers %>% group_by(cluster) %>% top_n(n=25, wt = avg_log2FC) -> top25
write.csv(top10,file = paste0(DataDir, "Mac_subset_top10.csv")) 
write.csv(top25,file = paste0(DataDir, "Mac_subset_top25.csv")) 
