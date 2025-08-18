#### prepare core&environment ####
rm(list=ls())
sink("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering02_CN1-23samples_PSCT_20pc_1-5res_DFprocessing_230426/Clustering02_CN1-23samples_PSCT_20pc_1-5res_DFprocessing_230426.txt", append=TRUE, split=TRUE )
pdf("/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/bashout/02Clustering/Clustering02_CN1-23samples_PSCT_20pc_1-5res_DFprocessing_230426/Clustering02_CN1-23samples_PSCT_20pc_1-5res_DFprocessing_230426.pdf") #
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
	#suppressPackageStartupMessages(library(clusterProfiler))
	#suppressPackageStartupMessages(library(org.Hs.eg.db))
	options(stringsAsFactors = F)
	filesDir<-'/share/home/sigrid/IPMN_PROGRAM/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/'
	setwd(filesDir)
	getwd()

    if(!dir.exists(paste0(filesDir, "Results/02Clustering/SaveData/DFprocessing/"))){dir.create(paste0(filesDir, "Results/02Clustering/SaveData/DFprocessing/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Figures/DFprocessing/"))){dir.create(paste0(filesDir, "Results/02Clustering/Figures/DFprocessing/"), recursive = TRUE)}
    if(!dir.exists(paste0(filesDir, "Results/02Clustering/Tables/DFprocessing/"))){dir.create(paste0(filesDir, "Results/02Clustering/Tables/DFprocessing/"), recursive = TRUE)}



		CN1_23_filter <- readRDS(paste0(filesDir, "Results/01PreProcess/SaveData/CN1_23_filter_order_230426.rds"))
		# load(paste0(filesDir, "Results/02Clustering/SaveData/CN1_23_filter.RData"))

		#### SCTranform tips ####
		#Before we run this for loop, we know that the output can generate large R objects/variables in terms of memory. If we have a large dataset, then we might need to adjust the limit for allowable object sizes within R (Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:
		#options(future.globals.maxSize = 4000 * 1024^2)
		# SCTncells <- ncol(sc_filter_sampledata)
		# sc_SCT_sampledata <- SCTransform(sc_filter_sampledata, ncells = SCTncells, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)

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

		names(CN1_23filter_SCT)
	#rm(sceList1_SCT1_1_preannotation_persample_DF, sceList1_SCT1_2_preannotation_persample_DF, sceList1_SCT1_3_preannotation_persample_DF)

	#saveRDS(CN1_23filter_SCT, paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_SCT.rds"))

	namelist <- names(CN1_23filter_SCT)

#### PCA&FINDCLUSTER  ####
	#### functions ####
	    ####  functions: RUNPCA step specific npcs   ####   
	        pca_persample <- function(sc_persampledata, samplename, method = c("RNA", "SCT"), npc_num = 50 ){
	            if(method == "SCT"){
	                DefaultAssay(sc_persampledata) <- "SCT"
	            }else{
	                DefaultAssay(sc_persampledata) <- "RNA"
	            }
	            sc_persampledata <- RunPCA(sc_persampledata, npcs = npc_num, verbose = FALSE)
	            return(sc_persampledata)
	        }
	  	#### functions: find clusters step for clustering each sample ####   
	        findclusters_persample <- function(sc_persampledata, samplename, method = c("RNA", "SCT"), npc_num = 50, res_num = seq(0.1, 0.5, by = 0.2)){
	            if(method == "SCT"){
	                DefaultAssay(sc_persampledata) <- "SCT"
	                method <- "SCT"
	            }else{
	                DefaultAssay(sc_persampledata) <- "RNA"
	                method <- "RNA"
	            }
	            sc_persampledata <- FindNeighbors(object = sc_persampledata , dims = 1:npc_num, verbose = FALSE)
	            sc_persampledata <- FindClusters(object = sc_persampledata , resolution = res_num, verbose = FALSE) 
	            sc_persampledata <- RunUMAP(sc_persampledata , dims = 1:npc_num)
	            Idents(sc_persampledata) <- paste0(method,"_snn_res.", as.character(res_num))  
	            return(sc_persampledata) 
	        }


#### total_structure  ####
    namelist <- names(CN1_23filter_SCT)
    namelist

    CN1_23filter_PSCT_20pc.5res_DF = lapply(namelist,function(samplename){
	    #### get the samplename parameters ####
	    	print(samplename)
	    	sc_persampledata <- CN1_23filter_SCT[[samplename]]
	    	print(head(sc_persampledata))
	    	method <- "SCT"
	    	npc_num <- 20
    		res_num <- 0.5

	    #### PCA & findclusters ####
	    	sc_persampledata <- pca_persample(sc_persampledata, samplename, method = method, npc_num = npc_num )
	    	sc_persampledata <- findclusters_persample(sc_persampledata, samplename, method = method, npc_num = npc_num, res_num = res_num)
			DefaultAssay(sc_persampledata) <- "SCT"
		#### pK Identification (no ground-truth) ####---------------------------------------------------------------------------------------
			sweep_res_list_sc_persampledata <- paramSweep_v3(sc_persampledata, PCs = 1:npc_num, sct = TRUE)
			sweep_stats_sc_persampledata <- summarizeSweep(sweep_res_list_sc_persampledata, GT = FALSE)
			bcmvn_sc_persampledata <- find.pK(sweep_stats_sc_persampledata)
			
		#### optimal parameter ####
			mpK<-as.numeric(as.vector(bcmvn_sc_persampledata$pK[which.max(bcmvn_sc_persampledata$BCmetric)]))
			print(paste0(samplename, "_mpK: ", mpK))
			class(mpK)
			#rm(sweep_stats_sc_persampledata, sweep_stats_sc_persampledata, bcmvn_sc_persampledata)

		#### Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
			## (3) Homotypic Doublet Proportion Estimate -------------------------------------
			#annotations <- sc_persampledata$seurat_clusters
			#head(annotations)
			#homotypic.prop <- modelHomotypic(annotations)  
			DoubletRate = ncol(sc_persampledata)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
			print(paste0(samplename, "_doubletrate: ", DoubletRate))
			#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
			nExp_poi <- round(DoubletRate*length(sc_persampledata$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
			print(paste0(samplename, "_nExp_poi: ", nExp_poi))
			# 计算双细胞比例
			#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
			#print(paste0(samplename, "_nExp_poi.adj: ", nExp_poi.adj))
	 		## ex: annotations <- sc_persampledata@meta.data$ClusteringResults
			## Assuming 7.5% doublet formation rate - tailor for your dataset


		#### Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
			sc_persampledata <- doubletFinder_v3(sc_persampledata, PCs = 1:npc_num, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
			label_nExp_poi <- paste0("DF.classifications_0.25", "_", mpK,'_', nExp_poi)
			#sc_persampledata <- doubletFinder_v3(sc_persampledata, PCs = 1:npc_num, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
			#label_nExp_poi_adj <- paste0("DF.classifications_0.25","_", mpK,"_", nExp_poi.adj)

	    #### DoubletFinder results
			#sc_persampledata@meta.data[,"DF_hi_lo_20pc_5res"] <- sc_persampledata$DF.classifications_0.25_0.03_56
			sc_persampledata@meta.data[,"DF_hi_lo_20pc_5res"] <- sc_persampledata@meta.data[,label_nExp_poi]
			sc_persampledata@meta.data$DF_hi_lo_20pc_5res[which(sc_persampledata@meta.data$DF_hi_lo_20pc_5res == "Doublet" )] <- "Doublet_High_Confidience"
			sc_persampledata@meta.data$DF_hi_lo_20pc_5res[which(sc_persampledata@meta.data$DF_hi_lo_20pc_5res == "Singlet" )] <- "Doublet_Low_Zero_Confidience"

			table(sc_persampledata@meta.data$DF_hi_lo_20pc_5res) 
			DoubletConfidience_df <- as.data.frame(table(sc_persampledata@meta.data[ ,"DF_hi_lo_20pc_5res"]))
			DoubletConfidience_df
			Doublet_High_Confidience_cells <- as.character(DoubletConfidience_df[DoubletConfidience_df$Var1=="Doublet_High_Confidience",2]) 
			Doublet_High_Confidience_cells
			Doublet_Low_Zero_Confidience_cells <- as.character(DoubletConfidience_df[DoubletConfidience_df$Var1=="Doublet_Low_Zero_Confidience",2])
			Doublet_Low_Zero_Confidience_cells
			rm(DoubletConfidience_df)
			colnames(sc_persampledata@meta.data)
			DF_hi_lo_20pc_5res_dimplot_subtitle <- labs(subtitle = paste0("TOTAL",as.character(ncol(sc_persampledata)), " /DH", 
													Doublet_High_Confidience_cells, " /S", Doublet_Low_Zero_Confidience_cells))
			mytitle <- labs(title = paste0(samplename,"_",method, "_", as.character(npc_num), "pc_", as.character(res_num), "res_DFprocessing"))  #+ theme(plot.subtitle = element_text(hjust = 0.5))
	       
			total_dimplot <- DimPlot(object = sc_persampledata, reduction = 'umap', label = TRUE, group.by = "seurat_clusters") + 
	                                                mytitle  + theme(plot.subtitle = element_text(hjust = 0.5),  plot.title = element_text(hjust = 0.5))  +
	                                                coord_fixed(ratio=1/1) 
			DF_hi_lo_20pc_5res_dimplot <- DimPlot(sc_persampledata, group.by ="DF_hi_lo_20pc_5res", cols =c("black","red","gold"),reduction = 'umap') +
													mytitle + DF_hi_lo_20pc_5res_dimplot_subtitle + 
													theme(plot.subtitle = element_text(hjust = 0.5),  plot.title = element_text(hjust = 0.5)) +
													coord_fixed(ratio=1/1)


			a<- plot_grid(total_dimplot,DF_hi_lo_20pc_5res_dimplot , ncol = 2, align = "vh") 

		#### Plot results ---------------------------------------------------------------------------
			if(!is.null(dev.list())){dev.off()}
	        pdf(file = paste0(filesDir, "Results/02Clustering/Figures/DFprocessing/", samplename, "_SCT_", as.character(npc_num), "pc_", as.character(res_num), "res_DFprocessing.pdf"), width = 20, height = 12)
	  		print(a)
	  		print(plot(x = as.numeric(as.character(bcmvn_sc_persampledata$pK)), y = bcmvn_sc_persampledata$BCmetric, pch = 16,type="b", col = "blue",lty=1) + 
								        abline(v=mpK,lwd=2,col='red',lty=2) + 
								        title(paste0("The BCmvn distributions mpK:", mpK, " nExp_poi:", nExp_poi)) + 
								        text(mpK,max(bcmvn_sc_persampledata$BCmetric),as.character(mpK),pos = 4,col = "red")
				)
	    	dev.off()
	      
	    return(sc_persampledata)
    })



	names(CN1_23filter_PSCT_20pc.5res_DF)
	names(CN1_23filter_PSCT_20pc.5res_DF) <- namelist
	names(CN1_23filter_PSCT_20pc.5res_DF)
	#saveRDS(CN1_23filter_PSCT_20pc.5res_DF, file = paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_PSCT_20pc.5res_DF_230426.rds"))





#### total_structure  40pc 0.1####
    namelist <- names(CN1_23filter_PSCT_20pc.5res_DF)
    namelist

    CN1_23filter_PSCT_20pc.1res_DF = lapply(namelist,function(samplename){
	    #### get the samplename parameters ####
	    	print(samplename)
	    	sc_persampledata <- CN1_23filter_PSCT_20pc.5res_DF[[samplename]]
	    	print(head(sc_persampledata))
	    	method <- "SCT"
	    	npc_num <- 20
    		res_num <- 0.1

	    #### PCA & findclusters ####
	    	#sc_persampledata <- pca_persample(sc_persampledata, samplename, method = method, npc_num = npc_num )
	    	sc_persampledata <- findclusters_persample(sc_persampledata, samplename, method = method, npc_num = npc_num, res_num = res_num)
			DefaultAssay(sc_persampledata) <- "SCT"
		#### pK Identification (no ground-truth) ####---------------------------------------------------------------------------------------
			sweep_res_list_sc_persampledata <- paramSweep_v3(sc_persampledata, PCs = 1:npc_num, sct = TRUE)
			sweep_stats_sc_persampledata <- summarizeSweep(sweep_res_list_sc_persampledata, GT = FALSE)
			bcmvn_sc_persampledata <- find.pK(sweep_stats_sc_persampledata)
			
		#### optimal parameter ####
			mpK<-as.numeric(as.vector(bcmvn_sc_persampledata$pK[which.max(bcmvn_sc_persampledata$BCmetric)]))
			print(paste0(samplename, "_mpK: ", mpK))
			class(mpK)
			#rm(sweep_stats_sc_persampledata, sweep_stats_sc_persampledata, bcmvn_sc_persampledata)

		#### Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
			## (3) Homotypic Doublet Proportion Estimate -------------------------------------
			#annotations <- sc_persampledata$seurat_clusters
			#head(annotations)
			#homotypic.prop <- modelHomotypic(annotations)  
			DoubletRate = ncol(sc_persampledata)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
			print(paste0(samplename, "_doubletrate: ", DoubletRate))
			#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
			nExp_poi <- round(DoubletRate*length(sc_persampledata$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
			print(paste0(samplename, "_nExp_poi: ", nExp_poi))
			# 计算双细胞比例
			#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
			#print(paste0(samplename, "_nExp_poi.adj: ", nExp_poi.adj))
	 		## ex: annotations <- sc_persampledata@meta.data$ClusteringResults
			## Assuming 7.5% doublet formation rate - tailor for your dataset


		#### Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
			sc_persampledata <- doubletFinder_v3(sc_persampledata, PCs = 1:npc_num, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
			label_nExp_poi <- paste0("DF.classifications_0.25", "_", mpK,'_', nExp_poi)
			#sc_persampledata <- doubletFinder_v3(sc_persampledata, PCs = 1:npc_num, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
			#label_nExp_poi_adj <- paste0("DF.classifications_0.25","_", mpK,"_", nExp_poi.adj)

	    #### DoubletFinder results
			#sc_persampledata@meta.data[,"DF_hi_lo_20pc_5res"] <- sc_persampledata$DF.classifications_0.25_0.03_56
			sc_persampledata@meta.data[,"DF_hi_lo_20pc_0.1res"] <- sc_persampledata@meta.data[,label_nExp_poi]
			sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res[which(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res == "Doublet" )] <- "Doublet_High_Confidience"
			sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res[which(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res == "Singlet" )] <- "Doublet_Low_Zero_Confidience"

			table(sc_persampledata@meta.data$DF_hi_lo_20pc_0.1res) 
			DoubletConfidience_df <- as.data.frame(table(sc_persampledata@meta.data[ ,"DF_hi_lo_20pc_0.1res"]))
			DoubletConfidience_df
			Doublet_High_Confidience_cells <- as.character(DoubletConfidience_df[DoubletConfidience_df$Var1=="Doublet_High_Confidience",2]) 
			Doublet_High_Confidience_cells
			Doublet_Low_Zero_Confidience_cells <- as.character(DoubletConfidience_df[DoubletConfidience_df$Var1=="Doublet_Low_Zero_Confidience",2])
			Doublet_Low_Zero_Confidience_cells
			rm(DoubletConfidience_df)
			colnames(sc_persampledata@meta.data)
			DF_hi_lo_20pc_5res_dimplot_subtitle <- labs(subtitle = paste0("TOTAL",as.character(ncol(sc_persampledata)), " /DH", 
													Doublet_High_Confidience_cells, " /S", Doublet_Low_Zero_Confidience_cells))
			mytitle <- labs(title = paste0(samplename,"_",method, "_", as.character(npc_num), "pc_", as.character(res_num), "res_DFprocessing"))  #+ theme(plot.subtitle = element_text(hjust = 0.5))
	       
			total_dimplot <- DimPlot(object = sc_persampledata, reduction = 'umap', label = TRUE, group.by = "seurat_clusters") + 
	                                                mytitle  + theme(plot.subtitle = element_text(hjust = 0.5),  plot.title = element_text(hjust = 0.5))  +
	                                                coord_fixed(ratio=1/1) 
			DF_hi_lo_20pc_5res_dimplot <- DimPlot(sc_persampledata, group.by ="DF_hi_lo_20pc_0.1res", cols =c("black","red","gold"),reduction = 'umap') +
													mytitle + DF_hi_lo_20pc_5res_dimplot_subtitle + 
													theme(plot.subtitle = element_text(hjust = 0.5),  plot.title = element_text(hjust = 0.5)) +
													coord_fixed(ratio=1/1)


			a<- plot_grid(total_dimplot,DF_hi_lo_20pc_5res_dimplot , ncol = 2, align = "vh") 

		#### Plot results ---------------------------------------------------------------------------
			if(!is.null(dev.list())){dev.off()}
	        pdf(file = paste0(filesDir, "Results/02Clustering/Figures/DFprocessing/", samplename, "_SCT_", as.character(npc_num), "pc_", as.character(res_num), "res_DFprocessing.pdf"), width = 20, height = 12)
	  		print(a)
	  		print(plot(x = as.numeric(as.character(bcmvn_sc_persampledata$pK)), y = bcmvn_sc_persampledata$BCmetric, pch = 16,type="b", col = "blue",lty=1) + 
								        abline(v=mpK,lwd=2,col='red',lty=2) + 
								        title(paste0("The BCmvn distributions mpK:", mpK, " nExp_poi:", nExp_poi)) + 
								        text(mpK,max(bcmvn_sc_persampledata$BCmetric),as.character(mpK),pos = 4,col = "red")
				)
	    	dev.off()
	      
	    return(sc_persampledata)
    })



	names(CN1_23filter_PSCT_20pc.1res_DF)
	names(CN1_23filter_PSCT_20pc.1res_DF) <- namelist
	names(CN1_23filter_PSCT_20pc.1res_DF)

	saveRDS(CN1_23filter_PSCT_20pc.1res_DF, file = paste0(filesDir, "Results/02Clustering/SaveData/CN1_23filter_PSCT_bothDF_230426.rds"))


	rm(list = ls())


