#### prepare core&environment ####
rm(list=ls())
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
library(ggalluvial)
library(RColorBrewer)

options(stringsAsFactors = F)
filesDir<-'I:/IPMN_PROGRAM_LINUX/analysis/sc_analysis/sc_analysis6_CN11UMCN2CN3/'
setwd(filesDir)
getwd()

method = "RNA"
npc_num <- 20
set_resolutions = seq(0.1, 1.5, by = 0.1)
method = "RNA"

reduction <- "umap"

identlevel <- "resfind_3"
samplename <- "CD4T"
method = "RNA"


filtercelldata <- readRDS(paste0("CD4T_20samples.rds"))
colnames(filtercelldata@meta.data)
table(filtercelldata$SubCellType1_CN11UMCN2)


####vlnplot ####
### datasave ####         
reduction <- "umap"
identlevel <- "DEGs_score_240527"
samplename <- "CD4T"
method = "RNA"


FigDir <- paste0(filesDir, "06TCellAnalysis/Figures/", samplename,"/",identlevel,"/")
TableDir <- paste0(filesDir, "06TCellAnalysis/Tables/", samplename,"/",identlevel,"/")
DataDir <- paste0(filesDir, "06TCellAnalysis/SaveData/", samplename,"/",identlevel,"/")


if(!dir.exists(paste0(FigDir))){dir.create(paste0(FigDir), recursive = TRUE)}
if(!dir.exists(paste0(TableDir))){dir.create(paste0(TableDir), recursive = TRUE)}
if(!dir.exists(paste0(DataDir))){dir.create(paste0(DataDir), recursive = TRUE)}

mycolor = c("CD4T_CCR7_SELL" = "#99749f", #深紫色
            "CD4T_CCR7_1" = "#d1aaeb", #浅紫色
            "CD4T_CCR7_2" = "#9c27b0", #紫粉色
            "CD4T_ANXA1_TRADD" = "#c2f4cf", #浅绿色
            "CD4T_CCL4_IFNG" = "#9cdfb3", #中绿色
            "CD4T_CCL4_NKG7" = "#3c9a68", #深绿色
            "CD4T_TNFSF13B_CCL20" = "#95b7d0", #深蓝色
            "CD4T_FOXP3_TNFRSF9" = "#eb8c30", #深黄色
            "CD4T_CXCL13_PDCD1" = "#efbb4b" #浅黄色 
)

#### score ####
 scorelist1 =list(
 exhaustion_score_gut = c("HAVCR2", "LAG3","ENTPD1"),
  cytotoxic_score_gut = c("PRF1","IFNG","GNLY","NKG7","GZMB","GZMA","GZMH","KLRK1","KLRB1","KLRD1","CTSW","CST7"),
  inflammatory_score_gut = c("IRF1","CD8A","CCL2","CCL3","CCL4","CXCL9","CXCL10","ICOS","GZMK"),
 effector_score_gut = c("GZMA","GZMB","PRF1","EOMES","IFNG","TNF","CXCL9","CXCL10","CD8A","CD4","FOXP3","ICOS","CTLA4"),
 naive_score_science = c("CCR7", "TCF7", "LEF1", "SELL", "MAL"),
exhuastion_TF_science = c("SOX4", "FOXP3", "TOX", "TOX2", "RBPJ", " ZBED2", "PRDM1", "IKZF4", "BATF", "STAT3", "IFI16")
)  


# "naive_score_science")] <- "naive"
# "cytotoxic_score_gut")] <- "cytotoxic"
# "effector_score_gut")] <- "effector"
# "exhaustion_score_gut")] <- "exhuastion"
# "exhuastion_TF_science")] <- "exhuastion_TF"


table(filtercelldata$SubCellTypeup_CN11UMCN2)
table(filtercelldata$SubCellType1_CN11UMCN2)
ncol(filtercelldata)
Idents(filtercelldata) <- filtercelldata$SubCellType1_CN11UMCN2
themefamily <- 
  theme(axis.text.x = element_text( color = "#000000"),
        axis.text.y = element_text( color = "#000000"))


scoregenenamelist <- names(scorelist)
genelist_filt0 = lapply(scoregenenamelist,function(genesetname){
  geneset <- scorelist[[genesetname]]
  geneset <- geneset[geneset%in%rownames(filtercelldata)]
  geneset <- unique(geneset)
  return(geneset)
})
scoregenenamelist
names(genelist_filt0) <- scoregenenamelist
names(genelist_filt0)

scoregenenamelist <- names(genelist_filt0)
genelistkegg_filt = lapply(scoregenenamelist,function(genesetname){
  geneset <- genelist_filt0[[genesetname]]
  geneexpr <-GetAssayData(object = filtercelldata[["RNA"]], slot = "data")[geneset,]
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
  filtercelldata <- AddModuleScore(
    filtercelldata,
    featureset,
    name = genesetname
  )
}

#### score vlnplot SubCellType1_CN11UMCN2 ####

medata  <-  filtercelldata@meta.data 
colnames(medata)[which(colnames(medata) == "naive_score_science1")] <- "naive_score"
colnames(medata)[which(colnames(medata) == "cytotoxic_score_gut1")] <- "cytotoxic_score"
colnames(medata)[which(colnames(medata) == "effector_score_gut1")] <- "effector_score"
colnames(medata)[which(colnames(medata) == "exhaustion_score_gut1")] <- "exhaustion_score"
colnames(medata)[which(colnames(medata) == "exhuastion_TF_science1")] <- "exhaustion_TF_score"
scorelevel <- c("naive_score","cytotoxic_score","effector_score","exhaustion_score","exhaustion_TF_score")
medata$SubCellType1_CN11UMCN2 <- factor(medata$SubCellType1_CN11UMCN2, levels = names(mycolor))
colnames(medata)

for (scorename in scorelevel) {
  
  p4 <- ggplot(data = medata, aes_string(x = "SubCellType1_CN11UMCN2", y = scorename)) +
    geom_point(position = 'jitter', color = 'grey',
               size = 0.5, alpha = 0.8) +
    geom_violin(aes(fill = SubCellType1_CN11UMCN2),
                color = 'grey', alpha = 1,
                scale = 'width',
                linewidth = 0.5, #外轮廓粗细
                trim = TRUE) +
    geom_boxplot(color = 'white',
                 outlier.color = 'black',
                 width = 0.4, #箱子宽度
                 size = 0.5, #外轮廓描边粗细
                 fill = NA) +
    theme_bw()+
    themefamily + 
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    ggtitle(paste(scorename," score")) +
    scale_fill_manual(values = mycolor) +
    scale_color_manual(values = mycolor)# + #+RotatedAxis() +

  p4
  ggsave(paste0(FigDir,"score_240527/",scorename,"_vlnplot2.pdf"),width = 4, height = 2.5)
  ggsave(paste0(FigDir,"score_240527/",scorename,"_vlnplot2.tiff"),width = 4, height = 2.5)
  
  
}
FigDir

####宽数据转成长数据 #####
scorelevel <- c("naive_score","cytotoxic_score","exhaustion_score","exhaustion_TF_score")
colnames(medata)
view(medata)
medata$cellbarcode <- rownames(medata)
medatalong <- medata %>% 
  pivot_longer(cols = scorelevel,
               names_to = "CD4T_score_type",
               values_to = "score",
               values_drop_na = T
  )
view(medatalong)
colnames(medatalong)

saveRDS(medatalong, paste0(DataDir,"CD4T_score_sourcedata.rds"))

medatalong <- readRDS(paste0(DataDir,"CD4T_score_sourcedata.rds"))
#### source data ####
colnames(medatalong)
#view(df)
sourcedata <- medatalong
colnames(sourcedata)
head(sourcedata)
sourcedata <-  sourcedata %>%  select(cellbarcode, orig.ident, SampleStage2,CD4T_score_type,score )
saveRDS(sourcedata, paste0(DataDir,"v1_fig1e_s2b_d_sourcedata.rds")) #### figs2b_d

#### score radater plot SampleStage2 ####
themefamily <- theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1,  color = "#000000",size = 14),
                     #axis.text.x = element_text(angle = 90, hjust = 1, color = "#000000"),
                     axis.text.y =  element_text(angle = 0, hjust = 0.5, color = "#000000",size = 14),
                     plot.title = element_text(hjust = 0.5,color = "#000000"),
                     axis.title.y =element_text( color="#000000",hjust=0.5,size = 14),
                     axis.title.x =element_text( color="#000000",hjust=0.5,size = 14),
                     panel.border = element_rect(color = "#000000", size = 0.5, linetype="solid", fill = NA),
                     legend.position = "none",
                     strip.background = element_rect(
                       fill="#efebbe", size= NULL, linetype= NULL, color =  NULL),
                     plot.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = 'pt'),
                     strip.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5),
                                                 size = 15,  hjust = 0.5, color = "#000000")
)


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
table(medatalong$SampleStage2)
medatalong$SampleStage2 <- factor(medatalong$SampleStage2,levels = c("UNIN","IPMN","PDAC"))
colnames(medatalong)

plist <- list()
for (scorelevelname in scorelevel) {
  scorename <- "score"
  
  #scorelevelname <- "exhaustion_score"
  medatalong1 <- medatalong %>% filter(CD4T_score_type %in% scorelevelname)
  
  maxnum <-  max(medatalong1[,scorename])
  y_positionlist <- positionfunc(maxnum)
  
  table(medatalong1$SampleStage2)
  p6 <- ggplot(data = medatalong1, aes_string(x = "SampleStage2", y = scorename)) +
    facet_wrap(vars(CD4T_score_type), ncol = 1, scales = "free_y") + 
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
    ) # 检验的类型
  
  p6
  
  plist[[scorelevelname]] <- p6
  ggsave(paste0(FigDir,"score_",scorelevelname,"_vlnplot_Samplestage2_sig_pvalue.pdf"),width = 2, height = 4)
  ggsave(paste0(FigDir,"score_",scorelevelname,"_vlnplot_Samplestage2_sig_pvalue.tiff"),width = 2, height = 4)
}

table(medatalong$SampleStage2)
plist$naive_score
plist$cytotoxic_score
