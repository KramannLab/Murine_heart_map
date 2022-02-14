#check and integrate Elivra Forte / Nadina Rosenthal MI heart scRNA
#https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7895/?s_page=2&s_pagesize=25
#all fastq files downloaded and run with cell ranger locally against mm10_gfp
library(Seurat)
library(cowplot)
library(ggplot2)
library(clustree)
library(genesorteR, quietly = TRUE)
library(writexl)
library(dplyr)
library(ggpubr)
library(dittoSeq)
library(reshape2)
library(readr)
library(stringr)
library(pals)
source("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/scripts_for_upload/helper_functions_DEG_GO_ploting.R")

setwd("~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/")
sample="Forte_MI_data"

#load all samples and merge
sc.sham1= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/Sham1_out/outs/filtered_feature_bc_matrix/")
sc.sham2= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/Sham2_out/outs/filtered_feature_bc_matrix/")
sc.none= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/none_out/outs/filtered_feature_bc_matrix/")
sc.MIday1= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday1_out/outs/filtered_feature_bc_matrix/")
sc.MIday3= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday3_out/outs/filtered_feature_bc_matrix/")
sc.MIday5= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday5_out/outs/filtered_feature_bc_matrix/")
sc.MIday7= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday7_out/outs/filtered_feature_bc_matrix/")
sc.MIday14= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday14_out/outs/filtered_feature_bc_matrix/")
sc.MIday28= Read10X(data.dir = "~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/cell ranger outputs/MIday28_out/outs/filtered_feature_bc_matrix/")

sc.list = list(sc.sham1,sc.sham2,sc.MIday1,sc.MIday3,sc.MIday5,sc.MIday7,sc.MIday14,sc.MIday28)
names(sc.list) = c("sc.sham1","sc.sham2","sc.MIday1","sc.MIday3","sc.MIday5","sc.MIday7","sc.MIday14","sc.MIday28")
for (sc.sample in names(sc.list)) {sc.list[[sc.sample]]=CreateSeuratObject(counts = sc.list[[sc.sample]], project = sc.sample, min.cells = 3, min.features = 200)}
sc.none <- CreateSeuratObject(counts = sc.none, project = "sc.none", min.cells = 3, min.features = 200)
samples.merged <- merge(x = sc.none, y = sc.list , project = "Forte_MI_merged")
rm(sc.sham1,sc.sham2,sc.MIday1,sc.MIday3,sc.MIday5,sc.MIday7,sc.MIday14,sc.MIday28)

#filter and process
samples.merged[['percent.mt']] = PercentageFeatureSet(samples.merged, pattern = '^mt-')
samples.merged = subset(samples.merged, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA<15000)
samples.merged = NormalizeData(samples.merged, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
samples.merged = FindVariableFeatures(samples.merged, selection.method = 'vst', nfeatures = 2000,y.cutoff=0.5, verbose = FALSE)
samples.merged <- ScaleData(samples.merged, verbose = FALSE)
samples.merged <- RunPCA(samples.merged, verbose = FALSE)
# t-SNE and Clustering
samples.merged <- RunTSNE(samples.merged, reduction = "pca", dims = 1:24)
samples.merged <- FindNeighbors(samples.merged, reduction = "pca", dims = 1:24)
samples.merged <- FindClusters(samples.merged, resolution = 0.5)
#save whole dataset
saveRDS(samples.merged,file = paste0(sample,"whole_dataset_first_clustering.rds"))

#loaded
samples.merged <- readRDS("~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/Forte_MI_datawhole_dataset_first_clustering.rds")
  
#plot
p1 <- DimPlot(samples.merged, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(samples.merged, reduction = "tsne", group.by = "orig.ident", pt.size = 1) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"_filtered_clusters_res0.5.jpeg"), width=20 , height = 10)


#markers of stromal subsets that were further subclustered in the paper
DotPlot(samples.merged,features = c("Wt1","Dkk3","Mt1","Gsn","Postn","Clu"))

#-------subset to only fibroblast and recluster as in the methods indicated--------------
sample = "Forte_MI_data_fibro.subset"
fibro.subset <- subset(samples.merged,idents = c(0,3,5,6,7,15,17))
fibro.subset = NormalizeData(fibro.subset, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
fibro.subset = FindVariableFeatures(fibro.subset, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
fibro.subset <- ScaleData(fibro.subset, verbose = FALSE)
fibro.subset <- RunPCA(fibro.subset, verbose = FALSE)
# t-SNE and Clustering
fibro.subset <- RunTSNE(fibro.subset, reduction = "pca", dims = 1:20)
fibro.subset <- FindNeighbors(fibro.subset, reduction = "pca", dims = 1:20)
fibro.subset <- FindClusters(fibro.subset, resolution = 0.3)

p1 <- DimPlot(fibro.subset, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(fibro.subset, reduction = "tsne", group.by = "orig.ident", pt.size = 1) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"_filtered_clusters_res0.3.jpeg"), width=20 , height = 10)


# ------subtype markers from the papers-----------
stromal.marker =c("Cxcl14","Dpep1","Sparcl1","Cd248","Pi16","Timp1","Angptl4","Mt2","Acta2","Cthrc1","Stmn1","H2afz","Lyz2","Notch2","Kcnk2","Adamtsl2","Col8a1","Postn","Cilp","Meox1","Thbs4","Wisp2","Sfrp2","Comp","Prg4","Lbr","Wif1","Dkk3","Clu","Wt1","Dmkn","Isg15","Ifit3","Cd52","Ptprc")
DotPlot(fibro.subset,group.by = "celltypes.short", features = stromal.marker,dot.scale = 10,assay = "RNA")+ coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = paste0(sample,"marker from the paper.jpeg"), width=6 , height = 10)

#----------short abroviation of annotation-----------
Idents(fibro.subset)="seurat_clusters"
types <- list("0"="HEpiD","1"="HEpiD","2"="PLS","3"="Myof","4"="MFC","5"="EndD","6"="IR","7"="Epi","8"="IFNr","9"="DC")
#Rename identities
fibro.subset <- RenameIdents(fibro.subset, types)
fibro.subset$celltypes.short <- Idents(fibro.subset)
cell.levels = c("HEpiD","PLS","Myof","MFC","EndD","IR","Epi","IFNr","DC")
fibro.subset$celltypes.short = factor(fibro.subset$celltypes.short,levels = cell.levels)
Idents(fibro.subset)="celltypes.short"
#long annotation names
Idents(fibro.subset)="seurat_clusters"
types <- list("0"="homeostatic epicardial derived","1"="homeostatic epicardial derived","2"="progenitor-like state","3"="myofibroblasts","4"="matrifibrocytes","5"="endocardial-derived","6"="injury response","7"="Epicardium","8"="interferon response","9"="dendritic-like")
#Rename identities
fibro.subset <- RenameIdents(fibro.subset, types)
fibro.subset$celltypes <- Idents(fibro.subset)
cell.levels = c("homeostatic epicardial derived","progenitor-like state","myofibroblasts","matrifibrocytes","endocardial-derived","injury response","Epicardium","interferon response","dendritic-like")
fibro.subset$celltypes = factor(fibro.subset$celltypes,levels = cell.levels)
#short abroviation of annotation
Idents(fibro.subset)="seurat_clusters"

DimPlot(fibro.subset,group.by = "celltypes.short", reduction = "tsne", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_annotation_short.jpeg"), width=10 , height = 10)

DimPlot(fibro.subset,group.by = "celltypes", reduction = "tsne",  pt.size = 0.5)
ggsave(filename = paste0(sample,"_annotation.jpeg"), width=8 , height = 5)


#-----add general sample groups------------
fibro.subset$stim=fibro.subset$orig.ident
MI.control.ls = fibro.subset$stim  %in% c("sc.sham1","sc.sham2","sc.none")
MI.early.ls = fibro.subset$stim  %in% c("sc.MIday1","sc.MIday3","sc.MIday5","sc.MIday7","sc.MIday14","sc.MIday28")
fibro.subset$stim[MI.control.ls]="control"
fibro.subset$stim[MI.early.ls]="MI"
fibro.subset$stim=factor(fibro.subset$stim,levels = c("control","MI"))
#-----order MI time points in orig.ident---------
fibro.subset$orig.ident=factor(fibro.subset$orig.ident,levels=c("sc.none","sc.sham1","sc.sham2","sc.MIday1","sc.MIday3","sc.MIday5","sc.MIday7","sc.MIday14","sc.MIday28"))

#-----save with all annotations-----------------
saveRDS(fibro.subset,file = paste0(sample,"fibroblast_subset_annotated.rds"))
fibro.subset <- readRDS("~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/Forte_MI_datafibroblast_subset_annotated.rds")

#----plot gene expression in fib subtype---------
FeaturePlot(fibro.subset,features = "Wif1",order = T,pt.size = 1.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste(sample,"_Wif1.jpeg"), width=10 , height = 10)

FeaturePlot(fibro.subset,features = "Pdgfrb",order = T,pt.size = 1.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste(sample,"_Pdgfrb.jpeg"), width=10 , height = 10)

#------subset to only control samples and recluster
Idents(fibro.subset)="stim"
fibro.subset.control = subset(fibro.subset,ident="control")
fibro.subset.control<-recluster_RNA(fibro.subset.control,0.5)
dp1<-DimPlot(fibro.subset.control,group.by = "celltypes.short",  label = TRUE, pt.size = 1, label.size = 8,repel = T) + NoLegend() + ggtitle("Control samples Forte et al")
f1<-FeaturePlot(fibro.subset.control,features = "Postn",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggarrange(dp1,f1)
DimPlot(fibro.subset.control,group.by = "celltypes.short",  label = TRUE, pt.size = 2, label.size = 12,repel = T) + NoLegend() + ggtitle("Control samples Forte et al")
ggsave(filename = paste0(sample,"_onyl_control.jpeg"), width=10 , height = 10) 

FeaturePlot(fibro.subset.control,features = "Postn",order = T,pt.size = 2)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste0(sample,"_onyl_control_postn2.jpeg"), width=10 , height = 10) 




