#check and integrate McLellan data AngII heart scRNA
#https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8810/files/
#https://pubmed.ncbi.nlm.nih.gov/32795101/
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

setwd("~/Documents/FP_scRNA/R stuff/McLellan paper data angII heart scRNA/")
sample="McLellan_data_"

#colorblind friendly color panels for all figures
#for conditions/stim #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
sham.col = Okabe_Ito[1]
tac.col = Okabe_Ito[2]
tac28.col = Okabe_Ito[3]
stim.col=c(sham.col,tac.col,tac28.col)
stim.col.light=c("#E69F0080","#56B4E980","#009E7380")
stim.col.Ang=c("#F0E442", "#0072B2", "#D55E00",sham.col,tac.col,tac28.col)
#for genotypes #From Paul Tol: https://personal.sron.nl/~pault/ '#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'
col.pdgfrb="#EE6677"
col.col1a1="#228833"
col.gli1="#4477AA"
geno.cl = c(col.col1a1,col.pdgfrb,col.gli1)
geno.cl.light = c("#22883380","#EE667780","#4477AA80")
#for subclustering #RColorBrewer
display.brewer.pal(name = "Paired",n = 12)
cluster.col=brewer.pal(name = "Paired",n=11)[1:6]
cluster.col.Ang = c(brewer.pal(name = "Paired",n=11)[1:6],colorBlindness::PairedColor12Steps)
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)

#load preprocessed data downloadad from EBI
data <- read_tsv("~/Documents/FP_scRNA/R stuff/McLellan paper data angII heart scRNA/processed_full_count_matrix.tsv") # processed data from the paper
data <- data.frame(data)
head(data)
#one cell was removed by not optmial data loading
rownames(data) <- data$AAACCTGAGGTGCTAG_2
data$AAACCTGAGGTGCTAG_2 <- NULL

scMacL <- CreateSeuratObject(data,project = sample, min.cells = 3, min.features = 200)

scMacL[['percent.mt']] = PercentageFeatureSet(scMacL, pattern = '^mt-')
VlnPlot(scMacL, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
ggsave(filename = paste0(sample,"_Gene&RNA count.jpeg"), width=10 , height = 10)
#qualtity filter
scMacL = subset(scMacL, subset = nFeature_RNA > 100 & nFeature_RNA < 15000 & nFeature_RNA < 50000 & percent.mt < 30) # in the original paper already done

#process as indicated in the methods of the paper
scMacL = NormalizeData(scMacL, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
scMacL = FindVariableFeatures(scMacL, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
scMacL = ScaleData(scMacL, verbose = FALSE, features = rownames(scMacL))
scMacL = RunPCA(scMacL, verbose = FALSE)
scMacL = FindNeighbors(scMacL, reduction = 'pca', dims = 1:30, verbose = FALSE)
scMacL = FindClusters(scMacL, resolution = 1.2, verbose = FALSE)
scMacL = RunTSNE(scMacL, reduction = 'pca', dims = 1:30, verbose = FALSE)

#plot
p1 <- DimPlot(scMacL, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(scMacL, reduction = "tsne", group.by = "orig.ident", pt.size = 1) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"_filtered_clusters_res1.2.jpeg"), width=20 , height = 10)

#add condition and sample to the cells
sample_id=c(scMacL$orig.ident)
stim=c(scMacL$orig.ident)
i=1
for (cell in colnames(scMacL)) {
  if (str_detect(string = cell, pattern = "2")) {
    sample_id[i] = "none1"
    stim[i]="none"}
  if (str_detect(string = cell, pattern = "3")) {
    sample_id[i] = "angII_1"
    stim[i]="angII"}
  if (str_detect(string = cell, pattern = "4")) {
    sample_id[i] = "saline1"
    stim[i]="saline"}
  if (str_detect(string = cell, pattern = "5")) {
    sample_id[i] = "angII_2"
    stim[i]="angII"}
  if (str_detect(string = cell, pattern = "6")) {
    sample_id[i] = "none2"
    stim[i]="none"}
  if (str_detect(string = cell, pattern = "7")) {
    sample_id[i] = "angII_3"
    stim[i]="angII"}
  if (str_detect(string = cell, pattern = "8")) {
    sample_id[i] = "saline2"
    stim[i]="saline"}
  if (str_detect(string = cell, pattern = "9")) {
    sample_id[i] = "angII_4"
    stim[i]="angII"}
  i=i+1
}
sc



$stim = factor(stim,levels = c("none","saline","angII"))
scMacL$orig.ident = factor(sample_id,levels = c("none1","none2","saline1","saline2","angII_1","angII_2","angII_3","angII_4"))
#save
saveRDS(scMacL,file = paste0(sample,"seurat_object.RDS"))
#load
scMacL<-readRDS("~/Documents/FP_scRNA/R stuff/MacLellan paper data angII heart scRNA/MacLellan_data_seurat_object.RDS")

##########------------------------------subcluster fibroblasts from mclellan and get cluster markers###############################################
scMacL.fib =subset(scMacL,idents = c("0","2","3","4","5","8","9","23","25"))
sample="McLellan_fibro_sub_"

#recluster to generate 9 subclusters
scMacL.fib = FindNeighbors(scMacL.fib, reduction = 'pca', dims = 1:30, verbose = FALSE)
scMacL.fib = FindClusters(scMacL.fib, resolution = 0.6, verbose = FALSE)

#plot
p1 <- DimPlot(scMacL.fib, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(scMacL.fib, reduction = "tsne", group.by = "orig.ident", pt.size = 1) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"_clusters_res0.6.jpeg"), width=20 , height = 10)

Idents(scMacL.fib)="seurat_clusters"
types <- list("0"="Fibroblast 1","1"="Fibroblast 2","2"="Fibroblast 3","3"="Fibroblast 4","4"="Fibroblast 5","5"="Fibroblast-Cilp","6"="Fibroblast 6","7"="Fibroblast-Wif1","8"="Fibroblast-Thbs4","9"="Fibroblast 7","10"="Fibroblast 8")
#Rename identities
scMacL.fib <- RenameIdents(scMacL.fib, types)
scMacL.fib$celltypes <- Idents(scMacL.fib)
cell.levels = c("Fibroblast 1","Fibroblast 2","Fibroblast 3","Fibroblast 4","Fibroblast 5","Fibroblast-Cilp","Fibroblast 6","Fibroblast-Wif1","Fibroblast-Thbs4","Fibroblast 7","Fibroblast 8")
scMacL.fib$celltypes = factor(scMacL.fib$celltypes,levels = cell.levels)
#short abroviation of annotation
Idents(scMacL.fib)="seurat_clusters"
types <- list("0"="Fib1","1"="Fib2","2"="Fib3","3"="Fib4","4"="Fib5","5"="FibCilp","6"="Fib6","7"="FibWif1","8"="FibThbs4","9"="Fib7","10"="Fib8")
#Rename identities
scMacL.fib <- RenameIdents(scMacL.fib, types)
scMacL.fib$celltypes.short <- Idents(scMacL.fib)
cell.levels = c("Fib1","Fib2","Fib3","Fib4","Fib5","FibCilp","Fib6","FibWif1","FibThbs4","Fib7","Fib8")
scMacL.fib$celltypes.short = factor(scMacL.fib$celltypes.short,levels = cell.levels)
Idents(scMacL.fib)="seurat_clusters"
#add celltypes + stim with levels
lvl=NULL
for (levels in levels(scMacL.fib$celltypes)) {lvl=c(lvl,paste0(levels,"_",levels(scMacL.fib$stim)))}
scMacL.fib$celltypes.stim=factor(paste0(scMacL.fib$celltypes,"_",scMacL.fib$stim),levels = lvl)

#-----save with all annotations-----------------
saveRDS(scMacL.fib,file = paste0(sample,"fibroblast_subset_annotated.rds"))
#load
scMacL.fib<-readRDS("~/Documents/FP_scRNA/R stuff/MacLellan paper data angII heart scRNA/McLellan_fibro_sub_fibroblast_subset_annotated.rds")

#-------plot annotation---------
DimPlot(scMacL.fib,group.by = "celltypes.short", reduction = "tsne", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_annotation_short.jpeg"), width=10 , height = 10)

DimPlot(scMacL.fib,group.by = "celltypes", reduction = "tsne",  pt.size = 0.5)
ggsave(filename = paste0(sample,"fibro.subset_annotation.jpeg"), width=8 , height = 5)

#------subset to only control samples and recluster
Idents(scMacL.fib)="stim"
scMacL.fib.control = subset(scMacL.fib,ident=c("none","saline"))
scMacL.fib.control<-recluster_RNA(scMacL.fib.control,0.5)
dp1<-DimPlot(scMacL.fib.control,group.by = "celltypes.short",  label = TRUE, pt.size = 1, label.size = 8,repel = T) + NoLegend() + ggtitle("Control samples mclellan et al")
f1<-FeaturePlot(scMacL.fib.control,features = "Postn",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggarrange(dp1,f1)
ggsave(filename = paste0(sample,"_onyl_control_postn.jpeg"), width=10 , height = 5) 
DimPlot(scMacL.fib.control,group.by = "celltypes.short",  label = TRUE, pt.size = 2, label.size = 12,repel = T) + NoLegend() + ggtitle("Control samples mclellan et al")
ggsave(filename = paste0(sample,"_onyl_control.jpeg"), width=10 , height = 10) 

FeaturePlot(scMacL.fib.control,features = "Postn",order = T,pt.size = 2)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste0(sample,"_onyl_control_postn2.jpeg"), width=10 , height = 10) 

