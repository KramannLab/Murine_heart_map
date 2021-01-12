#http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html
#script for the integration and analysis of all Sham and TAC samples based on only the basic filtered RDS files
library(harmony)
library(Seurat)
library(cowplot)
library(ggplot2)
library(clustree)
library(genesorteR, quietly = TRUE)
library(writexl)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(dittoSeq)
library(reshape)
setwd("~/Documents/FP_scRNA/R stuff/PVM harmony/")
set.seed(111)
sample="harmony_full_int"

#colors for genotypes
col.pdgfrb="#EC407A"
col.myh11="#FFEB3B"
col.cdh5 = "#9C27B0"
col.ng2="#795548"
col.col1a1="#2196F3"
col.gli1="#4CAF50"
genotype.colors = c(col.pdgfrb,col.cdh5,col.col1a1,col.myh11,col.gli1,col.ng2)

#our handpicked color code for clusters
color_code=c("#EC407A","#9C27B0" ,"#2196F3" ,"#4CAF50" ,"#FFEB3B" ,"#FF5722" ,"#009688" ,"#607D8B" ,"#303F9F" ,"#795548" ,"#FFB300" ,"#DCE775")
  
hellblau = "#2E8DCC"
organe ="#D95F02"
lila ="#7570B3"
violet ="#E7298A" 
gelb ="#E6AB02"
UKAgrün ="#669900"
grau = "#607D8B"

cluster.colors = c(hellblau, organe, lila, UKAgrün,violet ,gelb ,grau)

######################################-------load all prefiltered datasets (also log normalized)------#############################################################
datapath="~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering/"
#load for ------PDGFRb------ samples
PDGFRb_Sham = 'FP3_PDGFRb_Sham'
PDGFRb_TAC = 'FP14_PDGFRb_TAC'
PDGFRb_Sham=readRDS(file = paste0(datapath, PDGFRb_Sham , "_seurat_object_after_bacis_filters.rds"))
PDGFRb_TAC=readRDS(file = paste0(datapath, PDGFRb_TAC , "_seurat_object_after_bacis_filters.rds"))

#load for ------Gli1------ samples
Gli1_Sham = 'FP13_Gli1_Sham'
Gli1_TAC = 'FP12_Gli1_TAC'
Gli1_Sham=readRDS(file = paste0(datapath, Gli1_Sham , "_seurat_object_after_bacis_filters.rds"))
Gli1_TAC=readRDS(file = paste0(datapath, Gli1_TAC , "_seurat_object_after_bacis_filters.rds"))

#load for ------Col1a1------ samples
Col1a1_Sham = 'FP19_Col1a1_Sham'
Col1a1_TAC = 'FP18_Col1a1_TAC'
Col1a1_Sham=readRDS(file = paste0(datapath, Col1a1_Sham , "_seurat_object_after_bacis_filters.rds"))
Col1a1_TAC=readRDS(file = paste0(datapath, Col1a1_TAC , "_seurat_object_after_bacis_filters.rds"))

#load for ------Cdh5------ samples
Cdh5_Sham = 'FP21_Cdh5_Sham'
Cdh5_TAC = 'FP20_Cdh5_TAC'
Cdh5_Sham=readRDS(file = paste0(datapath, Cdh5_Sham , "_seurat_object_after_bacis_filters.rds"))
Cdh5_TAC=readRDS(file = paste0(datapath, Cdh5_TAC , "_seurat_object_after_bacis_filters.rds"))

#load for ------NG2------ samples
NG2_Sham = 'FP23_NG2_Sham'
NG2_TAC = 'FP22_NG2_TAC'
NG2_Sham=readRDS(file = paste0(datapath, NG2_Sham , "_seurat_object_after_bacis_filters.rds"))
NG2_TAC=readRDS(file = paste0(datapath, NG2_TAC , "_seurat_object_after_bacis_filters.rds"))

#load for ------Myh11------ samples
Myh11_Sham = 'FP56_Myh11_Sham'
Myh11_TAC = 'FP57_Myh11_TAC'
Myh11_Sham=readRDS(file = paste0(datapath, Myh11_Sham , "_seurat_object_after_bacis_filters.rds"))
Myh11_TAC=readRDS(file = paste0(datapath, Myh11_TAC , "_seurat_object_after_bacis_filters.rds"))

###########----merge for harmony

#prepare merged Seurat object for harmony
samples.merged <- merge(x = PDGFRb_Sham, y =  list(PDGFRb_TAC, Cdh5_TAC, Gli1_TAC, Col1a1_TAC, NG2_TAC, Myh11_TAC, Cdh5_Sham, Gli1_Sham, Col1a1_Sham, NG2_Sham, Myh11_Sham), project = "PVM")
sample = "full_integration_harmony"
#first processing
samples.merged <- samples.merged %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = samples.merged@var.genes, npcs = 20, verbose = FALSE)

###---subset the individual stim with just TAC or Sham--###
Samples.filtered=samples.merged
TACs = Samples.filtered@meta.data$stim %in% c("FP12_Gli1_TAC","FP22_NG2_TAC","FP14_PDGFRb_TAC","FP18_Col1a1_TAC","FP20_Cdh5_TAC","FP57_Myh11_TAC")
Samples.filtered@meta.data$stim[!TACs]="Sham"
Samples.filtered@meta.data$stim[TACs]="TAC" 
samples.merged=Samples.filtered
remove(Samples.filtered)

# add column for only the genotype
Samples.filtered=samples.merged
Gli1.ls = Samples.filtered@meta.data$orig.ident %in% c("FP12_Gli1_TAC","FP13_Gli1_Sham")
PDGFRb.ls = Samples.filtered@meta.data$orig.ident %in% c("FP14_PDGFRb_TAC","FP3_PDGFRb_Sham")
Cdh5.ls = Samples.filtered@meta.data$orig.ident %in% c("FP21_Cdh5_Sham","FP20_Cdh5_TAC")
Col1a1.ls = Samples.filtered@meta.data$orig.ident %in% c("FP18_Col1a1_TAC","FP19_Col1a1_Sham")
NG2.ls = Samples.filtered@meta.data$orig.ident %in% c("FP22_NG2_TAC","FP23_NG2_Sham")
Myh11.ls = Samples.filtered@meta.data$orig.ident %in% c("FP56_Myh11_Sham","FP57_Myh11_TAC")
Samples.filtered@meta.data$orig.geno[Gli1.ls]="Gli1"
Samples.filtered@meta.data$orig.geno[PDGFRb.ls]="PDGFRb"
Samples.filtered@meta.data$orig.geno[Cdh5.ls]="Cdh5"
Samples.filtered@meta.data$orig.geno[Col1a1.ls]="Col1a1"
Samples.filtered@meta.data$orig.geno[NG2.ls]="NG2"
Samples.filtered@meta.data$orig.geno[Myh11.ls]="Myh11"
samples.merged=Samples.filtered
remove(Samples.filtered)

#settle the order
samples.merged$orig.ident = factor(samples.merged$orig.ident,levels = c("FP19_Col1a1_Sham","FP18_Col1a1_TAC","FP3_PDGFRb_Sham","FP14_PDGFRb_TAC","FP13_Gli1_Sham","FP12_Gli1_TAC","FP21_Cdh5_Sham","FP20_Cdh5_TAC","FP56_Myh11_Sham","FP57_Myh11_TAC","FP22_NG2_TAC","FP23_NG2_Sham"))
samples.merged$orig.geno = factor(samples.merged$orig.geno, levels = c("Col1a1","PDGFRb","Gli1","Cdh5","Myh11","NG2"))

#see batch effect before harmony
samples.merged.raw <- samples.merged %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
samples.merged.raw$geno.stim = paste0(samples.merged.raw$orig.geno,"_",samples.merged.raw$stim)
DimPlot(object = samples.merged.raw, reduction = "umap", pt.size = .1, group.by = "geno.stim")+theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"_raw_data_just_merged.jpeg"),width = 10,height = 10)
DimPlot(object = samples.merged.raw, reduction = "umap", pt.size = .1, split.by = "geno.stim",ncol=3)
ggsave(filename = paste0(sample,"_raw_data_just_merged_split.jpeg"),width = 15,height = 20)
DimPlot(object = samples.merged.raw, reduction = "umap", pt.size = .1, split.by = "orig.geno",ncol=3)
ggsave(filename = paste0(sample,"_raw_data_just_merged_split_geno.jpeg"))
saveRDS(samples.merged.raw,file = paste0(sample,"_merged_raw.rds"))
samples.merged.raw <- readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM harmony/analysis PVM_full integration v10 harmony/RDS/harmony_full_int_merged_raw.rds")

#run harmony on the data
#start processing
samples.merged <- samples.merged %>% RunHarmony("orig.ident", plot_convergence = TRUE,epsilon.cluster = -Inf)
saveRDS(samples.merged,file = paste0(sample,"_merged_harmony.rds"))
samples.merged <- readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM harmony/analysis PVM_full integration v10 harmony/RDS/harmony_full_int_merged_harmony.rds")
#process based on harmony result
samples.merged <- samples.merged %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
#see batcheffect eliminated
samples.merged$geno.stim = paste0(samples.merged$orig.geno,"_",samples.merged$stim)
DimPlot(object = samples.merged, reduction = "umap", pt.size = .1, group.by = "geno.stim")+theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"_after_harmony.jpeg"),width = 10,height = 10)
DimPlot(object = samples.merged, reduction = "umap", pt.size = .1, split.by = "geno.stim",ncol=3)
ggsave(filename = paste0(sample,"_after_harmony_split.jpeg"),width = 15,height = 20)
DimPlot(object = samples.merged, reduction = "umap", pt.size = .1, split.by = "orig.geno",ncol=3)
ggsave(filename = paste0(sample,"_after_harmony_split_geno.jpeg"))
#-------------ploting only tdtom
#quality check on tdTom expression before filtering
FeaturePlot(samples.merged, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(sample,"_tdTom_before.jpeg"), width=10 , height =10)
VlnPlot(samples.merged, features = c("tdTomato-Green"), pt.size = 0)+NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_tdTom_Vln_before.jpeg"), width=10 , height =10,limitsize = F)
RidgePlot(samples.merged,features = "tdTomato-Green",group.by = "orig.ident",slot = "counts",assay = "RNA")
ggsave(filename = paste0(sample,"_tdTom_reads_before.jpeg"), width=10 , height = 10)

#------plot to explain filtering tdTom-----------------------------
Idents(samples.merged) = "orig.ident"
myh11.subset = subset(samples.merged,idents = "FP56_Myh11_Sham")
Idents(myh11.subset)="seurat_clusters"
DimPlot(myh11.subset, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend() + xlim(-10,12) + ylim(-14,10)
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove.jpeg", width=10 , height = 10)
FeaturePlot(myh11.subset, features = c("tdTomato-Green"),pt.size = 2) + xlim(-10,12) + ylim(-14,10)
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove_tdTom.jpeg", width=10 , height = 10)
VlnPlot(myh11.subset, features = c("tdTomato-Green"), pt.size = 0)+NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove_tdTom_vln.jpeg", width=10 , height = 10)
table(myh11.subset$seurat_clusters) #remove cluster with less than 20 cells 
myh11.subset = subset(myh11.subset,idents = c("0","1","2","3","6","18","21","22","23"),invert=T)
myh11.subset$seurat_clusters=Idents(myh11.subset)
DimPlot(myh11.subset, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()+ xlim(-10,12) + ylim(-14,10)
ggsave(filename = "FP56_Myh11_Sham_before_tdTom-remove_lower20removed.jpeg", width=10 , height = 10)
sg = sortGenes(myh11.subset@assays$RNA@data, Idents(myh11.subset))
myh11.subset$tdTom.score=1
i=1
for (cluster in levels(factor(myh11.subset$seurat_clusters))) {myh11.subset$tdTom.score[myh11.subset$seurat_clusters==cluster]=round((sg$condGeneProb)["tdTomato-Green",i],digits = 3);i=i+1}
Idents(myh11.subset)="tdTom.score"
DimPlot(myh11.subset, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12,label.box = T,repel = T) + NoLegend() + xlim(-10,12) + ylim(-14,10)
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove_td-score.jpeg", width=10 , height = 10)
Idents(myh11.subset)="seurat_clusters"
myh11.subset = subset(myh11.subset,idents = c("4","5","7","9","10","11","13","14","16","17"),invert=T)
FeaturePlot(myh11.subset, features = c("tdTomato-Green"),pt.size = 2) + xlim(-10,12) + ylim(-14,10)
ggsave(filename = "Only_FP56_Myh11_Sham_tdTom-removed.jpeg", width=10 , height = 10)


#---run specialized function to filter out lowtdtomato cluster contributions per sample, clusters with less than 20 cells and 0 tdtom reads cells
source("~/Documents/FP_scRNA/R stuff/scripts/PVM scripts of current relevant analysis/filterLowTdtomCl_function_harmony.R")
samples.merged.td.filter<- filterLowTdtomCl.harmony(samples.merged,sample)
saveRDS(samples.merged.td.filter,file = paste0(sample,"_merged_td_filter.rds"))
samples.merged.td.filter <- readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM harmony/analysis PVM_full integration v10 harmony/RDS/harmony_full_int_merged_td_filter.rds")

#-------------ploting only tdtom
#quality check on tdTom expression after filtering
FeaturePlot(samples.merged.td.filter, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(sample,"_tdTom_after.jpeg"), width=10 , height =10)
VlnPlot(samples.merged.td.filter, features = c("tdTomato-Green"), pt.size = 0) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_tdTom_Vln_after.jpeg"), width=10 , height =10,limitsize = F)
RidgePlot(samples.merged.td.filter,features = "tdTomato-Green",group.by = "orig.ident",slot = "counts",assay = "RNA")
ggsave(filename = paste0(sample,"_tdTom_reads_after.jpeg"), width=10 , height = 10)
#clusters after removing <20cell cluster contributions and 0-low tdtom cells/cluster contributions
samples.merged.td.filter$geno.stim = paste0(samples.merged.td.filter$orig.geno,"_",samples.merged.td.filter$stim)
DimPlot(object = samples.merged.td.filter, reduction = "umap", pt.size = .1, group.by = "geno.stim") +theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"_after_harmony_td.jpeg"),width = 10,height = 10)
DimPlot(object = samples.merged.td.filter, reduction = "umap", pt.size = .1, split.by = "geno.stim",ncol=3)
ggsave(filename = paste0(sample,"_after_harmony_td_split.jpeg"),width = 15,height = 20)
DimPlot(object = samples.merged.td.filter, reduction = "umap", pt.size = .1, split.by = "orig.geno",ncol=3)
ggsave(filename = paste0(sample,"_after_harmony_td_split_geno.jpeg"))

#------------basic nfeature and ncount plots------------------
#new clusters after tdtom filter
samples.merged.td.filter <- samples.merged.td.filter %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
#before
FeaturePlot(samples.merged.td.filter, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount_before.jpeg"), width=10 , height = 5)
VlnPlot(samples.merged.td.filter,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount_before_vln.jpeg"), width=10 , height =5)
VlnPlot(samples.merged.td.filter,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0.1) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_before_vln.jpeg"), width=10 , height =10)


#----------subsetting to filter out feature count low or high clusters (low quality cells and potential doublets)-----------
samples.filtered=subset(samples.merged.td.filter,idents = c("12","19"),invert=T) 
#after
FeaturePlot(samples.filtered, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount_after.jpeg"), width=10 , height = 5)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount_after_vln.jpeg"), width=10 , height =5)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0.1) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_after_vln.jpeg"), width=10 , height =10)

#new clusters after nfeature remove
samples.filtered <- samples.filtered %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  identity()

#set resolution for analysis
samples.filtered <- samples.filtered %>% FindClusters(resolution = 0.1)
DimPlot(object = samples.filtered, reduction = "umap", pt.size = .1,label=T,label.size = 10)+NoLegend()  
ggsave(filename = paste0(sample,"_clusters_res0.1.jpeg"), width=10 , height = 10)

#---------------checking for cells between EC and mural cells, if they are doublets
samples.filtered <- samples.filtered %>% FindClusters(resolution = 0.3)
DimPlot(samples.filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_cluster_res0.3 check for cl7.jpeg"), width=10 , height = 10)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_ncount_check for cl7_vln.jpeg"), width=10 , height =5)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_check for cl7_vln.jpeg"), width=5 , height =5)
#run all.markers and further analysis downstream

#after further analysis cluster 7 seemed to be doublets from endothelial and mural cells, therefore is removed
samples.filtered=subset(samples.filtered,idents = c("7"),invert=T) 
DimPlot(samples.filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_cluster_cl7_removed.jpeg"), width=10 , height = 10)

#new clusters after cl7 remove
samples.filtered <- samples.filtered %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()

#run all.markers and further analysis downstream
DimPlot(samples.filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_after_cl7_excl_before_cl5_excl.jpeg"), width=10 , height = 10)
#top 1-2 markers to show doublet features of cl5
all.markers = FindAllMarkers(samples.filtered, test.use = "MAST",min.pct = 0.3,assay = "RNA")
top.markers <-c("Dcn","Fabp4","Rgs5","Ifit1","Plp1","Kdr")
#vertikal
DotPlot(samples.filtered, features = top.markers,dot.scale = 10,assay = "RNA")+ coord_flip() + NoLegend()
ggsave(filename = paste0(sample,"_top_markers_per_cl_for_cl5-excl.jpeg"), width= 3.5, height = 2.7)

#----------after further analysis cluster 5 seemed to be doublets from endothelial and fibroblasts, therefore is removed-----------
samples.filtered=subset(samples.filtered,idents = c("5"),invert=T) 

#new clusters after cl5 remove
samples.filtered <- samples.filtered %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  identity()

#clustertree for optimazation of cluster resolution
sc.int=samples.filtered
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)
  
#-low resolution to only identify major celltypes
samples.filtered <- samples.filtered %>% FindClusters(resolution = 0.1)
DimPlot(samples.filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
    ggsave(filename = paste0(sample,"_clusters_res0.1.jpeg"), width=10 , height = 10)

#-------save filtered object -----------
saveRDS(samples.filtered,file = paste0(sample,"_filtered_final.rds"))
################################################----end of sample processing before analysis----##############################################################
##############################################################################################################################################################
########################################################----start of sample analysis----######################################################################
#----------load filtered object---------
Samples.combined <- readRDS(file ="~/Documents/FP_scRNA/R stuff/PVM harmony/analysis PVM_full integration v10 harmony/RDS/harmony_full_int_filtered_final.rds")

#-----------------------add sub_cell_type according to annotation-----------------------------------
#cl0 "Fibroblast", cl1 "Endothelial",cl2 "Mural",cl3 "Interferon",cl4 "Neuronal", cl5 "Lympathic EC", cl6 "Proliferating"
for (i in c(1:length(Samples.combined@meta.data$seurat_clusters))) {
    Samples.combined@meta.data$sub_cell_type[i]<-switch(as.numeric(Samples.combined@meta.data$seurat_clusters[i]), 
                                                          "Fibroblast", "Endothelial", "Mural", "Interferon", "Neuronal", "Lympathic EC", "Proliferating")
    }
Samples.combined$sub_cell_type = factor(Samples.combined$sub_cell_type,levels = c("Fibroblast", "Endothelial", "Mural", "Interferon", "Neuronal", "Lympathic EC", "Proliferating"))
DimPlot(Samples.combined,repel = T,group.by = "sub_cell_type",reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 10,cols = cluster.colors) + NoLegend()
ggsave(filename = paste0(sample,"_processed_clustering_res0.1_small.jpeg"), width=10 , height = 10)
    
#--------top10 & top5 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#result of findallmarkers
saveRDS(all.markers,file = paste0(sample,"all.markers_per_cluster_result.RDS"))
#load results if they were saved allready before
all.markers <-readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM harmony/analysis PVM_full integration v10 harmony/based on basic/harmony_full_intall.markers_per_cluster_result.RDS")
#get top10 markergenes per cluster
top10 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene))
#vertikal
DotPlot(Samples.combined, features = top10,dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"_top10_markers_per_cluster_v.jpeg"), width= 5, height = 20)
ggsave(filename = paste0(sample,"_top10_markers_per_cluster_v.svg"), width= 5, height = 20)
#horizontal
Idents(Samples.combined)= factor(Samples.combined@meta.data$seurat_clusters,levels(Idents(Samples.combined))[as.numeric(rev(levels(Idents(Samples.combined))))+1])
DotPlot(Samples.combined, features = rev(top10),dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave(filename = paste0(sample,"_top10_featureDOTplot_custom_cl_h.jpeg"), width=20 , height = 3)
ggsave(filename = paste0(sample,"_top10_featureDOTplot_custom_cl_h.svg"), width=20 , height = 3)
Idents(Samples.combined)=Samples.combined@meta.data$seurat_clusters
#top5 plot
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
top5[top5 %in% c("Gsn","Clec3b","Ly6c1","Egfl7")] = c("Col1a1","Pdgfra","Kdr","Pecam1") # exchanged some markers with well defined marker genes for clear annotation
#vertikal
DotPlot(Samples.combined,group.by = "sub_cell_type", features = top5,dot.scale = 10,assay = "RNA")+ coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.text = element_text(size = 15),axis.title.x = element_blank())
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.jpeg"), width=5 , height = 10.5)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.svg"), width=5 , height = 10.5)
#horizontal
Idents(Samples.combined)= factor(Samples.combined@meta.data$seurat_clusters,levels(Idents(Samples.combined))[as.numeric(rev(levels(Idents(Samples.combined))))+1])
DotPlot(Samples.combined, features = rev(top5),dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.jpeg"), width=12, height = 3.2)
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.svg"), width=12, height = 3.2)
Idents(Samples.combined)=Samples.combined@meta.data$seurat_clusters

# plot genotypes
DimPlot(Samples.combined, reduction = "umap", group.by = "orig.geno", pt.size = 1, cols = genotype.colors,order = c("NG2","Gli1","Myh11","Col1a1","Cdh5","PDGFRb")) + theme(legend.position = c(0.85, 0.9))
ggsave(filename = paste0(sample,"_processed_genotype.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_genotype.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", group.by = "orig.geno", pt.size = 0.5, cols = genotype.colors,order = c("NG2","Gli1","Myh11","Col1a1","Cdh5","PDGFRb")) + theme(legend.position = c(0.85, 0.9))
ggsave(filename = paste0(sample,"_processed_genotype_small.pt.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_genotype_small.pt.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap",group.by = "orig.geno",split.by = "orig.geno" , pt.size = 0.5, ncol = 3, cols = c(col.col1a1,col.pdgfrb,col.gli1,col.cdh5,col.myh11,col.ng2),order = c("NG2","Myh11","Cdh5","Gli1","PDGFRb","Col1a1")) + NoLegend()
ggsave(filename = paste0(sample,"_processed_genotype_split3x2.jpeg"), width=15 , height = 10)
ggsave(filename = paste0(sample,"_processed_genotype_split3x2.svg"), width=15 , height = 10)
DimPlot(Samples.combined, reduction = "umap",group.by = "orig.geno",split.by = "orig.geno" , pt.size = 0.5, ncol = 2, cols = c(col.col1a1,col.cdh5,col.pdgfrb,col.myh11,col.gli1,col.ng2),order = c("NG2","Gli1","Myh11","PDGFRb","Cdh5","Col1a1")) + NoLegend()
ggsave(filename = paste0(sample,"_processed_genotype_split2x3.jpeg"), width=10 , height = 15)
ggsave(filename = paste0(sample,"_processed_genotype_split2x3.svg"), width=10 , height = 15)
DimPlot(Samples.combined, reduction = "umap",group.by = "orig.geno",split.by = "orig.geno" , pt.size = 0.5, ncol = 6, cols = c(col.col1a1,col.pdgfrb,col.gli1,col.cdh5,col.myh11,col.ng2),order = c("NG2","Myh11","Cdh5","Gli1","PDGFRb","Col1a1")) + NoLegend()
ggsave(filename = paste0(sample,"_processed_genotype_split6x1.jpeg"), width=27 , height = 5)
ggsave(filename = paste0(sample,"_processed_genotype_split6x1.svg"), width=27 , height = 5)

# plot Sham vs TAC
DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 0.5, cols = c('#00C5C080','#FA756C80'),order = c("TAC","Sham")) + theme(legend.position = c(0.8, 0.95))
ggsave(filename = paste0(sample,"_processed_Sham_vs_TAC.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_Sham_vs_TAC.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 0.5, cols = c('#00C5C080','#FA756C80'),split.by = "stim") + NoLegend()
ggsave(filename = paste0(sample,"_processed_Sham_vs_TAC_split.jpeg"), width=20 , height = 10)
ggsave(filename = paste0(sample,"_processed_Sham_vs_TAC_split.svg"), width=20 , height = 10)

#plot clustering res0.1 with nice colours
DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12,cols = cluster.colors) + NoLegend()
ggsave(filename = paste0(sample,"_processed_clustering_res0.1.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_clustering_res0.1.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 12,cols = cluster.colors) + NoLegend()
ggsave(filename = paste0(sample,"_processed_clustering_res0.1_small.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_clustering_res0.1_small.svg"), width=10 , height = 10)

#ploting for tdTom in the final clusters
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5,order = T)
ggsave(filename = paste0(sample,"_tdTom.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_tdTom_Vln.jpeg"), width=10 , height =10,limitsize = F)

#-----bar graph per genotype
dittoBarPlot(Samples.combined, "orig.geno", group.by = "sub_cell_type",color.panel = genotype.colors, var.labels.reorder = c(6,1,2,4,3,5),x.reorder = c(2,1,5,3,6,4,7)) + theme(axis.text.x = element_text(angle = 45,hjust = 1),plot.title = element_blank(),axis.title.x = element_blank())
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster.jpeg"), width=5 , height =5)

#--------ploting of multiple featues and save individual------------------------------
features_VSMC = c("Myh11","Acta2","Tagln","Cnn1","Myocd","Vcl","Smtn","Ly6a")
features_peri = c("Cspg4","Pdgfrb","Pdgfra","Colec11","Kcnj8","Abcc9","Vtn","Notch3")
features_fibro = c("Dcn","Lum","Cxcl1","Col1a1","Apoe","Ogn","Postn","Gli1")
features_endo = c("Cdh5","Pecam1","Kdr","Fabp4","Cd36","Srgn")
features_all = c(features_endo,features_fibro,features_peri,features_VSMC,"end")
allfp = vector("list",1)
dir.create("multiple_marker_plots")
for (f in features_all) {
  if (f=="end") {
    print(ggarrange(plotlist = allfp, widths = 30, heights = 30))
    ggsave(filename = paste0("multiple_marker_plots/",sample,"all_features.jpeg"), width=30 , height =30)
  }
  FeaturePlot(Samples.combined, features = f, min.cutoff = "q9",pt.size = 0.5,order = T) 
  ggsave(filename = paste0("multiple_marker_plots/",sample,"_",f,".jpeg"), width=10 , height = 10)
  allfp[[f]]<- FeaturePlot(Samples.combined, features = f, min.cutoff = "q9",pt.size = 0.5,order = T)+NoLegend()
  
  VlnPlot(Samples.combined, features = f,pt.size = 0)
  ggsave(filename = paste0("multiple_marker_plots/",sample,"_",f,"_Vln.jpeg"), width=10 , height = 10)
}

#plot ECM between genotype (run ECM scoring first)
Samples.combined$orig.geno = factor(Samples.combined$orig.geno,levels = c("Col1a1","PDGFRb","Gli1","Cdh5","Myh11","NG2"))
print(VlnPlot(Samples.combined, features = "Core_matrisome1", pt.size = 0,group.by = "orig.geno", cols = c(col.col1a1,col.pdgfrb,col.gli1,col.cdh5,col.myh11,col.ng2)) +
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0),axis.text.x.bottom = element_text(hjust = 0.5))+ NoLegend())
ggsave(filename = paste0("ECM scoring output/",sample,"_ECMscore_genotypes_vln.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0("ECM scoring output/",sample,"_ECMscore_genotypes_vln.svg"), width=10 , height = 10) 
#plot ECM between subclusters
print(VlnPlot(Samples.combined, features = "Core_matrisome1", pt.size = 0,group.by = "seurat_clusters",cols = cluster.colors) +
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0),axis.text.x.bottom = element_text(hjust = 0.5))+ NoLegend())
ggsave(filename = paste0("ECM scoring output/",sample,"_ECMscore_seurat_clusters_vln.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0("ECM scoring output/svg/",sample,"_ECMscore_seurat_clusters_vln.svg"), width=10 , height = 10) 