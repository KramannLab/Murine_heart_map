#script for the integration and analysis of all Sham and TAC samples based on only the basic filtered RDS files with harmony
library(harmony)#http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html
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
library(reshape2)
source("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/scripts_for_upload/helper_functions_DEG_GO_ploting.R")
set.seed(111)

setwd("~/Documents/FP_scRNA/R stuff/PVM harmony/integration_new_28d/")
sample="harmony_full_int"


#colorblind friendly color panels for all figures
#for conditions/stim #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
sham.col = Okabe_Ito[1]
tac.col = Okabe_Ito[2]
tac28.col = Okabe_Ito[3]
stim.col=c(sham.col,tac.col,tac28.col)
stim.col.light=c("#E69F0080","#56B4E980","#009E7380")
#for genotypes #From Paul Tol: https://personal.sron.nl/~pault/ '#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'
col.pdgfrb="#EE6677"
col.myh11="#CCBB44"
col.ng2="#66CCEE"
col.col1a1="#228833"
col.gli1="#4477AA"
col.cdh5="#AA3377"
geno.cl = c(col.col1a1,col.pdgfrb,col.gli1,col.cdh5,col.myh11,col.ng2)
geno.cl.light = c("#22883380","#EE667780","#4477AA80","#AA337780","#CCBB4480","#66CCEE80")
#for subclustering #RColorBrewer
display.brewer.pal(name = "Paired",n = 12)
cluster.col=brewer.pal(name = "Paired",n=12)
cluster.col=cluster.col[c(2,4,6,8,10,12,1)]
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)


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
Gli1_TAC_28 = "FP71_Gli1_TAC_28d"
Gli1_Sham=readRDS(file = paste0(datapath, Gli1_Sham , "_seurat_object_after_bacis_filters.rds"))
Gli1_TAC=readRDS(file = paste0(datapath, Gli1_TAC , "_seurat_object_after_bacis_filters.rds"))
Gli1_TAC_28=readRDS(file = paste0(datapath, Gli1_TAC_28 , "_seurat_object_after_bacis_filters.rds"))

#load for ------Col1a1------ samples
Col1a1_Sham = 'FP19_Col1a1_Sham'
Col1a1_TAC = 'FP18_Col1a1_TAC'
Col1a1_Sham=readRDS(file = paste0(datapath, Col1a1_Sham , "_seurat_object_after_bacis_filters.rds"))
Col1a1_TAC=readRDS(file = paste0(datapath, Col1a1_TAC , "_seurat_object_after_bacis_filters.rds"))

#load for ------Cdh5------ samples
Cdh5_Sham = 'FP21_Cdh5_Sham'
Cdh5_TAC = 'FP20_Cdh5_TAC'
Cdh5_TAC_28 = 'FP72_Cdh5_TAC_28d'
Cdh5_Sham=readRDS(file = paste0(datapath, Cdh5_Sham , "_seurat_object_after_bacis_filters.rds"))
Cdh5_TAC=readRDS(file = paste0(datapath, Cdh5_TAC , "_seurat_object_after_bacis_filters.rds"))
Cdh5_TAC_28=readRDS(file = paste0(datapath, Cdh5_TAC_28 , "_seurat_object_after_bacis_filters.rds"))

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
samples.merged <- merge(x = PDGFRb_Sham, y =  list(PDGFRb_TAC, Cdh5_TAC,Cdh5_TAC_28, Gli1_TAC,Gli1_TAC_28, Col1a1_TAC, NG2_TAC, Myh11_TAC,Cdh5_Sham, Gli1_Sham, Col1a1_Sham, NG2_Sham, Myh11_Sham), project = "PVM")
samples.merged$orig.ident = factor(samples.merged$orig.ident,levels = c("FP19_Col1a1_Sham","FP18_Col1a1_TAC","FP3_PDGFRb_Sham","FP14_PDGFRb_TAC","FP13_Gli1_Sham","FP12_Gli1_TAC","FP71_Gli1_TAC_28d","FP21_Cdh5_Sham","FP20_Cdh5_TAC","FP72_Cdh5_TAC_28d","FP56_Myh11_Sham","FP57_Myh11_TAC","FP22_NG2_TAC","FP23_NG2_Sham"))

#first processing
samples.merged <- samples.merged %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = samples.merged@var.genes, npcs = 20, verbose = FALSE)

###---subset the individual stim with just TAC or Sham--###
samples.merged[['stim']] <- ifelse(samples.merged@meta.data$stim %in% c("FP12_Gli1_TAC","FP22_NG2_TAC","FP14_PDGFRb_TAC","FP18_Col1a1_TAC","FP20_Cdh5_TAC","FP57_Myh11_TAC"),'TAC',samples.merged@meta.data$stim)
samples.merged[['stim']] <- ifelse(samples.merged@meta.data$stim %in% c("FP3_PDGFRb_Sham","FP21_Cdh5_Sham", "FP13_Gli1_Sham", "FP19_Col1a1_Sham", "FP23_NG2_Sham", "FP56_Myh11_Sham"),'Sham',samples.merged@meta.data$stim)
samples.merged[['stim']] <- ifelse(samples.merged@meta.data$stim %in% c("FP72_Cdh5_TAC_28d","FP71_Gli1_TAC_28d"),'TAC_28',samples.merged@meta.data$stim)
samples.merged$stim = factor(samples.merged$stim ,levels = c("Sham","TAC","TAC_28"))

# add column for only the genotype
Gli1.ls = samples.merged@meta.data$orig.ident %in% c("FP12_Gli1_TAC","FP13_Gli1_Sham","FP71_Gli1_TAC_28d")
PDGFRb.ls = samples.merged@meta.data$orig.ident %in% c("FP14_PDGFRb_TAC","FP3_PDGFRb_Sham")
Cdh5.ls = samples.merged@meta.data$orig.ident %in% c("FP21_Cdh5_Sham","FP20_Cdh5_TAC","FP72_Cdh5_TAC_28d")
Col1a1.ls = samples.merged@meta.data$orig.ident %in% c("FP18_Col1a1_TAC","FP19_Col1a1_Sham")
NG2.ls = samples.merged@meta.data$orig.ident %in% c("FP22_NG2_TAC","FP23_NG2_Sham")
Myh11.ls = samples.merged@meta.data$orig.ident %in% c("FP56_Myh11_Sham","FP57_Myh11_TAC")
samples.merged@meta.data$orig.geno[Gli1.ls]="Gli1"
samples.merged@meta.data$orig.geno[PDGFRb.ls]="Pdgfrb"
samples.merged@meta.data$orig.geno[Cdh5.ls]="Cdh5"
samples.merged@meta.data$orig.geno[Col1a1.ls]="Col1a1"
samples.merged@meta.data$orig.geno[NG2.ls]="Ng2"
samples.merged@meta.data$orig.geno[Myh11.ls]="Myh11"
samples.merged$orig.geno = factor(samples.merged$orig.geno, levels = c("Col1a1","Pdgfrb","Gli1","Cdh5","Myh11","Ng2"))

#see batch effect before harmony
samples.merged.raw <- samples.merged %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
samples.merged.raw$geno.stim = paste0(samples.merged.raw$orig.geno,"_",samples.merged.raw$stim)
DimPlot(object = samples.merged.raw, reduction = "umap", pt.size = .1, group.by = "geno.stim")+NoLegend()
ggsave(filename = paste0(sample,"_raw_data_just_merged.jpeg"),width = 10,height = 10)
DimPlot(object = samples.merged.raw, reduction = "umap", pt.size = .1, group.by = "geno.stim")
ggsave(filename = paste0(sample,"_raw_data_just_merged_v2.jpeg"),width = 10,height = 10)
DimPlot(object = samples.merged.raw, reduction = "umap", pt.size = .1, split.by = "geno.stim",ncol=3)
ggsave(filename = paste0(sample,"_raw_data_just_merged_split.jpeg"),width = 15,height = 20)
DimPlot(object = samples.merged.raw, reduction = "umap", pt.size = .1, split.by = "orig.geno",ncol=3)
ggsave(filename = paste0(sample,"_raw_data_just_merged_split_geno.jpeg"))
saveRDS(samples.merged.raw,file = paste0(sample,"_merged_raw.rds"))
samples.merged.raw <- readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM harmony/integration_new_28d/RDS/harmony_full_int_merged_raw.rds")

#run harmony on the data
#start processing
samples.merged <- samples.merged %>% RunHarmony("orig.ident", plot_convergence = TRUE,epsilon.cluster = -Inf)
#process based on harmony result
samples.merged <- samples.merged %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
saveRDS(samples.merged,file = paste0(sample,"_merged_harmony.rds"))
samples.merged <- readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM harmony/integration_new_28d/RDS/harmony_full_int_merged_harmony.rds")
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
DimPlot(myh11.subset, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()+ xlim(-11,8) + ylim(-6,12)
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove.jpeg", width=10 , height = 10)
FeaturePlot(myh11.subset, features = c("tdTomato-Green"),pt.size = 2)+ xlim(-11,8) + ylim(-6,12)
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove_tdTom.jpeg", width=10 , height = 10)
VlnPlot(myh11.subset, features = c("tdTomato-Green"), pt.size = 0)+NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove_tdTom_vln.jpeg", width=10 , height = 10)
table(myh11.subset$seurat_clusters) #remove cluster with less than 20 cells 
myh11.subset = subset(myh11.subset,idents = c("0","1","3","4","6","7","15","18","20","21","22","23"),invert=T)
myh11.subset$seurat_clusters=Idents(myh11.subset)
DimPlot(myh11.subset, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()+ xlim(-11,8) + ylim(-6,12)
ggsave(filename = "FP56_Myh11_Sham_before_tdTom-remove_lower20removed.jpeg", width=10 , height = 10)
sg = sortGenes(myh11.subset@assays$RNA@data, Idents(myh11.subset))
myh11.subset$tdTom.score=1
i=1
for (cluster in levels(factor(myh11.subset$seurat_clusters))) {myh11.subset$tdTom.score[myh11.subset$seurat_clusters==cluster]=round((sg$condGeneProb)["tdTomato-Green",i],digits = 3);i=i+1}
Idents(myh11.subset)="tdTom.score"
DimPlot(myh11.subset, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12,label.box = T,repel = T) + NoLegend() + xlim(-11,8) + ylim(-6,12)
ggsave(filename = "Only_FP56_Myh11_Sham_before_tdTom-remove_td-score.jpeg", width=10 , height = 10)
Idents(myh11.subset)="seurat_clusters"
myh11.subset = subset(myh11.subset,idents = c("2","5","9","10","11","13","14","17"),invert=T)
FeaturePlot(myh11.subset, features = c("tdTomato-Green"),pt.size = 2,order = T) + xlim(-11,8) + ylim(-6,12)
ggsave(filename = "Only_FP56_Myh11_Sham_tdTom-removed.jpeg", width=10 , height = 10)

#---run specialized function to filter out lowtdtomato cluster contributions per sample, clusters with less than 20 cells and 0 tdtom reads cells
source("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/filterLowTdtomCl_function_harmony.R")
samples.merged.td.filter<- filterLowTdtomCl.harmony(samples.merged,sample)
saveRDS(samples.merged.td.filter,file = paste0(sample,"_merged_td_filter.rds"))
samples.merged.td.filter <- readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM harmony/integration_new_28d/RDS/harmony_full_int_merged_td_filter.rds")

#-------------ploting only tdtom
#quality check on tdTom expression after filtering
FeaturePlot(samples.merged.td.filter, features = c("tdTomato-Green"),pt.size = 0.5,order = T)
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
VlnPlot(samples.merged.td.filter,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_before_vln.jpeg"), width=10 , height =10)
FeaturePlot(samples.merged.td.filter, features = c('nFeature_RNA'),pt.size = .5)+theme(legend.text = element_text(size=20))
ggsave(filename = paste0(sample,"_nfeature_before.jpeg"), width=10 , height = 10)

#----------subsetting to filter out feature count low or high clusters (low quality cells and potential doublets)-----------
samples.filtered=subset(samples.merged.td.filter,idents = c("13","20"),invert=T) 
#after
FeaturePlot(samples.filtered, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount_after.jpeg"), width=10 , height = 5)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount_after_vln.jpeg"), width=10 , height =5)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_after_vln.jpeg"), width=10 , height =10)
FeaturePlot(samples.filtered, features = c('nFeature_RNA'),pt.size = .5)+theme(legend.text = element_text(size=20))
ggsave(filename = paste0(sample,"_nfeature_after.jpeg"), width=10 , height = 10)

#new clusters after nfeature remove
samples.filtered <- samples.filtered %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  identity()

#set resolution for analysis
samples.filtered <- samples.filtered %>% FindClusters(resolution = 0.1)
DimPlot(object = samples.filtered, reduction = "umap", pt.size = .1,label=T,label.size = 10)+NoLegend()  
ggsave(filename = paste0(sample,"_clusters_res0.1.jpeg"), width=10 , height = 10)

#----remove artifact cluster6------
samples.filtered <- subset(samples.filtered,ident=6,invert=T)

#---------------checking for cells between EC and mural cells, if they are doublets
samples.filtered <- samples.filtered %>% FindNeighbors(reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.6)
DimPlot(samples.filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_cluster_res0.6 check for cl11.jpeg"), width=10 , height = 10)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_ncount_check for cl11_vln.jpeg"), width=10 , height =5)
VlnPlot(samples.filtered,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0) +NoLegend()+theme(axis.text.x = element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"_nfeature_check for cl11_vln.jpeg"), width=5 , height =5)
#run all.markers and further analysis downstream

#after further analysis cluster 11 seemed to be doublets from endothelial and mural cells, therefore is removed
samples.filtered=subset(samples.filtered,idents = c("11"),invert=T)
DimPlot(samples.filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_cluster_cl11_removed.jpeg"), width=10 , height = 10)

#new clusters after cl11 remove
samples.filtered <- samples.filtered %>% 
  RunUMAP(reduction = "harmony", dims = 1:20,repulsion.strength = 5) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()

#-------normalize and scale RNA data------------------------
DefaultAssay(samples.filtered) <- "RNA"
Samples.combined <- NormalizeData(samples.filtered, verbose = FALSE)
Samples.combined <- ScaleData(samples.filtered, verbose = FALSE)

#-------save filtered object -----------
saveRDS(samples.filtered,file = paste0(sample,"_filtered_final.rds"))
################################################----end of sample processing before analysis----##############################################################
##############################################################################################################################################################
########################################################----start of sample analysis----######################################################################
#----------load filtered object---------
Samples.combined <- readRDS(file ="~/Documents/FP_scRNA/R stuff/PVM harmony/integration_new_28d/RDS/harmony_full_int_filtered_final.rds")

#--------annotation for the major celltypes based on marker genes and original genotype lineage tracing------
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Fibroblast","1"="Endothelial","2"="Mural","3"="Interferon","4"="Schwann cells","5"="Lymphatic EC","6"="Proliferating EC")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes <- Idents(Samples.combined)
cell.levels = c("Fibroblast","Endothelial","Mural","Interferon","Schwann cells","Lymphatic EC","Proliferating EC")
Samples.combined$celltypes = factor(Samples.combined$celltypes,levels = cell.levels)
#short abroviation of annotation
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Fib","1"="EC","2"="Mural","3"="Int","4"="Sw","5"="LymEC","6"="Prol")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes.short <- Idents(Samples.combined)
cell.levels = c("Fib","EC","Mural","Int","Sw","LymEC","Prol")
Samples.combined$celltypes.short = factor(Samples.combined$celltypes.short,levels = cell.levels)
Idents(Samples.combined)="seurat_clusters"
#save after renaming
saveRDS(Samples.combined,file = paste0(sample,"_filtered_final.rds"))

#--------top10 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#top5 plot
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
top5[top5 %in% c("Mfap5","Smoc2","Ly6c1","Gpihbp1","Steap4")] = c("Col1a1","Pdgfra","Kdr","Pecam1","Colec11") # exchanged some markers with well defined marker genes for clear annotation
#vertikal
DotPlot(Samples.combined,group.by = "celltypes", features = top5,dot.scale = 10,assay = "RNA")+ coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.text.y = element_text(face = "italic"),axis.text = element_text(size = 15),axis.title.x = element_blank()) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.jpeg"), width=5.5 , height = 12)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.svg"), width=5.5 , height = 12)

#---------------------------------umap plots---------------------------
# plot genotypes
DimPlot(Samples.combined, reduction = "umap", group.by = "orig.geno", pt.size = 1, cols = geno.cl) + theme(legend.position = c(0.85, 0.9))
ggsave(filename = paste0(sample,"_processed_genotype.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_genotype.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", group.by = "orig.geno", pt.size = 0.5, cols = geno.cl) + theme(legend.position = c(0.85, 0.9))
ggsave(filename = paste0(sample,"_processed_genotype_small.pt.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_genotype_small.pt.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap",group.by = "orig.geno",split.by = "orig.geno" , pt.size = 0.5, ncol = 2, cols = c(col.col1a1,col.cdh5,col.pdgfrb,col.myh11,col.gli1,col.ng2),order = c("Ng2","Gli1","Myh11","Pdgfrb","Cdh5","Col1a1")) + NoLegend()
ggsave(filename = paste0(sample,"_processed_genotype_split2x3.jpeg"), width=10 , height = 15)
ggsave(filename = paste0(sample,"_processed_genotype_split2x3.svg"), width=10 , height = 15)
# plot Sham vs TAC
DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 0.5, cols = stim.col.light) + theme(legend.position = c(0.8, 0.9))
ggsave(filename = paste0(sample,"_processed_Sham_vs_TAC.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_Sham_vs_TAC.svg"), width=10 , height = 10)
#plot clustering res0.1 with nice colours
DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12,cols = cluster.col) + NoLegend()
ggsave(filename = paste0(sample,"_processed_clustering_res0.1.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_processed_clustering_res0.1.svg"), width=10 , height = 10)
#ploting for tdTom in the final clusters
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5,order = T)+ scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_tdTom.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined,group.by = "celltypes.short", features = c("tdTomato-Green"), pt.size = 0,cols = cluster.col) +NoLegend()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = paste0(sample,"_tdTom_Vln.jpeg"), width=5 , height =5,limitsize = F)
#---plot only the mural cell cluster with pericyte and VSMCs markers
mural_cells = subset(Samples.combined,idents = 2)
#vsmc marker
FeaturePlot(mural_cells, features = c("Myh11"),pt.size = 1,order = T) + xlim(-9,2) + ylim(12,24)+ scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_mural-cluster_Myh11.jpeg"), width=10 , height =10,limitsize = F)
FeaturePlot(mural_cells, features = c("Acta2"),pt.size = 1,order = T) + xlim(-9,2) + ylim(12,24)+ scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_mural-cluster_Acta2.jpeg"), width=10 , height =10,limitsize = F)
FeaturePlot(mural_cells, features = c("Tagln"),pt.size = 1,order = T) + xlim(-9,2) + ylim(12,24)+ scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_mural-cluster_Tagln.jpeg"), width=10 , height =10,limitsize = F)
#peri marker
FeaturePlot(mural_cells, features = c("Kcnj8"),pt.size = 1,order = T) + xlim(-9,2) + ylim(12,24)+ scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_mural-cluster_Kcnj8.jpeg"), width=10 , height =10,limitsize = F)
FeaturePlot(mural_cells, features = c("Abcc9"),pt.size = 1,order = T) + xlim(-9,2) + ylim(12,24)+ scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_mural-cluster_Abcc9.jpeg"), width=10 , height =10,limitsize = F)
FeaturePlot(mural_cells, features = c("Colec11"),pt.size = 1,order = T) + xlim(-9,2) + ylim(12,24)+ scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_mural-cluster_Colec11.jpeg"), width=10 , height =10,limitsize = F)

#-----------------------bar graph per genotype--------------------------------
dittoBarPlot(Samples.combined, "orig.geno", group.by = "celltypes",color.panel = geno.cl,var.labels.reorder = c(2,6,3,1,4,5),x.reorder = c(2,1,5,3,7,4,6)) + theme(axis.text.x = element_text(angle = 45,hjust = 1),plot.title = element_blank(),axis.title.x = element_blank())
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster.jpeg"), width=5 , height =7.5)
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster.svg"))

#--------ploting of multiple featues and save individual------------------------------
FeaturePlot(Samples.combined, features = "Col1a1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste0("multiple_marker_plots/",sample,"_Col1a1.jpeg"), width=10 , height = 10)
FeaturePlot(Samples.combined,features = "Pdgfrb",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Pdgfrb.jpeg"), width=10 , height = 10)
FeaturePlot(Samples.combined,features = "Gli1",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Gli1.jpeg"), width=10 , height = 10)
FeaturePlot(Samples.combined,features = "Cdh5",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Cdh5.jpeg"), width=10 , height = 10)
FeaturePlot(Samples.combined,features = "Myh11",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Myh11.jpeg"), width=10 , height = 10)
FeaturePlot(Samples.combined,features = "Cspg4",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Cspg4.jpeg"), width=10 , height = 10)

#-------------ECM scoring-----------------------------
Samples.combined <- scoreECM(Samples.combined)
#plot ECM between genotype (run ECM scoring first)
VlnPlot(Samples.combined, features = "Core_matrisome1", pt.size = 0,group.by = "orig.geno", cols = c(col.col1a1,col.pdgfrb,col.gli1,col.cdh5,col.myh11,col.ng2)) +
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0),axis.text.x.bottom = element_text(hjust = 0.5))+ NoLegend()
ggsave(filename = paste0(sample,"_ECMscore_genotypes_vln.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0(sample,"_ECMscore_genotypes_vln.svg"), width=10 , height = 10) 
#plot ECM between subclusters split
VlnPlot(Samples.combined, features = "Core_matrisome1", pt.size = 0,group.by = "celltypes.short",split.by = "stim",cols = stim.col) +
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0),axis.text.x.bottom = element_text(hjust = 0.5))+ NoLegend()
ggsave(filename = paste0(sample,"_ECMscore_seurat_clusters_vln_split.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0(sample,"_ECMscore_seurat_clusters_vln_split.svg"), width=10 , height = 10) 
#feature plot matrisome score for fig 3
FeaturePlot(Samples.combined, features = "Core_matrisome1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_ECMscore_f.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0(sample,"_ECMscore_f.svg"), width=10 , height = 10) 
#plot collagen between subclusters split
VlnPlot(Samples.combined, features = "Collagens1", pt.size = 0,group.by = "celltypes.short",split.by = "stim",cols = stim.col) +
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0),axis.text.x.bottom = element_text(hjust = 0.5))+ NoLegend()
ggsave(filename = paste0(sample,"_Collagen_1score_seurat_clusters_vln_split.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0(sample,"_Collagen_1score_seurat_clusters_vln_split.svg"), width=10 , height = 10) 
#feature plot collagen score for supp 3
FeaturePlot(Samples.combined, features = "Collagens1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_Collagen_1score_f.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0(sample,"_Collagen_1score_f.svg"), width=10 , height = 10) 
#plot Proteoglycan between subclusters split
VlnPlot(Samples.combined, features = "Proteoglycans1", pt.size = 0,group.by = "celltypes.short",split.by = "stim",cols = stim.col) +
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0),axis.text.x.bottom = element_text(hjust = 0.5))+ NoLegend()
ggsave(filename = paste0(sample,"_Proteoglycans1score_seurat_clusters_vln_split.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0(sample,"_Proteoglycans1score_seurat_clusters_vln_split.svg"), width=10 , height = 10) 
#feature plot Proteoglycan score for supp 3
FeaturePlot(Samples.combined, features = "Proteoglycans1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_Proteoglycans1score_f.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0(sample,"_Proteoglycans1score_f.svg"), width=10 , height = 10) 

#statistic testing for the ECM scorings
Samples.combined$celltypes.stim = paste0(Samples.combined$celltypes.short,"_",Samples.combined$stim)
#test between celltypes/conditions
ECM.score.data <- FetchData(Samples.combined,vars = c("Core_matrisome1","celltypes.stim"))
#test overall
k.result <- kruskal.test(ECM.score.data$Core_matrisome1~ECM.score.data$celltypes.stim)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(ECM.score.data$Core_matrisome1,ECM.score.data$celltypes.stim,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result ECMscore.xlsx"))
#test between genotypes
ECM.score.data <- FetchData(Samples.combined,vars = c("Core_matrisome1","orig.geno"))
#test overall
k.result <- kruskal.test(ECM.score.data$Core_matrisome1~ECM.score.data$orig.geno)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(ECM.score.data$Core_matrisome1,ECM.score.data$orig.geno,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result ECMscore_genotypes.xlsx"))

col.score.data <- FetchData(Samples.combined,vars = c("Collagens1","celltypes.stim"))
#test overall
k.result <- kruskal.test(col.score.data$Collagens1~col.score.data$celltypes.stim)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(col.score.data$Collagens1,col.score.data$celltypes.stim,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result Collagens1score.xlsx"))

proteo.score.data <- FetchData(Samples.combined,vars = c("Proteoglycans1","celltypes.stim"))
#test overall
k.result <- kruskal.test(proteo.score.data$Proteoglycans1~proteo.score.data$celltypes.stim)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(proteo.score.data$Proteoglycans1,proteo.score.data$celltypes.stim,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result Proteoglycans1score.xlsx"))

