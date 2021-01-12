#script for analysis of Ng2 Sham and TAC integrated
library(Seurat)
library(cowplot)
library(ggplot2)
library(clustree)
library(genesorteR, quietly = TRUE)
library(writexl)
library(dplyr)
library(ggpubr)
setwd("~/Documents/FP_scRNA/R stuff/")
sample="Ng2_integrated_"
dir.create(paste0(sample,"pairwise"))

#----function to perfrom clustering after every filterstep------
recluster <- function(object,res){
  DefaultAssay(object) <- "integrated"
  # Run the standard workflow for visualization and clustering
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, verbose = FALSE)
  # t-SNE and Clustering
  object <- RunUMAP(object, reduction = "pca", dims = 1:30,repulsion.strength = 5)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object, resolution = res)
  return(object)
}
################################################----start of sample processing before analysis----##################################################
Samples.combined = readRDS(file = "~/Documents/FP_scRNA/R stuff/analysis PVM_v8/RDS_basic_filtering_Sham&TAC_integrated/NG2_Sham&TAC_integrated_seurat_object_after_bacis_filters.rds")
Samples.combined = recluster(Samples.combined,0.5)
dir.create(paste0(sample,"pairwise/processing_steps/"))

#plot before tdTom filtering steps
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 2, cols = c('#FA756C80','#00C5C080')) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"pairwise/processing_steps/",sample,"before_processing_clusters_res0.5.jpeg"), width=20 , height = 10)

#-------------ploting only tdtom
#quality check on tdTom expression
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(sample,"pairwise/processing_steps/",sample,"tdTom.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(sample,"pairwise/processing_steps/",sample,"tdTom_Vln.jpeg"), width=10 , height =10,limitsize = F)
RidgePlot(Samples.combined,features = "tdTomato-Green",group.by = "orig.ident",slot = "counts",assay = "RNA")
ggsave(filename = paste0(sample,"pairwise/processing_steps/",sample,"tdTom_reads.jpeg"), width=10 , height = 10)

#run Genesorter to get GeneProb of tdtom per cluster
sg = sortGenes(Samples.combined@assays$RNA@data, Idents(Samples.combined))
colnames(sg$condGeneProb) = paste0(levels(Idents(Samples.combined)))
Samples.combined_geneProb = as.data.frame(sg$condGeneProb)["tdTomato-Green",]
Samples.combined_geneProb
remove(sg)

#----------subsetting to filter out tdtom low clusters-----------
Samples.combined=subset(Samples.combined,idents = c("0","1","5"))
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined=subset(Samples.combined,cells = WhichCells(Samples.combined,slot = "counts",expression = `tdTomato-Green` == 0),invert=T)
Samples.combined = recluster(Samples.combined,0.5)

#plot after tdTom filtering step and before nfeature filter
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 2, cols = c('#FA756C80','#00C5C080')) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"pairwise/processing_steps/",sample,"after_tdTom_low_remove_clusters_res0.5.jpeg"), width=20 , height = 10)

#------------basic nfeature and ncount plots------------------
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"pairwise/processing_steps/",sample,"_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"pairwise/processing_steps/",sample,"_nfeature_ncount_vln.jpeg"), width=10 , height =5)

#----------subsetting to filter out feature count low clusters-----------
#no cluster with low features in this dataset

#plot after processing
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 2, cols = c('#FA756C80','#00C5C080')) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5.jpeg"), width=20 , height = 10)
#no label version
DimPlot(Samples.combined, reduction = "umap",  pt.size = 2) + NoLegend()
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5_noLabel.jpeg"), width=10 , height = 10)

#-------normalize and scale RNA data------------------------
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

###---------------seperating Sham and Tac bac to the individual samples-----------
Idents(Samples.combined)=Samples.combined@meta.data$orig.ident
Sham = subset(Samples.combined,idents = "FP23_NG2_Sham")
saveRDS(Sham,file = paste0(sample,"pairwise/",sample,"_Sham_indiv_after_pairwise_integration.rds"))
TAC = subset(Samples.combined,idents = "FP22_NG2_TAC")
saveRDS(TAC,file = paste0(sample,"pairwise/",sample,"_TAC_indiv_after_pairwise_integration.rds"))