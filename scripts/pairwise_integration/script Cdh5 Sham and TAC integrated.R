#script for analysis of Cdh5 Sham and TAC integrated
library(Seurat)
library(cowplot)
library(ggplot2)
library(clustree)
library(genesorteR, quietly = TRUE)
library(writexl)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(readr)
setwd("~/Documents/FP_scRNA/R stuff/")
sample="Cdh5_integrated_"
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
Samples.combined = readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering_Sham&TAC_integrated/Cdh5_Sham&TAC_integrated_seurat_object_after_bacis_filters.rds")
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
Samples.combined=subset(Samples.combined,idents = c("11","12"),invert=T)
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
DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 2, cols = c('#FA756C80','#00C5C080')) +theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5_tac_vs_sham.svg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5_tac_vs_sham.jpeg"), width=10 , height = 10)
#no label version
DimPlot(Samples.combined, reduction = "umap",  pt.size = 2) + NoLegend()
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5_noLabel.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"pairwise/",sample,"after_processing_clusters_res0.5_noLabel.svg"), width=10 , height = 10)

#-------normalize and scale RNA data------------------------
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

#--------save filtered object-----------
saveRDS(Samples.combined,file = paste0(sample,"pairwise/",sample,"_filtered_processed.rds"))
################################################----end of sample processing before analysis----##################################################

########################################################----start of sample analysis----###########################################################
#----------load filtered object---------
Samples.combined <- readRDS(file = "~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Cdh5_integrated_pairwise/Cdh5_integrated__filtered_processed.rds")

#in case of rename of stim is requiered
table(Samples.combined$stim)
Samples.combined$stim[Samples.combined$stim=="FP21_Cdh5_Sham"] <- "Sham"
Samples.combined$stim[Samples.combined$stim=="FP20_Cdh5_TAC"] <- "TAC"

#clustertree for optimazation of cluster resolution
sc.int=Samples.combined
DefaultAssay(sc.int)<-"integrated"
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"pairwise/",sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)

#--------top5 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
top5 <- c(top5,"Npr3","Ctsh") #added Npr3 and Ctsh for 5 markers for endo
#vertikal
DotPlot(Samples.combined, features = rev(top5),dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"pairwise/",sample,"_top5_markers_per_cluster.jpeg"), width=6 , height = 12)
ggsave(filename = paste0(sample,"pairwise/",sample,"_top5_markers_per_cluster.svg"), width=6 , height = 12)

#horizontal
Idents(Samples.combined)= factor(Samples.combined@meta.data$seurat_clusters,levels(Idents(Samples.combined))[as.numeric(rev(levels(Idents(Samples.combined))))+1])
Samples.combined$sub_cell_type_short = factor(Samples.combined$sub_cell_type_short,levels = rev(celltypes))
DotPlot(Samples.combined,group.by = "sub_cell_type_short", features =top5,dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.jpeg"), width=16 , height = 4.3)
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.svg"), width=16 , height = 4.3)
Idents(Samples.combined)=Samples.combined@meta.data$seurat_clusters

#vln plot for the endocard markers
v1 <- VlnPlot(Samples.combined,group.by = "sub_cell_type_short",features = c("Npr3"),pt.size = 0)+NoLegend() + theme(axis.text.x = element_text(angle = 45,hjust =1),axis.title.x = element_blank())
v2 <- VlnPlot(Samples.combined,group.by = "sub_cell_type_short",features = c("Vwf"),pt.size = 0)+NoLegend() + theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.title.x = element_blank())
ggarrange(v1,v2,ncol = 1)
ggsave(filename = paste0(sample,"_Npr3_and_Vwf_for_endocard_cluster10.jpeg"), width=5 , height = 5)
ggsave(filename = paste0(sample,"_Npr3_and_Vwf_for_endocard_cluster10.svg"), width=5 , height = 5)
#generating matrix with genesorteR
source("~/Documents/FP_scRNA/R stuff/scripts/methods/genesorteR_FP.R") #change date of scripts folder
doTheGenesorteR(Samples.combined,paste0(sample,"pairwise/",sample))

#------advaced bar garph----Sham vs. TAC--------------------------------------------
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$stim) 
breakdown[,1]= round(breakdown[,1]*(sum(breakdown[,2])/sum(breakdown[,1]))) #normalize cellnumber of TAC to Sham
breakdown=t(breakdown)
breakdown <- apply(breakdown, 2, function(x){x*100/sum(x)})
breakdown = round(as.data.frame(breakdown),digits = 2)
breakdown.df = as.data.frame(breakdown)
breakdown.df = melt(t(breakdown.df))
celltypes<-c("Cap","CapA","CapV","StrEC","Ang","IntEC","Art","Lym","Rep","Cycl","Endo")

ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=c(celltypes,celltypes),y=1),hjust = 0, color="black", size=10)+
  coord_flip() + 
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = rev(levels(factor(breakdown.df$Var1))),name="Cluster")+
  scale_y_continuous(name="Percentage of total cells [%]") +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15)) + 
  NoLegend()+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_hline(yintercept=50, linetype="dashed", color = "black")
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.jpeg"), width=5 , height = 10)
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.svg"), width=5 , height = 10)

#--------heatmap of EC-atlas paper
EC_atlas_markers <- read_csv("~/Documents/FP_scRNA/R stuff/references/EC-atlas Heart heatmap table.csv")
markers=as.data.frame(EC_atlas_markers)

#scoring for EC map subtype
DefaultAssay(Samples.combined)="RNA"
Idents(Samples.combined)="seurat_clusters"
ctrl_genes = 35 #important
for (gset in levels(factor(markers$cell_type))){
  gene.list = list(markers$gene_name[markers$cell_type==gset])
  Samples.combined = AddModuleScore(object = Samples.combined, features = gene.list, name = gset, ctrl = ctrl_genes)
}
VlnPlot(Samples.combined,features = colnames(Samples.combined@meta.data[8:15]),pt.size = 0)
ggsave(filename = paste0(sample,"pairwise/",sample,"_EC-atlas_markers_as_Module_score.jpeg"), width=20, height =20)

avgexp = AverageExpression(Samples.combined, return.seurat = T)
DoHeatmap(avgexp,features = markers$gene_name,draw.lines = F,assay = "RNA",angle = 0)+ 
  scale_y_discrete(labels=rev(markers$cell_type))
ggsave(filename = paste0(sample,"pairwise/",sample,"_EC-atlas_markers_heatmap_RNA_averages.jpeg"), width=10, height =10)
ggsave(filename = paste0(sample,"pairwise/",sample,"_EC-atlas_markers_heatmap_RNA_averages.svg"), width=10, height =10)

DoHeatmap(Samples.combined,features = markers$gene_name,label = FALSE,assay = "RNA",angle = 0) + 
  scale_y_discrete(name="Top 20 markers from EC-atlas",labels=rev(markers$cell_type))
ggsave(filename = paste0(sample,"pairwise/",sample,"_EC-atlas_markers_heatmap_RNA.jpeg"), width=25, height =5)

DotPlot(Samples.combined, features = rev(markers[markers$cell_type=="large vein",1]),dot.scale = 10,assay = "RNA")+ coord_flip() + NoLegend()
ggsave(filename = paste0(sample,"pairwise/",sample,"_top20_large_vein_marker_cl10_endo.jpeg"), width=4.5 , height = 6)

##-------bar graph for percentage of cells per cluster (of total cells)------------------
breakdown<-round(prop.table(table(Samples.combined@meta.data$seurat_clusters))*100,digits = 1)
breakdown<-(as.data.frame(breakdown))
celltypes<-c("capillary EC","capillary artery EC","capillary vein EC","stressed EC","angiogenic EC","interferon EC","artery EC","lymphatic EC","DNA replicating EC","cycling EC","endocardial EC")
ggplot(data=breakdown, aes(x=Var1, y=Freq,fill=Var1))+
  geom_bar(stat="identity")+
  geom_text(aes(label=celltypes,y=1),hjust = 0, color="black", size=10)+
  coord_flip() + 
  scale_x_discrete(limits=rev(breakdown[,1])) + 
  scale_y_continuous(name="Percentage of total cells [%]") +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  NoLegend()
ggsave(filename = paste0(sample,"_barplot_percentage_per_cluster.jpeg"), width=5 , height = 10)
ggsave(filename = paste0(sample,"_barplot_percentage_per_cluster.svg"), width=5 , height = 10)

###---------------seperating Sham and Tac bac to the individual samples-----------
Idents(Samples.combined)=Samples.combined@meta.data$orig.ident
Sham = subset(Samples.combined,idents = "FP21_Cdh5_Sham")
saveRDS(Sham,file = paste0(sample,"pairwise/",sample,"_Sham_indiv_after_pairwise_integration.rds"))
TAC = subset(Samples.combined,idents = "FP20_Cdh5_TAC")
saveRDS(TAC,file = paste0(sample,"pairwise/",sample,"_TAC_indiv_after_pairwise_integration.rds"))

#---vln of shared DE genes-------------
Idents(Samples.combined)
print(VlnPlot(Samples.combined,group.by = "sub_cell_type_short",features =c("Col4a1","Col4a2","Col15a1") ,stack = T,flip = T,split.by = "stim",pt.size = 0,cols = c('#00C5C0','#FA756C')) +
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
        ggtitle(label = paste0("Colagens upregulated in TAC"))+
        theme(axis.title.x = element_blank()) + NoLegend())
ggsave(filename = paste0(sample,"Colagens upregulated in TAC_vln.jpeg"), width=10 , height = 10)

print(VlnPlot(Samples.combined,group.by = "sub_cell_type_short",features =c("Ltbp1","Ltbp4") ,stack = T,flip = T,split.by = "stim",pt.size = 0,cols = c('#00C5C0','#FA756C')) +
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
        ggtitle(label = paste0("Colagens upregulated in TAC"))+
        theme(axis.title.x = element_blank()) + NoLegend())
ggsave(filename = paste0(sample,"Ltbp's upregulated in TAC_vln.jpeg"), width=10 , height = 10)

print(VlnPlot(Samples.combined,group.by = "sub_cell_type_short",features =c("Nes","Cd34","Sox18","Aplnr") ,stack = T,flip = T,split.by = "stim",pt.size = 0,cols = c('#00C5C0','#FA756C')) +
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
        ggtitle(label = paste0("Colagens upregulated in TAC"))+
        theme(axis.title.x = element_blank()) + NoLegend())
ggsave(filename = paste0(sample,"interessting shared upregulated in TAC_vln.jpeg"), width=10 , height = 10)

#---annotation for EC clusters
#c("capillary EC","capillary artery EC","capillary vein EC","stressed EC","angiogenic EC","interferon EC","artery EC","lymphatic EC","DNA replicating EC","cycling EC","endocardial EC")
# ------------------------------------------------------------------------------------------- #
# add sub_cell_type 
# ------------------------------------------------------------------------------------------- #
celltypes<-c("Cap","CapA","CapV","StrEC","Ang","IntEC","Art","Lym","Rep","Cycl","Endo")
for (i in c(1:length(Samples.combined@meta.data$seurat_clusters))) {
  Samples.combined@meta.data$sub_cell_type[i]<-switch(as.numeric(Samples.combined@meta.data$seurat_clusters[i]), 
                                                      "capillary EC","capillary artery EC","capillary vein EC","stressed EC","angiogenic EC","interferon EC","artery EC","lymphatic EC","DNA replicating EC","cycling EC","endocardial EC")
  Samples.combined@meta.data$sub_cell_type_short[i]<-switch(as.numeric(Samples.combined@meta.data$seurat_clusters[i]), 
                                                            "Cap","CapA","CapV","StrEC","Ang","IntEC","Art","Lym","Rep","Cycl","Endo" )
}
Samples.combined$sub_cell_type_short = factor(Samples.combined$sub_cell_type_short,levels = celltypes)
