##------integration of fibroblasts from Pairwise filered col1a1, Gli1 and Pdgfrb------------------------
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
setwd("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/res0.3")
sample="Fibroblast_Col1a1_Gli1_Pdgfrb_integrated_"
dir.create("Fibroblast_Col1a1_Gli1_Pdgfrb_integrated")

#colors for genotypes
col.pdgfrb="#EC407A"
col.col1a1="#2196F3"
col.gli1="#4CAF50"
color.cl = c(col.col1a1,col.pdgfrb,col.gli1)

#load processed datafiles from pairwise
col1a1.sham = readRDS(file = "~/Documents/FP_scRNA/R stuff/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/RDS/Col1a1_integrated__Sham_indiv_after_pairwise_integration.rds")
col1a1.tac = readRDS(file = "~/Documents/FP_scRNA/R stuff/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/RDS/Col1a1_integrated__TAC_indiv_after_pairwise_integration.rds")

gli1.sham = readRDS(file = "~/Documents/FP_scRNA/R stuff/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/RDS/Gli1_integrated_Sham_indiv_after_pairwise_integration.rds")
gli1.tac = readRDS(file = "~/Documents/FP_scRNA/R stuff/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/RDS/Gli1_integrated_TAC_indiv_after_pairwise_integration.rds")

pdgfrb.sham = readRDS(file = "~/Documents/FP_scRNA/R stuff/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/RDS/PDGFRb_integrated__Sham_indiv_after_pairwise_integration_ONLY_FIBRO.rds")
pdgfrb.tac = readRDS(file = "~/Documents/FP_scRNA/R stuff/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/RDS/PDGFRb_integrated__TAC_indiv_after_pairwise_integration_ONLY_FIBRO.rds")

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
####################################################---integration--#################################################################
Samples.anchors <- FindIntegrationAnchors(object.list = list(col1a1.tac,col1a1.sham,gli1.tac,gli1.sham,pdgfrb.tac,pdgfrb.sham), dims = 1:20)

#--------------------exclude tdTom from the anchorlist-------------------------
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list

#only keep variable features after integration
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined,file = paste0(sample,"after_integration.rds"))

################################################----start of sample processing before analysis----##################################################
# subset the individual stim with just TAC or Sham
Samples.filtered=Samples.combined
TACs = Samples.filtered@meta.data$stim %in% c("FP12_Gli1_TAC","FP14_PDGFRb_TAC","FP18_Col1a1_TAC")
Samples.filtered@meta.data$stim[!TACs]="Sham"
Samples.filtered@meta.data$stim[TACs]="TAC" 
RidgePlot(Samples.filtered, assay = "RNA", features = "tdTomato-Green",group.by = "stim",log = T) + NoLegend()
Samples.combined=Samples.filtered
remove(Samples.filtered)

# add colum for only the genotype
Samples.filtered=Samples.combined
Gli1.ls = Samples.filtered@meta.data$orig.ident %in% c("FP12_Gli1_TAC","FP13_Gli1_Sham")
PDGFRb.ls = Samples.filtered@meta.data$orig.ident %in% c("FP14_PDGFRb_TAC","FP3_PDGFRb_Sham")
Col1a1.ls = Samples.filtered@meta.data$orig.ident %in% c("FP18_Col1a1_TAC","FP19_Col1a1_Sham")
Samples.filtered@meta.data$orig.geno[Gli1.ls]="Gli1"
Samples.filtered@meta.data$orig.geno[PDGFRb.ls]="Pdgfrb"
Samples.filtered@meta.data$orig.geno[Col1a1.ls]="Col1a1"
Samples.combined=Samples.filtered
remove(Samples.filtered)

#settle the order
Samples.combined$orig.ident = factor(Samples.combined$orig.ident,levels = c("FP19_Col1a1_Sham","FP18_Col1a1_TAC","FP3_PDGFRb_Sham","FP14_PDGFRb_TAC","FP13_Gli1_Sham","FP12_Gli1_TAC"))
Samples.combined$orig.geno = factor(Samples.combined$orig.geno, levels = c("Col1a1","Pdgfrb","Gli1"))

##-------first clustering after integration-------------------------------------------
Samples.combined = recluster(Samples.combined,0.5)

#plot for quality check
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 2, cols = c('#FA756C80','#00C5C080')) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"quality_check_clusters_res0.5.jpeg"), width=20 , height = 10)

#-------------ploting only tdtom
#quality check on tdTom expression
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(sample,"_tdTom.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(sample,"_tdTom_Vln.jpeg"), width=10 , height =10,limitsize = F)
RidgePlot(Samples.combined,features = "tdTomato-Green",group.by = "orig.ident",slot = "counts",assay = "RNA")
ggsave(filename = paste0(sample,"_tdTom_reads_filtered.jpeg"), width=10 , height = 10)

#run Genesorter to get GeneProb of tdtom per cluster
sg = sortGenes(Samples.combined@assays$RNA@data, Idents(Samples.combined))
colnames(sg$condGeneProb) = paste0(levels(Idents(Samples.combined)))
Samples.combined_geneProb = as.data.frame(sg$condGeneProb)["tdTomato-Green",]
Samples.combined_geneProb
remove(sg)

#------------basic nfeature and ncount plots------------------
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount_vln.jpeg"), width=10 , height =5)

#-------normalize and scale RNA data------------------------
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

source(file = "~/Documents/FP_scRNA/R stuff/scripts/methods/batch_test_umap.R")
batchTestUmap(Samples.combined,sample,integration = "seurat")

#--------save filtered object-----------
saveRDS(Samples.combined,file = paste0(sample,"_filtered_processed_final.rds"))
################################################----end of sample processing before analysis----##################################################

########################################################----start of sample analysis----###########################################################
#----------load filtered object---------
Samples.combined <- readRDS(file = "~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/RDS/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated__filtered_processed_final.rds")

  #clustertree for optimazation of cluster resolution
sc.int=Samples.combined
DefaultAssay(sc.int)<-"integrated"
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)

#according to clustree res0.3 might be more resonable
DefaultAssay(Samples.combined)<-"integrated"
Samples.combined<-FindClusters(Samples.combined, resolution = 0.3)
Samples.combined <- RunUMAP(Samples.combined, reduction = "pca", dims = 1:30,local.connectivity = 20) #optimized UMAP and replotet quality plots
DefaultAssay(Samples.combined)<-"RNA"

#cluster plots - seurat cluster
DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_only_clusters_res0.3.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_only_clusters_res0.3.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", pt.size = .5) + NoLegend() +  xlim(-6,9) + ylim(-7,6)
ggsave(filename = paste0(sample,"_only_clusters_res0.3_shifted.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_only_clusters_res0.3_shifted.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", label = F, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_only_clusters_res0.3_nolab.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_only_clusters_res0.3_nolab.svg"), width=10 , height = 10)
#cluster plots - sham vs. tac
DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = c('#FA756C80','#00C5C080')) + theme(legend.position = c(0.8,0.1))
ggsave(filename = paste0(sample,"_shamvstac_clusters_res0.3.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_shamvstac_clusters_res0.3.svg"), width=10 , height = 10)
#cluster plots - genotype
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno",group.by = "orig.geno", pt.size = 1,cols = color.cl,order = c("Gli1","Pdgfrb","Col1a1")) + NoLegend()
ggsave(filename = paste0(sample,"_genotype_clusters_res0.3_split.jpeg"), width=30 , height = 10)
ggsave(filename = paste0(sample,"_genotype_clusters_res0.3_split.svg"), width=30 , height = 10)
DimPlot(Samples.combined, reduction = "umap",group.by = "orig.geno", pt.size = 0.5,cols = c("#2196F380","#EC407A80","#4CAF5080"),order = c("Gli1","Pdgfrb","Col1a1")) + theme(legend.position = c(0.8,0.1))
ggsave(filename = paste0(sample,"_genotype_clusters_res0.3.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_genotype_clusters_res0.3.svg"), width=10 , height = 10)

#--------top10 & top5 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#result of findallmarkers
saveRDS(all.markers,file = "all.markers_markergenes_per_cluster_result.RDS")
#load results if they were saved allready before
all.markers <-readRDS(file = "~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/res0.3/all.markers_markergenes_per_cluster_result.RDS")
#get top10 markergenes per cluster
top10 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene))
#vertikal
DotPlot(Samples.combined, features = top10,dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"_top10_markers_per_cluster_v.jpeg"), width=5.3 , height = 20)
ggsave(filename = paste0(sample,"_top10_markers_per_cluster_v.svg"), width=5.3 , height = 20)
#horizontal
Idents(Samples.combined)= factor(Samples.combined@meta.data$seurat_clusters,levels(Idents(Samples.combined))[as.numeric(rev(levels(Idents(Samples.combined))))+1])
DotPlot(Samples.combined, features = top10,dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank())
ggsave(filename = paste0(sample,"_top10_featureDOTplot_custom_cl_h.jpeg"), width=22 , height = 3.2)
ggsave(filename = paste0(sample,"_top10_featureDOTplot_custom_cl_h.svg"), width=22 , height = 3.2)
Idents(Samples.combined)=Samples.combined@meta.data$seurat_clusters
#get top5 markergenes per cluster
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
top5 <- c(top5,"Dkk3") # added Dkk3 for cl6 to get 5 specific markers
#vertikal
DotPlot(Samples.combined, features = top5,dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.jpeg"), width=5.3 , height = 11)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.svg"), width=5.3 , height = 11)
#horizontal
Samples.combined$sub_cell_type_short = factor(Samples.combined$sub_cell_type_short,levels = rev(celltypes_short))
DotPlot(Samples.combined,group.by = "sub_cell_type_short", features = top5,dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title = element_blank())
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.jpeg"), width=12 , height = 3)
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.svg"), width=12 , height = 3)


#-----bar graph per genotype
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$orig.geno) 
total.col1a1 = sum(breakdown[,1])
total.gli1 = sum(breakdown[,2])
total.pdgfrb = sum(breakdown[,3])
ratio.col1a1.gli1 = total.col1a1/total.gli1
ratio.col1a1.pdgfrb = total.col1a1/total.pdgfrb
breakdown[,2] = round(breakdown[,2]*ratio.col1a1.gli1)
breakdown[,3] = round(breakdown[,3]*ratio.col1a1.pdgfrb)
breakdown=t(breakdown)
breakdown <- round(apply(breakdown, 2, function(x){x*100/sum(x)}),2)
breakdown.df = as.data.frame(breakdown)
breakdown.df.new = breakdown.df[1,]
breakdown.df.new[2,] = breakdown.df[3,]
breakdown.df.new[3,] = breakdown.df[2,]
breakdown.df = breakdown.df.new
breakdown.df = melt(t(breakdown.df))
ggplot(data=breakdown.df, aes(x=Var1, y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = levels(factor(breakdown.df$Var1)),name="Cluster")+
  scale_y_continuous(name="Normalizes percentage of contribution [%]") +
  scale_fill_manual( values = c("Col1a1" = col.col1a1, "Gli1" = col.gli1, "Pdgfrb" = col.pdgfrb))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.y = element_text(size=20,face = "bold"),axis.title.y = element_text(size = 15),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15),legend.title = element_blank(),legend.position = c(0.92,0.92),legend.key.size = unit(0.5, "cm"),legend.box.background = element_rect(colour = "black"))
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster_normalized.jpeg"), width=10 , height =10)
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster_normalized.svg"), width=10 , height =10)

#------advaced bar garph----Sham vs. TAC--------------------------------------------
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$stim) 
breakdown[,1]= round(breakdown[,1]*(sum(breakdown[,2])/sum(breakdown[,1]))) #normalize cellnumber of TAC to Sham
breakdown=t(breakdown)
breakdown <- apply(breakdown, 2, function(x){x*100/sum(x)})
breakdown = round(as.data.frame(breakdown),digits = 2)
breakdown.df = as.data.frame(breakdown)
breakdown.df = melt(t(breakdown.df))

ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=rev(Var2))) +
  geom_bar(stat="identity")+
  geom_text(aes(label=c(celltypes_short,celltypes_short),y=1),hjust = 0, color="black", size=10)+
  coord_flip() + 
  scale_y_continuous(name="Percentage of total cells [%]") +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  NoLegend()+
  geom_hline(yintercept=50, linetype="dashed", color = "black")
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.jpeg"), width=5 , height = 10)
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.svg"), width=5 , height = 10)

###------sample specific plot for GOs-------------------------------
#selected two most interessting and representative GOs per cluster
go.picked = c("GO:0031012","GO:0005581","GO:0010942","GO:0043620","GO:0031589","GO:0050840","GO:0071953","GO:0009611","GO:0030198","GO:0071363","GO:0045087","GO:0035456","GO:0043292","GO:0001968")

#run GOs on top20 marker genes from Findallmarkers before ploting
go.results.up.top2 <- go.results.up[go.results.up$term_id %in% go.picked,]
order_desc_up <- unique(go.results.up.top2$term_name)

#-plot up vertikal
ggplot(go.results.up.top2, aes(x=factor(term_name, level = order_desc_up), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
  geom_point() + 
  scale_size(range = c(2,6)) +
  coord_flip() +
  ylab("Cluster") + xlab("") + ggtitle("Top2 repesentative GOs per cluster")+
  theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5),legend.key.size = unit(0.3, "cm"))
ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top2 representative GOs per cluster.jpeg"), width=6, height = 3)
ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top2 representative GOs per cluster.svg"), width=6, height = 3)

# ------------------------------------------------------------------------------------------- #
# add sub_cell_type according to annotation 
# ------------------------------------------------------------------------------------------- #
#for ------custom clusters of fibroblasts samples -----------------------------
#number of clusters: 6
#cl0 "Fibro-1",cl1 "Fibro-stressed",cl2 "Fibro-2",cl3 "Fibro-3",cl4 "Myofibroblast", cl5 "Fibro-Interferon", cl6 "Fibro-4"
celltypes_short<-c("Fib1", "StrFib","Fib2", "Fib3", "MyoF","IntF","Fib4")
sub_cell_type <- c("Fibro-1", "Fibro-stressed","Fibro-2", "Fibro-3", "Myofibroblast","Fibro-Interferon","Fibro-4")
for (i in c(1:length(Samples.combined@meta.data$seurat_clusters))) {
  Samples.combined@meta.data$sub_cell_type[i]<-switch(as.numeric(Samples.combined@meta.data$seurat_clusters[i]), 
                                                      "Fibro-1", "Fibro-stressed","Fibro-2", "Fibro-3", "Myofibroblast","Fibro-Interferon","Fibro-4")
  Samples.combined@meta.data$sub_cell_type_short[i]<-switch(as.numeric(Samples.combined@meta.data$seurat_clusters[i]), 
                                                            "Fib1", "StrFib","Fib2", "Fib3", "MyoF","IntF","Fib4")
}
Samples.combined$sub_cell_type = factor(Samples.combined$sub_cell_type,levels = c("Fibro-1", "Fibro-stressed","Fibro-2", "Fibro-3", "Myofibroblast","Fibro-Interferon","Fibro-4"))
Samples.combined$sub_cell_type_short = factor(Samples.combined$sub_cell_type_short,levels = c("Fib1", "StrFib","Fib2", "Fib3", "MyoF","IntF","Fib4"))

# ------------------------------------------------------------------------------------------- #
# special volcano plot for DE genes of cl4 and 5, run after "gProfileR2 for GO on DE" script
# ------------------------------------------------------------------------------------------- #
#adjusting cutoffs for plotting
cluster = 4 #ones run with 4 and 5
library(ggrepel)
DE.result.cl = stim.response[stim.response$cluster==cluster,]
DE.result.cl$p_val_adj[DE.result.cl$p_val_adj==0] = 1.000000e-300 #genes with significe of 0 are set to 1.000000e-300
DE.result.cl$p_val_adj[DE.result.cl$p_val_adj<1.000000e-150] = 1.000000e-150 #genes with very high pval_adjust were set to 1.0e-150 for visualisation
#look for genes about threshhold (cutoff for avg_logFC is deactivated)(cutoff for the p_val_adj ist dynamically changed for each sample to the man p_val_adj of all samples)
for (gene in DE.result.cl$gene_name) {
  #if (DE.result.cl$avg_logFC[DE.result.cl$gene_name == gene]>0.5 | DE.result.cl$avg_logFC[DE.result.cl$gene_name == gene]< (-0.5) ) {
  if (-log10(DE.result.cl$p_val_adj[DE.result.cl$gene_name == gene])>(mean(-log10(DE.result.cl$p_val_adj))) & DE.result.cl$unique[DE.result.cl$gene_name == gene] !="shared") { 
    DE.result.cl$high[DE.result.cl$gene_name == gene] = gene 
    DE.result.cl$overth[DE.result.cl$gene_name == gene] = TRUE
  } else {DE.result.cl$overth[DE.result.cl$gene_name == gene] = FALSE
  DE.result.cl$high[DE.result.cl$gene_name == gene] = ""}
  #} else {
  #  DE.result.cl$high[DE.result.cl$gene_name == gene] = ""
  #  DE.result.cl$overth[DE.result.cl$gene_name == gene] = FALSE
  #}
}
DE.result.cl <- arrange(DE.result.cl,p_val_adj)

top5.up.and.down = c(head(DE.result.cl$gene_name[DE.result.cl$avg_logFC>0],n=5),head(DE.result.cl$gene_name[DE.result.cl$avg_logFC<0],n=5))
for (gene in top5.up.and.down) {
  DE.result.cl$high[DE.result.cl$gene_name == gene] = gene 
  DE.result.cl$overth[DE.result.cl$gene_name == gene] = TRUE
}

#add if alo unique or shared
for (gene in DE.result.cl$gene_name) {
  if (DE.result.cl$overth[DE.result.cl$gene_name == gene] == FALSE) {
    DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "shared"
    if (DE.result.cl$unique[DE.result.cl$gene_name == gene] != "shared") {DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "unique"}
  } else {
    if (DE.result.cl$unique[DE.result.cl$gene_name == gene] == "shared") {
      DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "shared"
    } else {
      DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "unique"
    }
  }
}
#volcanoplot
ggplot(DE.result.cl) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj),colour=high.unique),size=2.5) +
  ggtitle(paste0("Cluster ",cluster, " overexpression in ",treatment_sample, "top5 up&down + unique")) +
  xlab("avg_logFC") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept=(mean(-log10(DE.result.cl$p_val_adj))+1), linetype="dashed", color = "grey")+
  geom_vline(xintercept=0.5, linetype="dashed", color = "grey")+
  geom_vline(xintercept=-0.5, linetype="dashed", color = "grey")+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = high),size = 5,min.segment.length = 0,segment.size = 1) +
  scale_color_manual(values = c("shared"="blue","unique"="red")) +
  theme(legend.title = element_blank(),legend.position = c(0.95,0.1),plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))  
ggsave(filename = paste0("gprofiler2_out_put/volcano_plots_of_DE/",sample,"volcano_plot_cluster_custom",cluster,".jpeg"), width=8, height =5)
ggsave(filename = paste0("gprofiler2_out_put/volcano_plots_of_DE/svg/",sample,"volcano_plot_cluster_custom",cluster,".svg"), width=8, height =5)

#------------------------calculate ratios of ECM between conditions (per cluster) and add to violin----------------------------------------------
#core matrisome
Samples.combined$Core_matrisome1_delta = Samples.combined$Core_matrisome1
for (cluster.nr in levels(factor(Samples.combined$seurat_clusters))) {
  median.sham = median(Samples.combined$Core_matrisome1[Samples.combined$seurat_clusters==cluster.nr&Samples.combined$stim=="Sham"])
  median.TAC = median(Samples.combined$Core_matrisome1[Samples.combined$seurat_clusters==cluster.nr&Samples.combined$stim=="TAC"])
  Samples.combined$Core_matrisome1_delta[Samples.combined$seurat_clusters==cluster.nr] = round(median.TAC-median.sham,digits = 3)
}

VlnPlot(Samples.combined, features = "Core_matrisome1", pt.size = 0,split.by = "stim",split.plot= T,cols = c('#00C5C0','#FA756C')) +
  geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
  geom_text(aes(label=Samples.combined$Core_matrisome1_delta, y=c(rep(0.45,times=length(Samples.combined$orig.ident)))),hjust = -0.3, color="black", size=4)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0("ECM scoring output/",sample,"core_matrisome_vln_delta.jpeg"), width=10 , height = 10)
ggsave(filename = paste0("ECM scoring output/",sample,"core_matrisome_vln_delta.svg"), width=10 , height = 10)

#collagene
Samples.combined$Core_Collagens1_delta = Samples.combined$Collagens1
for (cluster.nr in levels(factor(Samples.combined$seurat_clusters))) {
  median.sham = median(Samples.combined$Collagens1[Samples.combined$seurat_clusters==cluster.nr&Samples.combined$stim=="Sham"])
  median.TAC = median(Samples.combined$Collagens1[Samples.combined$seurat_clusters==cluster.nr&Samples.combined$stim=="TAC"])
  Samples.combined$Core_Collagens1_delta[Samples.combined$seurat_clusters==cluster.nr] = round(median.TAC-median.sham,digits = 3)
}

VlnPlot(Samples.combined, features = "Collagens1", pt.size = 0,split.by = "stim",split.plot= T,cols = c('#00C5C0','#FA756C')) +
  geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
  geom_text(aes(label=Samples.combined$Core_Collagens1_delta, y=c(rep(0.5,times=length(Samples.combined$orig.ident)))),hjust = -0.3, color="black", size=4)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0("ECM scoring output/",sample,"Collagens1_vln_delta.jpeg"), width=10 , height = 10)
ggsave(filename = paste0("ECM scoring output/",sample,"Collagens1_vln_delta.svg"), width=10 , height = 10)

##---------reduction of progeny result------------------------------
progeny_hmap_v = pheatmap(t(summarized_progeny_scores_df[c(4,9),]),fontsize=14,
                        fontsize_row = 10,
                        color=myColor, breaks = progenyBreaks,
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA,cluster_cols = F,cluster_rows = F,filename = paste0(sample,"_TF_PROGENy_heatmap_reduced_v.jpeg"),width = 2,height = 6)

progeny_hmap_h = pheatmap(summarized_progeny_scores_df[c(4,9),],fontsize=14,
                          fontsize_row = 10,
                          color=myColor, breaks = progenyBreaks,
                          main = "PROGENy (500)", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,cluster_cols = F,cluster_rows = F,filename = paste0(sample,"_TF_PROGENy_heatmap_reduced_h.jpeg"),width = 8,height = 2)

# --------------check results for DE saved in gene.counter.df against the genes of Jak-stat / tgf-beta signaling from progeny
library(progeny)
mouse.prog.df <- progeny::getModel("Mouse")
genes.of.int.shared <- as.character(gene.counter.df$genes[gene.counter.df$count>1])
genes.of.int.unique <- as.character(gene.counter.df$genes[gene.counter.df$count==1])
jak.stat.genes <- rownames(mouse.prog.df[mouse.prog.df$`JAK-STAT`>0,])
tgfb.genes <- rownames(mouse.prog.df[mouse.prog.df$TGFb>0,])

#genes that are shared and also in the jak.stat dataset
intersect(genes.of.int.shared,jak.stat.genes)
gene.counter.df[gene.counter.df$genes %in% intersect(genes.of.int.shared,jak.stat.genes),]
#genes that are unique and also in the jak.stat dataset -> 0 -> gobal effect
intersect(genes.of.int.unique,jak.stat.genes)
DE.heatmap[DE.heatmap$gene_name %in% intersect(genes.of.int.unique,jak.stat.genes),]

#genes that are shared and also in the tgfb dataset
intersect(genes.of.int.shared,tgfb.genes)
gene.counter.df[gene.counter.df$genes %in% intersect(genes.of.int.shared,tgfb.genes),]
#genes that are unique and also in the tgfb dataset
intersect(genes.of.int.unique,tgfb.genes)
DE.heatmap[DE.heatmap$gene_name  %in% intersect(genes.of.int.unique,tgfb.genes),]
