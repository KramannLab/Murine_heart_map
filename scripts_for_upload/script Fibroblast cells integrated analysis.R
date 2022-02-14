#script for analysis of Fibroblast Sham and TAC integrated
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
library(pheatmap)
library(scProportionTest)
source("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/scripts_for_upload/helper_functions_DEG_GO_ploting.R")


setwd("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_integrated_all/")
sample="Fibroblast_integrated_"

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
col.col1a1="#228833"
col.gli1="#4477AA"
geno.cl = c(col.col1a1,col.pdgfrb,col.gli1)
geno.cl.light = c("#22883380","#EE667780","#4477AA80")
#for subclustering #RColorBrewer
display.brewer.pal(name = "Paired",n = 12)
cluster.col=brewer.pal(name = "Paired",n=11)
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)

################################################----start of sample processing before analysis----##################################################
Samples.combined = readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering_Sham&TAC_integrated/Fibroblast_datasets_integrated_after_bacis_filters+TAC28d.rds")
#add stim based on sha, tac or tac28d
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP12_Gli1_TAC','TAC',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP14_PDGFRb_TAC','TAC',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP18_Col1a1_TAC','TAC',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP13_Gli1_Sham','Sham',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP3_PDGFRb_Sham','Sham',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP19_Col1a1_Sham','Sham',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP71_Gli1_TAC_28d','TAC_28',Samples.combined@meta.data$stim)
Samples.combined$stim <- factor(Samples.combined$stim,levels=c("Sham","TAC","TAC_28"))
#add genotype to metadata
# add colum for only the genotype
Samples.combined$orig.geno=Samples.combined$orig.ident
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("FP12_Gli1_TAC","FP13_Gli1_Sham","FP71_Gli1_TAC_28d")]="Gli1"
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("FP14_PDGFRb_TAC","FP3_PDGFRb_Sham")]="Pdgfrb"
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("FP18_Col1a1_TAC","FP19_Col1a1_Sham")]="Col1a1"
Samples.combined$orig.geno=factor(Samples.combined$orig.geno,levels = c("Col1a1","Pdgfrb","Gli1"))
#first clustering
Samples.combined = recluster(Samples.combined,0.5)

#plot before tdTom filtering steps
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = stim.col.light) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"before_processing_clusters_res0.5.jpeg"), width=20 , height = 10)
#cluster plots - genotype
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno",group.by = "orig.geno", pt.size = 1,cols = geno.cl) + NoLegend()
ggsave(filename = paste0(sample,"_genotype_split.jpeg"), width=30 , height = 10)
#cluster plots - sample
#add celltypes + stim with levels
lvl=NULL
for (levels in levels(Samples.combined$orig.geno)) {lvl=c(lvl,paste0(levels,"_",levels(Samples.combined$stim)))}
Samples.combined$celltypes.stim=factor(paste0(Samples.combined$orig.geno,"_",Samples.combined$stim),levels = lvl)
DimPlot(Samples.combined,split.by = "celltypes.stim", reduction = "umap", pt.size = 0.1,ncol = 2) + NoLegend()
ggsave(filename = paste0(sample,"_sample_split.jpeg"), width=10 , height = 20)

#-------------ploting only tdtom
#quality check on tdTom expression
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(sample,"tdTom_unfiltered.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(sample,"tdTom_Vln_unfiltered.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)+theme(axis.text.y = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.title.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(sample,"tdTom_Vln_unfiltered_v2.jpeg"), width=10 , height =10)
RidgePlot(Samples.combined,features = "tdTomato-Green",group.by = "orig.ident",slot = "counts",assay = "RNA")
ggsave(filename = paste0(sample,"tdTom_reads_unfiltered.jpeg"), width=10 , height = 10)

#run Genesorter to get GeneProb of tdtom per cluster
sg = sortGenes(Samples.combined@assays$RNA@data, Idents(Samples.combined))
colnames(sg$condGeneProb) = paste0(levels(Idents(Samples.combined)))
Samples.combined_geneProb = as.data.frame(sg$condGeneProb)["tdTomato-Green",]
Samples.combined_geneProb
remove(sg)

#------------basic nfeature and ncount plots
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount_unfiltered.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount__unfilteredvln.jpeg"), width=10 , height =5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.text.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(sample,"_nfeature_unfilteredvln_v2.jpeg"), width=10 , height =10)

#--------top10 markers per cluster
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#get top10 markergenes per cluster
top10 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene))
#vertikal
DotPlot(Samples.combined, features = top10,dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"_top10_markers_unfiltered.jpeg"), width=7 , height = 25)

#----------subsetting to filter out tdtom low and nfeature low clusters, also exclude mural cells
Samples.combined=subset(Samples.combined,idents = c("6","8","9","12","13","14","15"),invert=T)
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined=subset(Samples.combined,cells = WhichCells(Samples.combined,slot = "counts",expression = `tdTomato-Green` == 0),invert=T)
Samples.combined = recluster(Samples.combined,0.2,rep.str = 1,loc.con = 20)
#------cluster6 was highly similar to cluster1 in further analysis and therefore merged manualy
types <- list("6"="1")
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$seurat_clusters <- Idents(Samples.combined)
Samples.combined$seurat_clusters = factor(Samples.combined$seurat_clusters,levels = c(0:5))
Idents(Samples.combined)="seurat_clusters"

#plot after tdTom filtering step and before nfeature filter
p1 <- DimPlot(Samples.combined, reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = stim.col.light) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"_filtered_clusters_res0.3.jpeg"), width=20 , height = 10)

#------------basic nfeature ncount plots and tdtom
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount_vln.jpeg"), width=10 , height =5)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(sample,"tdTom_Vln.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0,cols=cluster.col)+theme(axis.text.y = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.title.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(sample,"tdTom_Vln_v2.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0,cols=cluster.col)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,size = 20),axis.title.x = element_blank(),axis.text.y = element_text(size = 25))+NoLegend()
ggsave(filename = paste0(sample,"_nfeature_v2.jpeg"), width=10 , height =10)

#plot after processing  
DimPlot(Samples.combined, reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"after_processing_clusters_res0.2.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.2.svg"), width=10 , height = 10)
#no label version
DimPlot(Samples.combined, reduction = "umap",cols = cluster.col,  pt.size = 1) + NoLegend()
ggsave(filename = paste0(sample,"after_processing_clusters_res0.2_noLabel.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.2_noLabel.svg"), width=10 , height = 10)
#cluster plots - genotype
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno",group.by = "orig.geno", pt.size = 1,cols = geno.cl) + NoLegend()
ggsave(filename = paste0(sample,"_genotype_clusters_res0.2_split.jpeg"), width=30 , height = 10)
ggsave(filename = paste0(sample,"_genotype_clusters_res0.2_split.svg"), width=30 , height = 10)

#-------normalize and scale RNA data
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

#--------save filtered object
saveRDS(Samples.combined,file = paste0(sample,"_filtered_processed.rds"))
################################################----end of sample processing before analysis----##################################################
########################################################----start of sample analysis----###########################################################
#----------load filtered object---------
Samples.combined <- readRDS(file ="~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_integrated_all/Fibroblast_integrated__filtered_processed.rds")

#-------clustertree for optimazation of cluster resolution-----
sc.int=Samples.combined
DefaultAssay(sc.int)<-"integrated"
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)

#-------annotation for Fib clusters according to marker genes, GO terms and ECM scoring-------
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Fibroblast 1","1"="Stressed Fibroblast","2"="Fibroblast 2","3"="Fibroblast 3","4"="ECM-Fibroblast","5"="Interferon Fibroblast")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes <- Idents(Samples.combined)
cell.levels = c("Fibroblast 1","Stressed Fibroblast","Fibroblast 2","Fibroblast 3","ECM-Fibroblast","Interferon Fibroblast")
Samples.combined$celltypes = factor(Samples.combined$celltypes,levels = cell.levels)
#short abroviation of annotation
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Fib1","1"="StrFib","2"="Fib2","3"="Fib3","4"="ECM-Fib","5"="IntFib")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes.short <- Idents(Samples.combined)
cell.levels = c("Fib1","StrFib","Fib2","Fib3","ECM-Fib","IntFib")
Samples.combined$celltypes.short = factor(Samples.combined$celltypes.short,levels = cell.levels)
Idents(Samples.combined)="seurat_clusters"
#add celltypes + stim with levels
lvl=NULL
for (levels in levels(Samples.combined$celltypes)) {lvl=c(lvl,paste0(levels,"_",levels(Samples.combined$stim)))}
Samples.combined$celltypes.stim=factor(paste0(Samples.combined$celltypes,"_",Samples.combined$stim),levels = lvl)
#save after renaming
saveRDS(Samples.combined,file = paste0(sample,"_filtered_processed.rds"))
#dimplot short labels
DimPlot(Samples.combined,group.by = "celltypes.short", reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"after_processing_clusters_res0.2_labeled_short.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.2_labeled_short.svg"), width=10 , height = 10)
#dimplot stim split
DimPlot(Samples.combined,group.by = "stim",split.by="stim", reduction = "umap",cols = stim.col, pt.size = 1) + NoLegend()
ggsave(filename = paste0(sample,"split_stim.jpeg"), width=30 , height = 10)


#-------top10 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#get top5 markergenes per cluster
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
#vertikal
DotPlot(Samples.combined, group.by = "celltypes.short",features = top5,dot.scale = 10,assay = "RNA")+ coord_flip()+theme(axis.text.x = element_text(angle = 45,hjust=1),axis.text.y = element_text(face = "italic")) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.jpeg"), width=5.2, height = 9.5)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.svg"), width=5.2 , height = 9.5)

#---expression of specific genes
FeaturePlot(Samples.combined,features = "Thbs4",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Thbs4.jpeg"), width=10 , height = 10)
VlnPlot(Samples.combined,group.by = "celltypes.short",features = "Thbs4",split.by = "stim",cols = stim.col,pt.size = 0) + theme(axis.title.x = element_blank())
ggsave(filename = paste("multiple_marker_plots/",sample,"_Thbs4_Vln.jpeg"), width=10 , height = 4)

FeaturePlot(Samples.combined,features = "Acta2",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Acta2.jpeg"), width=10 , height = 10)

FeaturePlot(Samples.combined,features = "Wif1",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Wif1.jpeg"), width=10 , height = 10)

FeaturePlot(Samples.combined,features = "Dkk3",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Dkk3.jpeg"), width=10 , height = 10)

FeaturePlot(Samples.combined,features = "Mical2",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Mical2.jpeg"), width=10 , height = 10)
VlnPlot(Samples.combined,features = "Mical2",split.by = "stim",cols = stim.col,pt.size = 0)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Mical2_Vln_split.jpeg"), width=10 , height = 10)

FeaturePlot(Samples.combined,features = "Palld",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Palld.jpeg"), width=10 , height = 10)
VlnPlot(Samples.combined,features = "Palld",split.by = "stim",cols = stim.col,pt.size = 0)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Palld_Vln_split.jpeg"), width=10 , height = 10)

#----barplot of genotype contribution per cluster (normalized to even input per genotype)-----
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
breakdown.df$Var2 = factor(breakdown.df$Var2,levels = levels(Samples.combined$orig.geno))
ggplot(data=breakdown.df, aes(x=Var1, y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = levels(factor(breakdown.df$Var1)),name="Cluster")+
  scale_y_continuous(name="Normalizes percentage of contribution [%]") +
  scale_fill_manual( values = c("Col1a1" = col.col1a1,  "Pdgfrb" = col.pdgfrb,"Gli1" = col.gli1))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.y = element_text(size=20,face = "bold"),axis.title.y = element_text(size = 15),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15),legend.title = element_blank(),legend.position = c(0.92,0.92),legend.key.size = unit(0.5, "cm"),legend.box.background = element_rect(colour = "black"))
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster_normalized.jpeg"), width=10 , height =10)
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster_normalized.svg"), width=10 , height =10)

#-------barplot of condition contribution per cluster (input normalized)-----
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$stim) 
breakdown[,1]=100*(breakdown[,1]/sum(breakdown[,1]))
breakdown[,2]=100*(breakdown[,2]/sum(breakdown[,2]))
breakdown[,3]=100*(breakdown[,3]/sum(breakdown[,3]))
breakdown=t(breakdown)
breakdown <- apply(breakdown, 2, function(x){x*100/sum(x)})
breakdown = round(as.data.frame(breakdown),digits = 2)
breakdown.df = as.data.frame(breakdown)
breakdown.df = melt(t(breakdown.df))
breakdown.df$Var2=factor(breakdown.df$Var2,levels = c("TAC_28","TAC","Sham"))
celltypes.short = c("Fib1","StrFib","Fib2","Fib3","ECM-Fib","IntFib")

ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=c(celltypes.short,celltypes.short,celltypes.short),y=1),hjust = 0, color="black", size=10)+
  coord_flip() + 
  scale_y_continuous(name="Percentage of total cells [%]") +
  scale_fill_manual(values = rev(stim.col))+
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_hline(yintercept=50, linetype="dashed", color = "black")
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.jpeg"), width=5 , height = 10)
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.svg"), width=5 , height = 10)

#-------proportion analysis-------
prop_test <- sc_utils(Samples.combined)
## Sham vs. TAC
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltypes.short",
  sample_1 = "Sham", sample_2 = "TAC",
  sample_identity = "stim",
  n_permutations = 10000
)
permutation_plot(prop_test) + theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"_proportionTest_sham_vs_tac.jpeg"), width=5 , height = 5)
results.sham.vs.tac <-prop_test@results$permutation
results.sham.vs.tac$TAC_28 = NA
## Sham vs. TAC28
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltypes.short",
  sample_1 = "Sham", sample_2 = "TAC_28",
  sample_identity = "stim",
  n_permutations = 10000
)
permutation_plot(prop_test) + theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"_proportionTest_sham_vs_tac28.jpeg"), width=5 , height = 5)
results.sham.vs.tac28 <-prop_test@results$permutation
results.sham.vs.tac28$TAC = NA
## TAC vs. TAC28
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltypes.short",
  sample_1 = "TAC", sample_2 = "TAC_28",
  sample_identity = "stim",
  n_permutations = 10000
)
permutation_plot(prop_test) + theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"_proportionTest_tac_vs_tac28.jpeg"), width=5 , height = 5)
results.tac.vs.tac28 <-prop_test@results$permutation
results.tac.vs.tac28$Sham = NA
## save summary
results.all <- rbind(results.sham.vs.tac,results.sham.vs.tac28,results.tac.vs.tac28,fill=T)
write_xlsx(results.all,path = paste0(sample,"_proportionTest_results.xlsx"))

#-------get GO-terms of markergenes per cluster-------------
types <- list("0"="Fib1","1"="StrFib","2"="Fib2","3"="Fib3","4"="ECM-Fib","5"="IntFib")
go.of.marker.genes <- get.markergenes.GOs(Samples.combined,sample,create.folder = T,all.markers = all.markers)
plot.markergenes.GOs(sample,go.of.marker.genes)
plot.markergenes.GOs(sample,go.of.marker.genes,topGOs.number = 3,custom.h = 2.4,custom.w = 4.8)

#--------scoring for cell cycle--------------------------------
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
#perform scoring for each cell
Samples.combined <- CellCycleScoring(Samples.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#plot results
v1<-VlnPlot(Samples.combined,features = "S.Score",group.by = "celltypes.short",cols = stim.col,split.by = "stim",pt.size = 0) + theme(axis.text.x = element_text(hjust = 1,angle = 45),axis.title.x = element_blank())
v1
ggsave(filename = paste0(sample,"_S.score_vln.jpeg"), width=10 , height = 10)
v2<-VlnPlot(Samples.combined,features = "G2M.Score",group.by = "celltypes.short",cols = stim.col,split.by = "stim",pt.size = 0) + theme(axis.text.x = element_text(hjust = 1,angle = 45),axis.title.x = element_blank())
v2
ggsave(filename = paste0(sample,"_G2M.score_vln.jpeg"), width=10 , height = 10)
v1 <- v1 + theme(axis.text.x = element_blank())
ggarrange(v1,v2,ncol = 1)
ggsave(filename = paste0(sample,"_S&G2M_score_stack_vln.jpeg"), width=5 , height = 5)


#-------DEG testing-------------------------------------------------
Idents(Samples.combined)="stim"
#Sham vs TAC_14
sample="Fibroblast_integrated_shamvstac14_"
sham.vs.tac14 <- subset(Samples.combined,idents=c("Sham","TAC"))
sham.vs.tac14.de <- get.DEG.between.conditions.per.cluster(sham.vs.tac14,sample,control_sample = "Sham",treatment_sample = "TAC")
sham.vs.tac14.de <- check.for.unique.DE(sham.vs.tac14.de)

#Sham vs TAC_28
sample="Fibroblast_integrated_shamvstac28_"
sham.vs.tac28 <-subset(Samples.combined,idents=c("Sham","TAC_28"))
sham.vs.tac28.de <- get.DEG.between.conditions.per.cluster(sham.vs.tac28,sample,control_sample = "Sham",treatment_sample = "TAC_28")
sham.vs.tac28.de <- check.for.unique.DE(sham.vs.tac28.de)

#TAC_14 vs TAC_28
sample="Fibroblast_integrated_tac14vstac28_"
tac14.vs.tac28 <-subset(Samples.combined,idents=c("TAC","TAC_28"))
tac14.vs.tac28.de <- get.DEG.between.conditions.per.cluster(tac14.vs.tac28,sample,control_sample = "TAC",treatment_sample = "TAC_28")
tac14.vs.tac28.de <- check.for.unique.DE(tac14.vs.tac28.de)


#volcano plots for ECM-Fib of Sham vs. TAC & Sham vs. TAC28 as most interesting cluster
## sham vs tac14
sham.vs.tac14.de.cl4 <- sham.vs.tac14.de[sham.vs.tac14.de$cluster==4,]
sham.vs.tac14.de.cl4$high <-""
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,desc(avg_logFC))
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,desc(unique))
sham.vs.tac14.de.cl4$high[c(1:4)] <-sham.vs.tac14.de.cl4$gene_name[c(1:4)]
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,avg_logFC)
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,desc(unique))
sham.vs.tac14.de.cl4$high[c(1:5)] <-sham.vs.tac14.de.cl4$gene_name[c(1:5)]
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,p_val_adj)
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,desc(unique))
sham.vs.tac14.de.cl4$high[c(1:10)] <-sham.vs.tac14.de.cl4$gene_name[c(1:10)]
sham.vs.tac14.de.cl4$unique<-factor(sham.vs.tac14.de.cl4$unique,levels = c("shared","unique"))
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,unique)

ggplot(sham.vs.tac14.de.cl4) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj),colour=unique),size=2) +
  ggtitle(paste0("ECM-Fib DEG Sham vs TAC14")) +
  xlab("avg_logFC") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept=(mean(-log10(sham.vs.tac14.de.cl4$p_val_adj))+1), linetype="dashed", color = "grey")+
  geom_vline(xintercept=0.5, linetype="dashed", color = "grey")+
  geom_vline(xintercept=-0.5, linetype="dashed", color = "grey")+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = high),size = 7,min.segment.length = 0,segment.size = 1,force = 2,box.padding = unit(0.5, "cm")) +
  scale_color_manual(values = c("shared"="blue","unique"="red")) +
  theme(legend.title = element_blank(),legend.position = c(0.95,0.1),plot.title = element_text(size = 20, hjust = 0.5),axis.title = element_text(size = 20),axis.text = element_text(size=15))  +
  NoLegend()
ggsave(filename = paste0(sample,"volcano_plot_ECM-Fib DEG Sham vs TAC14.jpeg"), width=8, height =8)
ggsave(filename = paste0(sample,"volcano_plot_ECM-Fib DEG Sham vs TAC14.svg"), width=8, height =8)
## sham vs tac28
sham.vs.tac28.de.cl4 <- sham.vs.tac28.de[sham.vs.tac28.de$cluster==4,]
sham.vs.tac28.de.cl4$high <-""
sham.vs.tac28.de.cl4 <- arrange(sham.vs.tac28.de.cl4,desc(avg_logFC))
sham.vs.tac28.de.cl4 <- arrange(sham.vs.tac28.de.cl4,desc(unique))
sham.vs.tac28.de.cl4$high[c(1:4)] <-sham.vs.tac28.de.cl4$gene_name[c(1:4)]
sham.vs.tac28.de.cl4 <- arrange(sham.vs.tac28.de.cl4,avg_logFC)
sham.vs.tac28.de.cl4 <- arrange(sham.vs.tac28.de.cl4,desc(unique))
sham.vs.tac28.de.cl4$high[c(1:5)] <-sham.vs.tac28.de.cl4$gene_name[c(1:5)]
sham.vs.tac28.de.cl4 <- arrange(sham.vs.tac28.de.cl4,p_val_adj)
sham.vs.tac28.de.cl4 <- arrange(sham.vs.tac28.de.cl4,desc(unique))
sham.vs.tac28.de.cl4$high[c(1:10)] <-sham.vs.tac28.de.cl4$gene_name[c(1:10)]
sham.vs.tac28.de.cl4$unique<-factor(sham.vs.tac28.de.cl4$unique,levels = c("shared","unique"))
sham.vs.tac28.de.cl4 <- arrange(sham.vs.tac28.de.cl4,unique)

ggplot(sham.vs.tac28.de.cl4) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj),colour=unique),size=2) +
  ggtitle(paste0("ECM-Fib DEG Sham vs TAC28")) +
  xlab("avg_logFC") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept=(mean(-log10(sham.vs.tac28.de.cl4$p_val_adj))+1), linetype="dashed", color = "grey")+
  geom_vline(xintercept=0.5, linetype="dashed", color = "grey")+
  geom_vline(xintercept=-0.5, linetype="dashed", color = "grey")+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = high),size = 7,min.segment.length = 0,segment.size = 1,force = 2,box.padding = unit(0.5, "cm")) +
  scale_color_manual(values = c("shared"="blue","unique"="red")) +
  theme(legend.title = element_blank(),legend.position = c(0.95,0.1),plot.title = element_text(size = 20, hjust = 0.5),axis.title = element_text(size = 20),axis.text = element_text(size=15))  +
  NoLegend()
ggsave(filename = paste0(sample,"volcano_plot_ECM-Fib DEG Sham vs TAC28.jpeg"), width=8, height =8)
ggsave(filename = paste0(sample,"volcano_plot_ECM-Fib DEG Sham vs TAC28.svg"), width=8, height =8)

#-------test for overlapping DE genes for individual clusters----------------------------------------------------------
sham.vs.tac14.de$comp = "sham.vs.tac14.de"
sham.vs.tac28.de$comp = "sham.vs.tac28.de"
tac14.vs.tac28.de$comp = "tac14.vs.tac28.de"
#joined list of DEG
joined.DEG = rbind(sham.vs.tac14.de,sham.vs.tac28.de,tac14.vs.tac28.de)
#filter for only logFC up
joined.DEG.up <- joined.DEG[joined.DEG$avg_logFC>0,]
joined.DEG.up[["order"]]<-ifelse(is.na(c(parse_number(joined.DEG.up$unique))),0,c(parse_number(joined.DEG.up$unique)))
joined.DEG.up <- arrange(joined.DEG.up,comp)
joined.DEG.up <- arrange(joined.DEG.up,cluster)
joined.DEG.up <- arrange(joined.DEG.up,-order)


#--------test genes clusters of DEG from the hclust function used in pheatmap for GOs----------------
#see the heatmap of All DEG form Sham vs TAC and Sham vs Tac28 (All DEG from both TAC conditions)
#cluster genes with hclust
Idents(Samples.combined)="celltypes.stim"
avgexp = AverageExpression(Samples.combined, return.seurat = T)
joined.DEG.up.values <-FetchData(avgexp,vars = unique(joined.DEG.up$gene_name),slot = "scale.data")
joined.DEG.up.values <- cluster.pseudotime.genes(joined.DEG.up.values,keep.hclust=F,k.treecut = 7)
joined.DEG.up.heatmap <- pheatmap(joined.DEG.up.values,main = "All up DEG form all comparisons (k=7 clusters)",cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_heatmap of all up DEG form all comparisons (k=7 clusters).jpeg"), width=1400 , height =1400,quality = 100,res = 200)
print(joined.DEG.up.heatmap)
dev.off()
#save gene as list
write_xlsx(data.frame("gene_name"=rownames(joined.DEG.up.values)),path = paste0(sample,"_list of all DEG corresponding to the k=7 cluster heatmap.xlsx"))

#GO terms of gene clusters
joined.DEG.up.values <-FetchData(avgexp,vars = unique(joined.DEG.up$gene_name),slot = "scale.data")
joined.DEG.up.values <- cluster.pseudotime.genes(joined.DEG.up.values,keep.hclust=T,k.treecut =7)
GO.gene.clusters(joined.DEG.up.values,5)
ggsave(filename = paste0(sample,"GO-terms associated with all up DEG form all comparisons (k=7 clusters).jpeg"), width=4.8, height = 6)
ggsave(filename = paste0(sample,"GO-terms associated with all up DEG form all comparisons (k=7 clusters).svg"), width=4.8, height = 6)

#-------plots after ECM scoring----------------
Samples.combined <- scoreECM(Samples.combined)

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "Core_matrisome1",pt.size = 0,cols = stim.col,split.by = "stim")+
  geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.shape = NA,coef=0,lwd=0.3) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust =1,size = 15),axis.text.y = element_text(size = 10))
ggsave(filename = paste0(sample,"core_matrisome_split_vln.jpeg"), width=5, height = 5)
ggsave(filename = paste0(sample,"core_matrisome_split_vln.svg"), width=5, height = 5)

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "Collagens1",pt.size = 0,cols = stim.col,split.by = "stim")+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust = 0.5))
ggsave(filename = paste0(sample,"collagens_split_vln.jpeg"), width=10, height = 10)

FeaturePlot(Samples.combined, features = "Core_matrisome1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_ECMscore_f.jpeg"), width=10 , height = 10) 
FeaturePlot(Samples.combined, features = "Collagens1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_Collagens1score_f.jpeg"), width=10 , height = 10)

#statistic testing for the ECM scorings
ECM.score.data <- FetchData(Samples.combined,vars = c("Core_matrisome1","celltypes.stim"))
#test overall
k.result <- kruskal.test(ECM.score.data$Core_matrisome1~ECM.score.data$celltypes.stim)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(ECM.score.data$Core_matrisome1,ECM.score.data$celltypes.stim,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result Core_matrisome1score.xlsx"))

col.score.data <- FetchData(Samples.combined,vars = c("Collagens1","celltypes.stim"))
#test overall
k.result <- kruskal.test(col.score.data$Collagens1~col.score.data$celltypes.stim)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(col.score.data$Collagens1,col.score.data$celltypes.stim,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result Collagens1score.xlsx"))

#-------plots after Collagen subgroup scoring----------------------------------
Samples.combined <- scoreCollagens(Samples.combined)
Idents(Samples.combined)=factor(Samples.combined$celltypes.stim,levels = rev(levels(Samples.combined$celltypes.stim)))
DotPlot(Samples.combined,dot.scale = 10,features =  c("fibrillar.col1","network.forming.col1","FACIT.col1","transmembrane.col1","multiplexin.col1"))+theme(axis.text.x = element_text(angle = 45,hjust = 1))+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste0(sample,"Collagen-subgroup_scoring_dp.jpeg"), width=6.3, height = 7.5)
ggsave(filename = paste0(sample,"Collagen-subgroup_scoring_dp.svg"), width=6.3, height = 7.5)

#list of collagens from the scoreCollagens function
col.genes = unlist(col.functional.groups)
collagen.label=col.genes
for (element in collagen.label) {
  if (element %in%fibrillar.col) {collagen.label[collagen.label==element] <- paste0("fibrillar_",element)}
  if (element %in%network.forming.col) {collagen.label[collagen.label==element] <- paste0("network_",element)} 
  if (element %in%FACIT.col) {collagen.label[collagen.label==element] <- paste0("FACIT_",element)}
  if (element %in%transmembrane.col) {collagen.label[collagen.label==element] <- paste0("transmem_",element)}
  if (element %in%multiplexin.col) {collagen.label[collagen.label==element] <- paste0("multiplexin_",element)}
}
names(collagen.label)<-col.genes

#plot expression of collagen subtypes by group or clustered 
Idents(Samples.combined)="celltypes.stim"
avgexp = AverageExpression(Samples.combined, return.seurat = T)

#collagen list reduced to the one that were expressed in the dataset at a relevant amount
col.genes.rel = c("Col12a1","Col14a1","Col15a1","Col16a1","Col1a1","Col1a2","Col27a1","Col3a1","Col4a1","Col4a2","Col4a3","Col4a4","Col4a5","Col5a1","Col5a2","Col5a3","Col6a1","Col6a2","Col6a3","Col6a6","Col8a1","Col8a2")
col.genes = col.genes[col.genes %in% col.genes.rel]
unique.combined.values <-FetchData(avgexp,vars = col.genes,slot = "scale.data")
renamed.df <- as.data.frame(t(unique.combined.values))
renamed.df$names = rownames(renamed.df)
collagen.label = collagen.label[renamed.df$names]
renamed.df[["names"]]=ifelse(renamed.df$names==names(collagen.label),collagen.label,"")
rownames(renamed.df)=renamed.df$names
renamed.df$names=NULL
collagen.subtypes.bytype <- pheatmap(renamed.df,cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_collagen-subgroups_heat_bytype_reduced.jpeg"), width=1200 , height =1200,quality = 100,res = 200)
print(collagen.subtypes.bytype)
dev.off()

#-------plots for MMPs and TIMPs-----------------
#source: https://fibrogenesis.biomedcentral.com/articles/10.1186/1755-1536-5-15
Idents(Samples.combined)="celltypes.stim"
avgexp = AverageExpression(Samples.combined, return.seurat = T)
## TIMPs
TIMPs = c("Timp1","Timp2","Timp3","Timp4")
TIMPs.values <-FetchData(avgexp,vars = TIMPs,slot = "scale.data")
TIMPs.heatmap <- pheatmap(t(TIMPs.values),cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100))
jpeg(filename = paste0(sample,"_heatmap of TIMPs.jpeg"), width=1000 , height =500,quality = 100,res = 200)
print(TIMPs.heatmap)
dev.off()
svg(filename = paste0(sample,"_heatmap of TIMPs.svg"), width=5 , height =2.5)
print(TIMPs.heatmap)
dev.off()

VlnPlot(Samples.combined,features = TIMPs,stack = T,group.by = "celltypes.short",split.by = "stim",flip = T,cols = stim.col)
ggsave(filename = paste0(sample,"timp_stack_vln.jpeg"), width=7, height = 5)

## MMP2
MMPs = c("Mmp1","Mmp3","Mmp8","Mmp13","Mmp2","Mmp9","Mmp12","Mmp28","Mmp14")#"Mmp1" not found
MMPs.values <-FetchData(avgexp,vars = MMPs,slot = "scale.data")
MMPs.heatmap <- pheatmap(t(MMPs.values),cluster_cols = F,cluster_rows = T,angle_col = 45,border_color = 0,color = col.ramp(100))
jpeg(filename = paste0(sample,"_heatmap of MMPs.jpeg"), width=1000 , height =700,quality = 100,res = 200)
print(MMPs.heatmap)
dev.off()
svg(filename = paste0(sample,"_heatmap of MMPs.svg"), width=5 , height =3.5)
print(MMPs.heatmap)
dev.off()

VlnPlot(Samples.combined,features = MMPs,stack = T,group.by = "celltypes.short",split.by = "stim",flip = T,cols = stim.col)
ggsave(filename = paste0(sample,"mmp_stack_vln.jpeg"), width=7, height = 8)

## GO_metalloendopeptidase_activity
#source: http://www.informatics.jax.org/go/term/GO:0004222#myDataTable=results%3D100%26startIndex%3D0%26sort%3Dterm%26dir%3Dasc
GO_metalloendopeptidase_activity <- read_delim("~/Documents/FP_scRNA/R stuff/references/GO_term_metalloendopeptidase activity_GO_0004222", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
met.endo.genes=GO_metalloendopeptidase_activity$Symbol
Samples.combined = AddModuleScore(object = Samples.combined, features = list(met.endo.genes), name = "met.endo.genes", ctrl = 35)
VlnPlot(Samples.combined,group.by = "celltypes.short",features = "met.endo.genes1",pt.size = 0,cols = stim.col,split.by = "stim")+
  geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.shape = NA,coef=0,lwd=0.3) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust =1,size = 15),axis.text.y = element_text(size = 10))+
  ggtitle("GO_metalloendopeptidase_activity")
ggsave(filename = paste0(sample,"met.endo.genes_split_vln.jpeg"), width=5, height = 5)

#-------plots after Stress scoring-----------------------------------
Samples.combined<-scoreStress(Samples.combined)

FeaturePlot(Samples.combined, features = "STRESSscore1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_StressScore_f.jpeg"), width=10 , height = 10) 

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "STRESSscore1",pt.size = 0,cols = cluster.col)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust = 1)) + NoLegend()
ggsave(filename = paste0(sample,"StressScore_vln.jpeg"), width=5, height = 5)


#-------correlation of ECM score to gene expression and dorothea TF activity scores------
#calculated by the script: ECMscore correlation
#selected Tead1 for futher analysis since it had a high score + correlation
#run the ECMscore correlation script first
#get tead1 regulon from dorothea
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))
regulon <- dorothea_regulon_mouse %>%  dplyr::filter(confidence %in% c("A","B","C"))
regulon <- as.data.frame(regulon)
regulon.tead1 <- regulon[regulon$tf=="Tead1",]

#get gene expression from the tead1 target genes that correlate significant with with ECMscore
cors <- cors[-1,]
expression.values <-FetchData(avgexp,vars = rownames(cors)[rownames(cors) %in%regulon.tead1$target & cors$pvalue<0.01],slot = "scale.data")
expression.values=cbind(MeanExprECM, expression.values)
expression.values=arrange(expression.values,Core_matrisome_score)
expression.values.heatmap <- pheatmap(t(expression.values)[1:30,],main = "top30 Tead1 target genes sig. corr. with Core matrisome score",cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20,labels_row = italic_heatmap_labels(t(expression.values)[1:30,]))
jpeg(filename = paste0(sample,"_top30 Tead1 target genes sig. corr. with Core matrisome score.jpeg"), width=1200 , height =1200,quality = 100,res = 200)
print(expression.values.heatmap)
dev.off()

#####################--------seperate integration/clustering of conditions------sham------##############################
# seperate the 3 sham samples from the integration back to individual samples and then reintegrate only the 3 shams
sample = "fib_int_only_sham_"
setwd("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_integrated_all/only_Sham_integrated/")

Idents(Samples.combined)="orig.ident"
col1a1.sham = subset(Samples.combined,idents = "FP19_Col1a1_Sham")
gli1.sham = subset(Samples.combined,idents = "FP13_Gli1_Sham") 
pdgfrb.sham = subset(Samples.combined,idents = "FP3_PDGFRb_Sham") 

#integration
Samples.anchors <- FindIntegrationAnchors(object.list = list(col1a1.sham,gli1.sham,pdgfrb.sham), dims = 1:20)

#exclude tdTom from the anchorlist
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list

#only keep variable features after integration
Samples.combined.sham <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
#normalize and scale RNA data
DefaultAssay(Samples.combined.sham) <- "RNA"
Samples.combined.sham <- NormalizeData(Samples.combined.sham, verbose = FALSE)
Samples.combined.sham <- ScaleData(Samples.combined.sham, verbose = FALSE)
#save
saveRDS(Samples.combined.sham,file = paste0(sample,"only_sham_integration.rds"))

##first clustering after integration
Samples.combined.sham = recluster(Samples.combined.sham,0.3)

#plot for quality check
DimPlot(Samples.combined.sham, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12,cols = cluster.col) + NoLegend()
ggsave(filename = paste0(sample,"quality_check_clusters_res0.5.jpeg"), width=10 , height = 10)

#genes
FeaturePlot(Samples.combined.sham,features = "Postn",order=T,pt.size = 1) + scale_colour_gradientn(colours = gradient.col)
ggsave(filename = paste0(sample,"_postn.jpeg"), width=10, height = 10)

#look were the clusters from the integration of tac and sham end up in only sham
meta <-FetchData(Samples.combined,"celltypes.short",cells = names(Samples.combined$stim[Samples.combined$stim=="Sham"]))
colnames(meta) <-"old_clusters"
Samples.combined.sham <- AddMetaData(Samples.combined.sham,metadata = meta)
DimPlot(Samples.combined.sham,group.by = "old_clusters", reduction = "umap", label = TRUE, pt.size = 1, label.size = 12,cols = cluster.col) + NoLegend()+ggtitle("Clusters from full integration")
ggsave(filename = paste0(sample,"_clusters from full integration projected on only sham.jpeg"), width=10 , height = 10)

#####################--------run velo on the seurat integration of tac and sham--------
setwd("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_integrated_all/RNAvelocity_scanpy/")
#http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html
library(SeuratWrappers)

# If you don't have velocyto's example mouse bone marrow dataset, download with the CURL command
# curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mousesc.forVelo/SCG71.loom', destfile
# = '~/Downloads/SCG71.loom')
ldat.pdgfrb.tac <- ReadVelocity(file = "~/Documents/FP_scRNA/velocyto/loom_files_PVM/FP14_PDGFRb_TAC.loom")
ldat.pdgfrb.sham <- ReadVelocity(file = "~/Documents/FP_scRNA/velocyto/loom_files_PVM/FP3_PDGFRb_Sham.loom")
ldat.gli1.tac <- ReadVelocity(file = "~/Documents/FP_scRNA/velocyto/loom_files_PVM/FP12_Gli1_TAC.loom")
ldat.gli1.sham<- ReadVelocity(file = "~/Documents/FP_scRNA/velocyto/loom_files_PVM/FP13_Gli1_Sham.loom")
ldat.gli1.tac28<- ReadVelocity(file = "~/Documents/FP_scRNA/velocyto/loom_files_PVM/FP71_Gli1_TAC_28d.loom")
ldat.col1a1.tac <- ReadVelocity(file = "~/Documents/FP_scRNA/velocyto/loom_files_PVM/FP18_Col1a1_TAC.loom")
ldat.col1a1.sham <- ReadVelocity(file = "~/Documents/FP_scRNA/velocyto/loom_files_PVM/FP19_Col1a1_Sham.loom")

#we excluded all genes that are known to be influenced by the tissue dissociation, to reduce a potential bias from these genes on the RNA velo analysis
#for the gene set source see also the scoreStress function
stress.markers = c("Fosb", "Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1", "Hspa1a", "Hspa1b", "Hsp90ab1", "Hspa8", "Hspb1", "Ier3", "Ier2", "Btg1", "Btg2", "Dusp1","Brd2", "Dnaja1", "Dnajb1", "Egr1", "Hsp90aa1", "Hsp90ab1", "Hspe1",  "Ier3", "Nr4a1","Junb", "Fosb", "Fos", "Jun", "Zfp36", "Egr1", "Hspa1a", "Hspa1b", "Hspa8", "Hspb1", "Cebpd", "Cebpb", "Atf3", "Socs3", "Jund", "Hspe1", "Hsp90ab1")
stress.markers = unique(stress.markers)

ldat.pdgfrb.tac$ambiguous=NULL
ldat.pdgfrb.tac$spliced <- ldat.pdgfrb.tac$spliced[!(rownames(ldat.pdgfrb.tac$spliced) %in% stress.markers),]
ldat.pdgfrb.tac$unspliced <- ldat.pdgfrb.tac$unspliced[!(rownames(ldat.pdgfrb.tac$unspliced) %in% stress.markers),]

ldat.pdgfrb.sham$ambiguous=NULL
ldat.pdgfrb.sham$spliced <- ldat.pdgfrb.sham$spliced[!(rownames(ldat.pdgfrb.sham$spliced) %in% stress.markers),]
ldat.pdgfrb.sham$unspliced <- ldat.pdgfrb.sham$unspliced[!(rownames(ldat.pdgfrb.sham$unspliced) %in% stress.markers),]

ldat.gli1.tac$ambiguous=NULL
ldat.gli1.tac$spliced <- ldat.gli1.tac$spliced[!(rownames(ldat.gli1.tac$spliced) %in% stress.markers),]
ldat.gli1.tac$unspliced <- ldat.gli1.tac$unspliced[!(rownames(ldat.gli1.tac$unspliced) %in% stress.markers),]

ldat.gli1.sham$ambiguous=NULL
ldat.gli1.sham$spliced <- ldat.gli1.sham$spliced[!(rownames(ldat.gli1.sham$spliced) %in% stress.markers),]
ldat.gli1.sham$unspliced <- ldat.gli1.sham$unspliced[!(rownames(ldat.gli1.sham$unspliced) %in% stress.markers),]

ldat.gli1.tac28$ambiguous=NULL
ldat.gli1.tac28$spliced <- ldat.gli1.tac28$spliced[!(rownames(ldat.gli1.tac28$spliced) %in% stress.markers),]
ldat.gli1.tac28$unspliced <- ldat.gli1.tac28$unspliced[!(rownames(ldat.gli1.tac28$unspliced) %in% stress.markers),]

ldat.col1a1.tac$ambiguous=NULL
ldat.col1a1.tac$spliced <- ldat.col1a1.tac$spliced[!(rownames(ldat.col1a1.tac$spliced) %in% stress.markers),]
ldat.col1a1.tac$unspliced <- ldat.col1a1.tac$unspliced[!(rownames(ldat.col1a1.tac$unspliced) %in% stress.markers),]

ldat.col1a1.sham$ambiguous=NULL
ldat.col1a1.sham$spliced <- ldat.col1a1.sham$spliced[!(rownames(ldat.col1a1.sham$spliced) %in% stress.markers),]
ldat.col1a1.sham$unspliced <- ldat.col1a1.sham$unspliced[!(rownames(ldat.col1a1.sham$unspliced) %in% stress.markers),]

#create seurat and subest the loom file to only the cells that where kept in the privious analysis
Idents(Samples.combined)="orig.ident"

ldat.pdgfrb.tac<- as.Seurat(x = ldat.pdgfrb.tac)
cells_umap = Embeddings(subset(Samples.combined,idents = "FP14_PDGFRb_TAC"),reduction = "umap")
rownames(cells_umap) = gsub("-1_2","",rownames(cells_umap))
rownames(cells_umap) = paste0("FP14_PDGFRb_TAC:",rownames(cells_umap),"x")
ldat.pdgfrb.tac = subset(ldat.pdgfrb.tac,cells = rownames(cells_umap))

ldat.pdgfrb.sham<- as.Seurat(x = ldat.pdgfrb.sham)
cells_umap = Embeddings(subset(Samples.combined,idents = "FP3_PDGFRb_Sham"),reduction = "umap")
rownames(cells_umap) = gsub("-1_1","",rownames(cells_umap))
rownames(cells_umap) = paste0("FP3_PDGFRb_Sham:",rownames(cells_umap),"x")
ldat.pdgfrb.sham = subset(ldat.pdgfrb.sham,cells = rownames(cells_umap))

ldat.gli1.tac<- as.Seurat(x = ldat.gli1.tac)
cells_umap = Embeddings(subset(Samples.combined,idents = "FP12_Gli1_TAC"),reduction = "umap")
rownames(cells_umap) = gsub("-1_6","",rownames(cells_umap))
rownames(cells_umap) = paste0("FP12_Gli1_TAC:",rownames(cells_umap),"x")
ldat.gli1.tac = subset(ldat.gli1.tac,cells = rownames(cells_umap))

ldat.gli1.sham<- as.Seurat(x = ldat.gli1.sham)
cells_umap = Embeddings(subset(Samples.combined,idents = "FP13_Gli1_Sham"),reduction = "umap")
rownames(cells_umap) = gsub("-1_5","",rownames(cells_umap))
rownames(cells_umap) = paste0("FP13_Gli1_Sham:",rownames(cells_umap),"x")
ldat.gli1.sham = subset(ldat.gli1.sham,cells = rownames(cells_umap))

ldat.gli1.tac28<- as.Seurat(x = ldat.gli1.tac28)
cells_umap = Embeddings(subset(Samples.combined,idents = "FP71_Gli1_TAC_28d"),reduction = "umap")
rownames(cells_umap) = gsub("-1_7","",rownames(cells_umap))
rownames(cells_umap) = paste0("21Oct37-DL012-FP71:",rownames(cells_umap),"x")
ldat.gli1.tac28 = subset(ldat.gli1.tac28,cells = rownames(cells_umap))

ldat.col1a1.tac<- as.Seurat(x = ldat.col1a1.tac)
cells_umap = Embeddings(subset(Samples.combined,idents = "FP18_Col1a1_TAC"),reduction = "umap")
rownames(cells_umap) = gsub("-1_4","",rownames(cells_umap))
rownames(cells_umap) = paste0("FP18_Col1a1_TAC:",rownames(cells_umap),"x")
ldat.col1a1.tac = subset(ldat.col1a1.tac,cells = rownames(cells_umap))

ldat.col1a1.sham<- as.Seurat(x = ldat.col1a1.sham)
cells_umap = Embeddings(subset(Samples.combined,idents = "FP19_Col1a1_Sham"),reduction = "umap")
rownames(cells_umap) = gsub("-1_3","",rownames(cells_umap))
rownames(cells_umap) = paste0("FP19_Col1a1_Sham:",rownames(cells_umap),"x")
ldat.col1a1.sham = subset(ldat.col1a1.sham,cells = rownames(cells_umap))


sc.forVelo <-merge(x = ldat.pdgfrb.tac,y = c(ldat.pdgfrb.sham,ldat.gli1.tac,ldat.gli1.sham,ldat.gli1.tac28,ldat.col1a1.tac,ldat.col1a1.sham))

sc.forVelo[["RNA"]] <- sc.forVelo[["spliced"]] #warning is ok

DefaultAssay(sc.forVelo) <- "RNA"

# add the UMAP from the privious analysis
cells_umap = Embeddings(Samples.combined,reduction = "umap")
cell_names = rownames(cells_umap)
cell_name_new = NULL
for (cell_name in cell_names) {
  if (str_detect(cell_name,"-1_2")) {cell_name = gsub("-1_2","",cell_name);cell_name = paste0("FP14_PDGFRb_TAC:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_1")) {cell_name = gsub("-1_1","",cell_name);cell_name = paste0("FP3_PDGFRb_Sham:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_6")) {cell_name = gsub("-1_6","",cell_name);cell_name = paste0("FP12_Gli1_TAC:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_5")) {cell_name = gsub("-1_5","",cell_name);cell_name = paste0("FP13_Gli1_Sham:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_7")) {cell_name = gsub("-1_7","",cell_name);cell_name = paste0("21Oct37-DL012-FP71:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_4")) {cell_name = gsub("-1_4","",cell_name);cell_name = paste0("FP18_Col1a1_TAC:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_3")) {cell_name = gsub("-1_3","",cell_name);cell_name = paste0("FP19_Col1a1_Sham:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
}
rownames(cells_umap) = cell_name_new
sc.forVelo[['umap']] = CreateDimReducObject(embeddings = cells_umap, key = 'UMAP_', assay = DefaultAssay(sc.forVelo))

#add seurat clusters from before
seurat_cluster = data.frame(Samples.combined$seurat_clusters)
cell_names = rownames(seurat_cluster)
cell_name_new = NULL
for (cell_name in cell_names) {
  if (str_detect(cell_name,"-1_2")) {cell_name = gsub("-1_2","",cell_name);cell_name = paste0("FP14_PDGFRb_TAC:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_1")) {cell_name = gsub("-1_1","",cell_name);cell_name = paste0("FP3_PDGFRb_Sham:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_6")) {cell_name = gsub("-1_6","",cell_name);cell_name = paste0("FP12_Gli1_TAC:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_5")) {cell_name = gsub("-1_5","",cell_name);cell_name = paste0("FP13_Gli1_Sham:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_7")) {cell_name = gsub("-1_7","",cell_name);cell_name = paste0("21Oct37-DL012-FP71:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_4")) {cell_name = gsub("-1_4","",cell_name);cell_name = paste0("FP18_Col1a1_TAC:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
  if (str_detect(cell_name,"-1_3")) {cell_name = gsub("-1_3","",cell_name);cell_name = paste0("FP19_Col1a1_Sham:",cell_name,"x");cell_name_new = c(cell_name_new,cell_name)}
}
rownames(seurat_cluster) = cell_name_new
sc.forVelo <- AddMetaData(sc.forVelo,metadata = seurat_cluster,col.name = "seurat_clusters")
Idents(sc.forVelo)="seurat_clusters"
DimPlot(sc.forVelo)


#generat subset for the conditions
table(sc.forVelo$orig.ident)
Idents(sc.forVelo)="orig.ident"
sc.forVelo.sham = subset(sc.forVelo,idents = c("FP3","FP13","FP19"))
sc.forVelo.tac = subset(sc.forVelo,idents = c("FP14","FP12","FP18"))
sc.forVelo.tac28 = subset(sc.forVelo,idents = c("SeuratProject"))

library(SeuratDisk)
#save for further processing in python
#whole integration
SaveH5Seurat(sc.forVelo, filename = "fibroblast_seurat_int_forVelo.h5Seurat")
Convert("fibroblast_seurat_int_forVelo.h5Seurat", dest = "h5ad")
#save tac
SaveH5Seurat(sc.forVelo.tac, filename = "fibroblast_seurat_int_forVelo_TAC.h5Seurat")
Convert("fibroblast_seurat_int_forVelo_TAC.h5Seurat", dest = "h5ad")
#save sham
SaveH5Seurat(sc.forVelo.sham, filename = "fibroblast_seurat_int_forVelo_Sham.h5Seurat")
Convert("fibroblast_seurat_int_forVelo_Sham.h5Seurat", dest = "h5ad")
#save tac28
SaveH5Seurat(sc.forVelo.tac28, filename = "fibroblast_seurat_int_forVelo_TAC28.h5Seurat")
Convert("fibroblast_seurat_int_forVelo_TAC28.h5Seurat", dest = "h5ad")


###########################################################################################################################
#----separating back to individual samples to show batch effect and result of harmony integration----------
#---batch effect
DefaultAssay(Samples.combined)="RNA"
Samples.combined <- FindVariableFeatures(Samples.combined,selection.method = "vst", nfeatures = 2000)
Samples.combined <- recluster_RNA(Samples.combined,res = 0.5)

Samples.combined$orig.geno_stim = paste0(Samples.combined$stim,"_",Samples.combined$orig.geno)
Samples.combined$orig.geno_stim=factor(Samples.combined$orig.geno_stim,levels = c("Sham_Col1a1","Sham_Gli1","Sham_Pdgfrb","TAC_Col1a1","TAC_Gli1","TAC_Pdgfrb","TAC_28_Gli1"))
#cluster plots - genotype
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno_stim",group.by = "orig.geno_stim",pt.size = 1,ncol = 3) + NoLegend()
ggsave(filename = paste0(sample,"_genotype_stim_split.jpeg"), width=20 , height = 20)

#-----run harmony on the data
#start processing
Samples.combined <- Samples.combined %>% RunHarmony("orig.ident", plot_convergence = TRUE,epsilon.cluster = -Inf)
#process based on harmony result
Samples.combined <- Samples.combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"after_harmony_res0.3.jpeg"), width=10 , height = 10)

DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno_stim",group.by = "orig.geno_stim",pt.size = 1,ncol = 3) + NoLegend()
ggsave(filename = paste0(sample,"_genotype_stim_split_harmony.jpeg"), width=20 , height = 20)

BarNormPlot(Samples.combined,"harmony_fib",group2 = "stim")
VlnPlot(Samples.combined,features = "Thbs4",split.by = "stim",pt.size = 0,cols = stim.col)+theme(axis.text.x = element_text(angle = 0))
ggsave(filename = paste0(sample,"after_harmony_Thbs4.jpeg"), width=10 , height = 10)
