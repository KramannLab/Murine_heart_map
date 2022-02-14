#integrate our heart map with MacLellan AngII data with harmony
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
library(reshape2)
library(pheatmap)
library(colorBlindness)
source("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/scripts_for_upload/helper_functions_DEG_GO_ploting.R")

setwd("~/Documents/FP_scRNA/R stuff/MacLellan paper data angII heart scRNA/MacLellan_PVM_Fibro_int")
sample = "MacLellan_PVM_Fibro_int"

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
cluster.col=brewer.pal(name = "Paired",n=12)
cluster.col.Ang = c(brewer.pal(name = "Paired",n=11)[1:6],colorBlindness::PairedColor12Steps)
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)


#----------load fibroblast dataset and merge-----------------
Samples.combined <- readRDS(file ="~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_integrated_all/Fibroblast_integrated__filtered_processed.rds")
fibro.subset <- readRDS("~/Documents/FP_scRNA/R stuff/MacLellan paper data angII heart scRNA/McLellan_fibro_sub_fibroblast_subset_annotated.rds")

#merge Seurat object for harmony
samples.merged <- merge(x = fibro.subset, y = Samples.combined, project = sample)

#first processing
samples.merged <- samples.merged %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = samples.merged@var.genes, npcs = 20, verbose = FALSE)

#run harmony on the data
samples.merged <- samples.merged %>% RunHarmony("orig.ident", plot_convergence = TRUE,epsilon.cluster = -Inf)
#process based on harmony result
samples.merged <- samples.merged %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

#--------save harmony integration----
#factor for annotation
samples.merged$celltypes.short[samples.merged$celltypes.stim %in% c("Fibroblast 1_Sham" , "Fibroblast 1_TAC" , "Fibroblast 1_TAC_28")]="PVM_Fib1"
samples.merged$celltypes.short[samples.merged$celltypes.stim %in% c("Fibroblast 2_Sham" , "Fibroblast 2_TAC" , "Fibroblast 2_TAC_28")]="PVM_Fib2"
samples.merged$celltypes.short[samples.merged$celltypes.stim %in% c("Fibroblast 3_Sham" , "Fibroblast 3_TAC" , "Fibroblast 3_TAC_28")]="PVM_Fib3"
cell.levels = c("PVM_Fib1","StrFib","PVM_Fib2","PVM_Fib3","ECM-Fib","IntFib","Fib1","Fib2","Fib3","Fib4","Fib5","FibCilp","Fib6","FibWif1","FibThbs4","Fib7","Fib8")
samples.merged$celltypes.short<-factor(samples.merged$celltypes.short,levels = cell.levels)
#add stim levels
samples.merged$stim=factor(samples.merged$stim,levels=c("none","saline","angII","Sham","TAC","TAC_28"))
#add celltype + stim
lvl=NULL
for (levels in levels(samples.merged$celltypes.short)) {lvl=c(lvl,paste0(levels,"_",levels(samples.merged$stim)))}
samples.merged$celltypes.stim=factor(paste0(samples.merged$celltypes.short,"_",samples.merged$stim),levels = lvl)
#save
saveRDS(samples.merged,file = paste0(sample,"_after_harmony.rds"))
#load
samples.merged <- readRDS("~/Documents/FP_scRNA/R stuff/MacLellan paper data angII heart scRNA/MacLellan_PVM_Fibro_int/MacLellan_PVM_Fibro_int_after_harmony.rds")

#--------plot umaps----------------
DimPlot(samples.merged,group.by = "seurat_clusters", reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 12,cols = cluster.col) + NoLegend()
ggsave(filename = paste0(sample,"_cluster0.5.jpeg"),width = 10,height = 10)

#-------ploting of multiple featues and save individual------------------------------
#6pack plot
fp.list=NULL
Idents(samples.merged)="stim"
for (stim in levels(samples.merged$stim)) {
  sc.stim <- subset(samples.merged,ident=stim)
  fp.list[[stim]]<-FeaturePlot(sc.stim,features = "Thbs4",order = T,pt.size = 1)+ scale_colour_gradientn(colours =gradient.col)+theme(axis.title = element_blank(),axis.text = element_blank())+ggtitle(stim)
}
Idents(samples.merged)="seurat_clusters"
ggarrange(plotlist = fp.list,common.legend = T,legend = "right")
ggsave(filename = paste("multiple_marker_plots/",sample,"_Thbs4_split.jpeg"), width=15 , height = 10)
#6pack plot
fp.list=NULL
Idents(samples.merged)="stim"
for (stim in levels(samples.merged$stim)) {
  sc.stim <- subset(samples.merged,ident=stim)
  fp.list[[stim]]<-FeaturePlot(sc.stim,features = "Acta2",order = T,pt.size = 1)+ scale_colour_gradientn(colours =gradient.col)+theme(axis.title = element_blank(),axis.text = element_blank())+ggtitle(stim)
}
Idents(samples.merged)="seurat_clusters"
ggarrange(plotlist = fp.list,common.legend = T,legend = "right")
ggsave(filename = paste("multiple_marker_plots/",sample,"_Acta2_split.jpeg"), width=15 , height = 10)
# stack vln of major GOI
VlnPlot(samples.merged,group.by = "celltypes.short",features = c("Thbs4","Postn","Acta2","Wif1","Dkk3","Fgl2","Atf3","Cd248"),split.by = "stim",cols = stim.col.Ang,pt.size = 0,stack = T,flip = T) + theme(axis.title.x = element_blank())
ggsave(filename = paste(sample,"_majorGOI_vln.jpeg"), width=8 , height =6)

FeaturePlot(samples.merged,features = "Wif1",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Wif1.jpeg"), width=10 , height = 10)

FeaturePlot(samples.merged,features = "Dkk3",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Dkk3.jpeg"), width=10 , height = 10)

#-------barplots-----------------------
dittoBarPlot(samples.merged,var = "celltypes.short", group.by = "seurat_clusters",var.labels.reorder = c(14,17,15,16,1,13,2,3,4,5,6,10,7,12,11,8,9),color.panel = cluster.col.Ang,x.labels.rotate = F,main = "Cells of condition per cluster (not normalized)")+theme(axis.text = element_text(size=15))
ggsave(filename = paste0(sample,"barplot_not_norm_anno.jpeg"),width = 5,height = 5)


#-----pearson correlation of subclusters based on the 500 most variable genes used in clustering----------------------
Idents(samples.merged)="celltypes.short"
avgexp = AverageExpression(samples.merged, return.seurat = T)

expression.values <-FetchData(avgexp,vars =samples.merged@assays$RNA@var.features[1:500] ,slot = "scale.data")
expression.values <-t(expression.values)

avgexp.cor = cor(x = expression.values[,1:6],y = expression.values[,7:17])
avgexp.cor.heatmap <- pheatmap(avgexp.cor,main = "pearson_correlation",cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_pearson_correlation_of_subclusters.jpeg"), width=800 , height =500,quality = 100,res = 200)
print(avgexp.cor.heatmap)
dev.off()


