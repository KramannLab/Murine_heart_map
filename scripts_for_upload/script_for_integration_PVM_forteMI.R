#integrate our heart map with forte MI data with harmony
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

setwd("~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/Forte_PVM_Fibro_int/")
sample = "Forte_PVM_Fibro_int_"

#colorblind friendly color panels for all figures
#for conditions/stim #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
sham.col = Okabe_Ito[1]
tac.col = Okabe_Ito[2]
tac28.col = Okabe_Ito[3]
stim.col=c(sham.col,tac.col,tac28.col)
stim.col.light=c("#E69F0080","#56B4E980","#009E7380")
stim.col.MI=c("#F0E442", "#0072B2", "#D55E00",sham.col,tac.col,tac28.col)
#for genotypes #From Paul Tol: https://personal.sron.nl/~pault/ '#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'
col.pdgfrb="#EE6677"
col.col1a1="#228833"
col.gli1="#4477AA"
geno.cl = c(col.col1a1,col.pdgfrb,col.gli1)
geno.cl.light = c("#22883380","#EE667780","#4477AA80")
#for subclustering #RColorBrewer
display.brewer.pal(name = "Paired",n = 12)
cluster.col=brewer.pal(name = "Paired",n=12)
cluster.col.MI = c(brewer.pal(name = "Paired",n=12)[1:6],colorBlindness::PairedColor12Steps)
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)


#----------load fibroblast dataset and merge-----------------
Samples.combined <- readRDS(file ="~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Fibroblast_integrated_all/Fibroblast_integrated__filtered_processed.rds")
fibro.subset <- readRDS("~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/Forte_MI_datafibroblast_subset_annotated.rds")

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
#add general sample groups
# add column for only the genotype
samples.merged$stim=samples.merged$orig.ident
Sham.ls = samples.merged$stim  %in% c("FP56_Myh11_Sham","FP13_Gli1_Sham","FP3_PDGFRb_Sham","FP19_Col1a1_Sham","FP21_Cdh5_Sham","FP23_NG2_Sham")
TAC.ls = samples.merged$stim  %in% c("FP12_Gli1_TAC","FP14_PDGFRb_TAC","FP20_Cdh5_TAC","FP18_Col1a1_TAC","FP22_NG2_TAC","FP57_Myh11_TAC")
TAC.28.ls = samples.merged$stim  %in% c("FP71_Gli1_TAC_28d")
MI.control.ls = samples.merged$stim  %in% c("sc.sham1","sc.sham2","sc.none")
MI.early.ls = samples.merged$stim  %in% c("sc.MIday1","sc.MIday3","sc.MIday5","sc.MIday7")
MI.late.ls = samples.merged$stim  %in% c("sc.MIday14","sc.MIday28")
samples.merged$stim[Sham.ls]="Sham"
samples.merged$stim[TAC.ls]="TAC"
samples.merged$stim[TAC.28.ls]="TAC_28"
samples.merged$stim[MI.control.ls]="MI.control"
samples.merged$stim[MI.early.ls]="MI.early"
samples.merged$stim[MI.late.ls]="MI.late"
samples.merged$stim<-factor(samples.merged$stim,levels = c("MI.control","MI.early","MI.late","Sham","TAC","TAC_28"))
#factor for annotation
cell.levels = c("Fib1","StrFib","Fib2","Fib3","ECM-Fib","IntFib","HEpiD","PLS","Myof","MFC","EndD","IR","Epi","IFNr","DC")
samples.merged$celltypes.short<-factor(samples.merged$celltypes.short,levels = cell.levels)
#add celltype + stim
lvl=NULL
for (levels in levels(samples.merged$celltypes.short)) {lvl=c(lvl,paste0(levels,"_",levels(samples.merged$stim)))}
samples.merged$celltypes.stim=factor(paste0(samples.merged$celltypes.short,"_",samples.merged$stim),levels = lvl)
#save
saveRDS(samples.merged,file = paste0(sample,"_after_harmony.rds"))
#load
samples.merged <- readRDS("~/Documents/FP_scRNA/R stuff/ElviraForte_NRosenthal_MI_scRNA/Forte_PVM_Fibro_int/Forte_PVM_Fibro_int__after_harmony.rds")

#--------plot umap----------------
DimPlot(samples.merged,group.by = "seurat_clusters", reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 12,cols = cluster.col) + NoLegend()
ggsave(filename = paste0(sample,"_cluster0.5.jpeg"),width = 10,height = 10)

#-------ploting of multiple featues and save individual------------------------------
# 6pack plot
fp.list=NULL
Idents(samples.merged)="stim"
for (stim in levels(samples.merged$stim)) {
  sc.stim <- subset(samples.merged,ident=stim)
  fp.list[[stim]]<-FeaturePlot(sc.stim,features = "Thbs4",order = T,pt.size = 1)+ scale_colour_gradientn(colours =gradient.col)+theme(axis.title = element_blank(),axis.text = element_blank())+ggtitle(stim)
}
Idents(samples.merged)="seurat_clusters"
ggarrange(plotlist = fp.list,common.legend = T,legend = "right")
ggsave(filename = paste("multiple_marker_plots/",sample,"_Thbs4_split.jpeg"), width=15 , height = 10)
# 6pack plot
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
VlnPlot(samples.merged,group.by = "celltypes.short",features = c("Thbs4","Postn","Acta2","Wif1","Dkk3","Fgl2","Atf3","Cd248"),split.by = "stim",cols = stim.col.MI,pt.size = 0,stack = T,flip = T) + theme(axis.title.x = element_blank())
ggsave(filename = paste(sample,"_majorGOI_vln.jpeg"), width=8 , height =6)

FeaturePlot(samples.merged,features = "Wif1",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Wif1.jpeg"), width=10 , height = 10)

FeaturePlot(samples.merged,features = "Dkk3",order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Dkk3.jpeg"), width=10 , height = 10)

#-------barplots-----------------------
dittoBarPlot(samples.merged,var = "celltypes.short", group.by = "seurat_clusters",x.reorder =c(1,2,4,5,6,7,8,9,10,11,3),var.labels.reorder = c(5,15,6,7,2,10,8,14,13,12,3,11,4,9,1),color.panel = cluster.col.MI,x.labels.rotate = F,main = "Cells of condition per cluster (not normalized)")+theme(axis.text = element_text(size=15))
ggsave(filename = paste0(sample,"barplot_not_norm_anno.jpeg"),width = 5,height = 5)

#-----pearson correlation of subclusters based on the 500 most variable genes used in clustering----------------------
Idents(samples.merged)="celltypes.short"

avgexp = AverageExpression(samples.merged, return.seurat = T)
expression.values <-FetchData(avgexp,vars =samples.merged@assays$RNA@var.features[1:500] ,slot = "scale.data")
expression.values <-t(expression.values)

avgexp.cor = cor(x = expression.values[,1:6],y = expression.values[,7:15])
avgexp.cor.heatmap <- pheatmap(avgexp.cor,main = "pearson_correlation",cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_pearson_correlation_of_subclusters.jpeg"), width=800 , height =500,quality = 100,res = 200)
print(avgexp.cor.heatmap)
dev.off()


