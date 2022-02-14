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
library(pheatmap)
library(scProportionTest)
source("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/scripts_for_upload/helper_functions_DEG_GO_ploting.R")

setwd("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Cdh5_integrated_all/")
sample="Cdh5_integrated_"

#colorblind friendly color panels for all figures
#for conditions/stim #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
sham.col = Okabe_Ito[1]
tac.col = Okabe_Ito[2]
tac28.col = Okabe_Ito[3]
stim.col=c(sham.col,tac.col,tac28.col)
stim.col.light=c("#E69F0080","#56B4E980","#009E7380")
#for genotypes #From Paul Tol: https://personal.sron.nl/~pault/ '#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'
#for subclustering #RColorBrewer
display.brewer.pal(name = "Paired",n = 12)
cluster.col=brewer.pal(name = "Paired",n=11)
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)
  
################################################----start of sample processing before analysis----##################################################
Samples.combined = readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering_Sham&TAC_integrated/Endothelial_datasets_integrated_after_bacis_filters+TAC28d.rds")
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP20_Cdh5_TAC','TAC',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP21_Cdh5_Sham','Sham',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP72_Cdh5_TAC_28d','TAC_28',Samples.combined@meta.data$stim)
Samples.combined$stim <- factor(Samples.combined$stim,levels=c("Sham","TAC","TAC_28"))
Samples.combined = recluster(Samples.combined,0.5)

#plot before tdTom filtering steps
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = stim.col.light) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"before_processing_clusters_res0.5.jpeg"), width=20 , height = 10)

#ploting only tdtom
#quality check on tdTom expression
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(sample,"tdTom.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(sample,"tdTom_Vln.jpeg"), width=10 , height =10)
RidgePlot(Samples.combined,features = "tdTomato-Green",group.by = "orig.ident",slot = "counts",assay = "RNA")
ggsave(filename = paste0(sample,"tdTom_reads.jpeg"), width=10 , height = 10)

#run Genesorter to get GeneProb of tdtom per cluster
sg = sortGenes(Samples.combined@assays$RNA@data, Idents(Samples.combined))
colnames(sg$condGeneProb) = paste0(levels(Idents(Samples.combined)))
Samples.combined_geneProb = as.data.frame(sg$condGeneProb)["tdTomato-Green",]
Samples.combined_geneProb
remove(sg)

#subsetting to filter out tdtom low and nfeature low clusters
Samples.combined=subset(Samples.combined,idents = c("7","12","13"),invert=T)
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined=subset(Samples.combined,cells = WhichCells(Samples.combined,slot = "counts",expression = `tdTomato-Green` == 0),invert=T)
Samples.combined = recluster(Samples.combined,0.5)

#cluster8 was highly similar to cluster2 in further analysis and therefore merged manualy
types <- list("8"="2",
              "9"="8",
              "10"="9",
              "11"="10")
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$seurat_clusters <- Idents(Samples.combined)
Samples.combined$seurat_clusters = factor(Samples.combined$seurat_clusters,levels = c(0:10))
Idents(Samples.combined)="seurat_clusters"

#plot after tdTom filtering step and before nfeature filter
p1 <- DimPlot(Samples.combined, reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = stim.col.light) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"after_tdTom_low_remove_clusters_res0.5.jpeg"), width=20 , height = 10)

#basic nfeature and ncount plots
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount_vln.jpeg"), width=10 , height =5)

#plot after processing  
DimPlot(Samples.combined, reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"after_processing_clusters_res0.5.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.5.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = stim.col.light) +theme(legend.position = "bottom")
ggsave(filename = paste0(sample,"after_processing_clusters_res0.5_tac_vs_sham.svg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.5_tac_vs_sham.jpeg"), width=10 , height = 10)
#no label version
DimPlot(Samples.combined, reduction = "umap",cols = cluster.col,  pt.size = 1) + NoLegend()
ggsave(filename = paste0(sample,"after_processing_clusters_res0.5_noLabel.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.5_noLabel.svg"), width=10 , height = 10)

#normalize and scale RNA data
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

#save filtered object
saveRDS(Samples.combined,file = paste0(sample,"_filtered_processed.rds"))
################################################----end of sample processing before analysis----##################################################

########################################################----start of sample analysis----###########################################################
#----------load filtered object---------
Samples.combined <- readRDS(file = "~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Cdh5_integrated_all/Cdh5_integrated__filtered_processed.rds")

#--------clustertree for optimazation of cluster resolution------
sc.int=Samples.combined
DefaultAssay(sc.int)<-"integrated"
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)

#--------top10 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
top5 <- setdiff(top5,c("Junb","Ier3","Atf3"))
top5 <- c(top5[1:48],"Hells","Lig1",top5[49:52],"Npr3")
#horizontal top5 after annotation below was added
Idents(Samples.combined)=factor(Samples.combined$celltypes.short,levels = rev(levels(Samples.combined$celltypes.short)))
DotPlot(Samples.combined, features =top5,dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1))+ scale_colour_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(face = "italic"))
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.jpeg"), width=16 , height = 4.3)
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.svg"), width=16 , height = 4.3)
Idents(Samples.combined)=Samples.combined$seurat_clusters


#--------heatmap of EC-atlas paper-------------
#https://www.sciencedirect.com/science/article/pii/S0092867420300623?via%3Dihub
EC_atlas_markers <- read_csv("~/Documents/FP_scRNA/R stuff/references/EC-atlas Heart heatmap table.csv")
markers=as.data.frame(EC_atlas_markers)

avgexp = AverageExpression(Samples.combined, return.seurat = T)
DoHeatmap(avgexp,features = markers$gene_name,draw.lines = F,assay = "RNA",angle = 0) +
  scale_y_discrete(labels=rev(markers$cell_type))  + scale_fill_gradientn(colors = gradient.col)
ggsave(filename = paste0(sample,"_EC-atlas_markers_heatmap_RNA_averages.jpeg"), width=10, height =10)
ggsave(filename = paste0(sample,"_EC-atlas_markers_heatmap_RNA_averages.svg"), width=10, height =10)

#--------annotation for EC clusters according to EC atlas, GO terms and specific markergene expression----------------
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Capillary EC","1"="Capillary artery EC","2"="Capillary vein EC","3"="Stressed EC","4"="Interferon EC",
              "5"="Angiogenic EC","6"="Artery EC","7"="Lymphatic EC","8"="Cycling EC","9"="DNA replicating EC","10"="Endocardial EC")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes <- Idents(Samples.combined)
cell.levels = c("Capillary EC","Capillary artery EC","Capillary vein EC","Stressed EC","Interferon EC","Angiogenic EC","Artery EC","Lymphatic EC","Cycling EC","DNA replicating EC","Endocardial EC")
Samples.combined$celltypes = factor(Samples.combined$celltypes,levels = cell.levels)
#short abroviation of annotation
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="CapEC","1"="CapA-EC","2"="CapV-EC","3"="StrEC","4"="IntEC","5"="AngEC","6"="ArtEC","7"="LymEC","8"="CyclEC","9"="RepEC","10"="EndoEC")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes.short <- Idents(Samples.combined)
cell.levels = c("CapEC","CapA-EC","CapV-EC","StrEC","IntEC","AngEC","ArtEC","LymEC","CyclEC","RepEC","EndoEC")
Samples.combined$celltypes.short = factor(Samples.combined$celltypes.short,levels = cell.levels)
Idents(Samples.combined)="seurat_clusters"
#add celltypes + stim with levels
lvl=NULL
for (levels in levels(Samples.combined$celltypes)) {lvl=c(lvl,paste0(levels,"_",levels(Samples.combined$stim)))}
Samples.combined$celltypes.stim=factor(paste0(Samples.combined$celltypes,"_",Samples.combined$stim),levels = lvl)
#save after renaming
saveRDS(Samples.combined,file = paste0(sample,"_filtered_processed.rds"))

#--------feature plots of interessting genes----------------------------
#vln plot for the endocard markers
v1 <- VlnPlot(Samples.combined,group.by = "celltypes.short",cols = cluster.col,features = c("Npr3"),pt.size = 0)+NoLegend() + theme(axis.text.x = element_text(angle = 45,hjust =1),axis.title.x = element_blank())
v2 <- VlnPlot(Samples.combined,group.by = "celltypes.short",cols = cluster.col,features = c("Vwf"),pt.size = 0)+NoLegend() + theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.title.x = element_blank())
ggarrange(v1,v2,ncol = 1)
ggsave(filename = paste0(sample,"_Npr3_and_Vwf_for_endocard_cluster10.jpeg"), width=5 , height = 5)
ggsave(filename = paste0(sample,"_Npr3_and_Vwf_for_endocard_cluster10.svg"), width=5 , height = 5)
#feature plot for the endocard markers
FeaturePlot(Samples.combined, features = "Npr3", min.cutoff = "q9",pt.size = 1,order = T)  + scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste0(sample,"Npr3.jpeg"), width=10 , height = 10)


VlnPlot(Samples.combined,features = c("Hif1a","Eno1","Epas1"),split.by = "stim",cols = stim.col,pt.size = 0,stack = T,flip = T)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Eno1_Hif1a_Hi2a_Vln_split.jpeg"), width=10 , height = 10)

VlnPlot(Samples.combined,features = "Tinagl1",split.by = "stim",cols = stim.col,pt.size = 0)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Tinagl1_Vln_split.jpeg"), width=10 , height = 10)

VlnPlot(Samples.combined,features = "Trpv4",split.by = "stim",cols = stim.col,pt.size = 0)
ggsave(filename = paste("multiple_marker_plots/",sample,"_Trpv4_Vln_split.jpeg"), width=10 , height = 10)

VlnPlot(Samples.combined,features = "tdTomato-Green",cols = cluster.col,pt.size = 0)+NoLegend()
ggsave(filename = paste("multiple_marker_plots/",sample,"_tdTomato-Green_Vln.jpeg"), width=10 , height = 10)


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
cell.levels =  c("CapEC","CapA-EC","CapV-EC","StrEC","IntEC","AngEC","ArtEC","LymEC","CyclEC","RepEC","EndoEC")

ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=c(cell.levels,cell.levels,cell.levels),y=1),hjust = 0, color="black", size=10)+
  scale_fill_manual(values = rev(stim.col))+
  coord_flip() + 
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = rev(levels(factor(breakdown.df$Var1))),name="Cluster")+
  scale_y_continuous(name="Percentage of total cells [%]") +
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

#--------bar graph for percentage of cells per cluster (of total cells)------------------
breakdown<-round(prop.table(table(Samples.combined@meta.data$seurat_clusters))*100,digits = 1)
breakdown<-(as.data.frame(breakdown))
cell.levels = c("Capillary EC","Capillary artery EC","Capillary vein EC","Stressed EC","Interferon EC","Angiogenic EC","Artery EC","Lymphatic EC","Cycling EC","DNA replicating EC","Endocardial EC")
ggplot(data=breakdown, aes(x=Var1, y=Freq,fill=Var1))+
  geom_bar(stat="identity")+
  geom_text(aes(label=cell.levels,y=1),hjust = 0, color="black", size=10)+
  scale_fill_manual(values = cluster.col)+
  coord_flip() + 
  scale_x_discrete(limits=rev(breakdown[,1])) + 
  scale_y_continuous(name="Percentage of total cells [%]") +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  NoLegend()
ggsave(filename = paste0(sample,"_barplot_percentage_per_cluster.jpeg"), width=5 , height = 10)
ggsave(filename = paste0(sample,"_barplot_percentage_per_cluster.svg"), width=5 , height = 10)



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

DimPlot(Samples.combined,group.by = "Phase",pt.size = 1)
ggsave(filename = paste0(sample,"_cc_Phase.jpeg"), width=10 , height = 10)

#--------get GO-terms of the top100 markergenes per cluster-------------
go.of.marker.genes.100 <- get.markergenes.GOs(Samples.combined,sample,create.folder = F,marker.genes.nr = 100,all.markers = all.markers)
plot.markergenes.GOs(sample,go.of.marker.genes.100,marker.genes.nr = 100,provide.label = types)
plot.markergenes.GOs(sample,go.of.marker.genes.100,marker.genes.nr = 100,topGOs.number = 2,custom.h = 4,custom.w = 5.6,provide.label = types)

#--------DEG testing-------------------------------------------------
Idents(Samples.combined)="stim"
#Sham vs TAC_14
sample="Cdh5_integrated_shamvstac14_"
sham.vs.tac14 <-subset(Samples.combined,idents=c("Sham","TAC"))
sham.vs.tac14.de <- get.DEG.between.conditions.per.cluster(sham.vs.tac14,sample,control_sample = "Sham",treatment_sample = "TAC")
sham.vs.tac14.de <- check.for.unique.DE(sham.vs.tac14.de)

#Sham vs TAC_28
sample="Cdh5_integrated_shamvstac28_"
sham.vs.tac28 <-subset(Samples.combined,idents=c("Sham","TAC_28"))
sham.vs.tac28.de <- get.DEG.between.conditions.per.cluster(sham.vs.tac28,sample,control_sample = "Sham",treatment_sample = "TAC_28")
sham.vs.tac28.de <- check.for.unique.DE(sham.vs.tac28.de)

#TAC_14 vs TAC_28
sample="Cdh5_integrated_tac14vstac28_"
tac14.vs.tac28 <-subset(Samples.combined,idents=c("TAC","TAC_28"))
tac14.vs.tac28.de <- get.DEG.between.conditions.per.cluster(tac14.vs.tac28,sample,control_sample = "TAC",treatment_sample = "TAC_28")
tac14.vs.tac28.de <- check.for.unique.DE(tac14.vs.tac28.de)

#--------barplot number of DEG-------
number.of.DE.up = NULL
number.of.DE.down = NULL
number.of.DE.df=NULL
for (cluster.nr in levels(factor(sham.vs.tac14.de$cluster))) {
  number.of.DE.up=length(sham.vs.tac14.de$gene_name[sham.vs.tac14.de$cluster==cluster.nr & sham.vs.tac14.de$avg_logFC>0 & sham.vs.tac14.de$p_val_adj<0.01])
  number.of.DE.down=length(sham.vs.tac14.de$gene_name[sham.vs.tac14.de$cluster==cluster.nr & sham.vs.tac14.de$avg_logFC<0 & sham.vs.tac14.de$p_val_adj<0.01])
  number.of.DE.df.new= data.frame("cluster"=cluster.nr, "DE.up"=number.of.DE.up,"DE.down"=number.of.DE.down)
  number.of.DE.df=rbind(number.of.DE.df.new,number.of.DE.df)
}
number.of.DE.df$cond="TAC_14"
df.14 <-number.of.DE.df

number.of.DE.up = NULL
number.of.DE.down = NULL
number.of.DE.df=NULL
for (cluster.nr in levels(factor(sham.vs.tac28.de$cluster))) {
  number.of.DE.up=length(sham.vs.tac28.de$gene_name[sham.vs.tac28.de$cluster==cluster.nr & sham.vs.tac28.de$avg_logFC>0 & sham.vs.tac28.de$p_val_adj<0.01])
  number.of.DE.down=length(sham.vs.tac28.de$gene_name[sham.vs.tac28.de$cluster==cluster.nr & sham.vs.tac28.de$avg_logFC<0 & sham.vs.tac28.de$p_val_adj<0.01])
  number.of.DE.df.new= data.frame("cluster"=cluster.nr, "DE.up"=number.of.DE.up,"DE.down"=number.of.DE.down)
  number.of.DE.df=rbind(number.of.DE.df.new,number.of.DE.df)
}
number.of.DE.df$cond="TAC_28"
df.28 <-number.of.DE.df

number.of.DE.df<-rbind(df.14,df.28)
number.of.DE.df<-arrange(number.of.DE.df,cluster)
number.of.DE.df<-number.of.DE.df[c(1:4,7:22,5,6),]
number.of.DE.df$cluster<-factor(number.of.DE.df$cluster,levels = c(0:10))

ggplot(number.of.DE.df, aes(x = cluster)) + 
  geom_bar(aes(y = DE.up,fill="up"), stat = "identity",position = "dodge2") + 
  geom_bar(aes(y = -DE.down,fill="down"), stat = "identity",position =  "dodge2") +
  scale_fill_manual(values =c("#3086e1","#e14f30")) +
  scale_y_continuous(labels = c(rev(seq(0,max(number.of.DE.df$DE.down)+5,50)),seq(0,max(number.of.DE.df$DE.up)+5,50)),
                     breaks = c(rev(seq(0,max(number.of.DE.df$DE.down)+5,50))*-1,seq(0,max(number.of.DE.df$DE.up)+5,50)))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.title = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("Cluster") +ylab("number of DE genes") + ggtitle("Number of DEG (up/down) per cluster (Sham vs. TAC_14 or TAC_28")+theme(title = element_text(size = 8))

ggsave(filename = paste0(sample,"number of Up&Down DE genes per cluster.jpeg"), width=5, height = 5)
ggsave(filename = paste0(sample,"number of Up&Down DE genes per cluster.svg"), width=5, height = 5)



#--------test for overlapping DE genes----------------------------------------------------------
#All DEG form Sham vs TAC and Sham vs Tac28 (All DEG from both TAC conditions)
joined.DEG.tacs = rbind(sham.vs.tac14.de,sham.vs.tac28.de)
Idents(Samples.combined)="celltypes.stim"
avgexp.joined.DEG.tacs = AverageExpression(Samples.combined, return.seurat = T)
joined.DEG.tacs.values <-FetchData(avgexp.joined.DEG.tacs,vars = unique(joined.DEG.tacs$gene_name),slot = "scale.data")
joined.DEG.tacs.heatmap <- pheatmap(t(joined.DEG.tacs.values),main ="All DEG from both TAC conditions" ,cluster_cols = F,angle_col = 45,border_color = 0,color = col.ramp(100))
jpeg(filename = paste0(sample,"_heatmap of all DEG from both TAC conditions.jpeg"), width=1800 , height =1400,quality = 100,res = 200)
print(joined.DEG.tacs.heatmap)
dev.off()


#--------test genes clusters of DEG from the hclust function used in pheatmap for GOs----------------
#see the heatmap of All DEG form Sham vs TAC and Sham vs Tac28 (All DEG from both TAC conditions)
#cluster genes with hclust
joined.DEG.tacs.values <-FetchData(avgexp.joined.DEG.tacs,vars = unique(joined.DEG.tacs$gene_name),slot = "scale.data")
joined.DEG.tacs.values <- cluster.pseudotime.genes(joined.DEG.tacs.values,keep.hclust=F,k.treecut = 5)
joined.DEG.tacs.heatmap <- pheatmap(joined.DEG.tacs.values,main = "All DEG from both TAC conditions (k=5 clusters)",cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_heatmap of all DEG from both TAC conditions_(k=5 clusters).jpeg"), width=1800 , height =1400,quality = 100,res = 200)
print(joined.DEG.tacs.heatmap)
dev.off()
#GO terms of gene clusters
joined.DEG.tacs.values <-FetchData(avgexp.joined.DEG.tacs,vars = unique(joined.DEG.tacs$gene_name),slot = "scale.data")
joined.DEG.tacs.values <- cluster.pseudotime.genes(joined.DEG.tacs.values,keep.hclust=T,k.treecut = 5)
GO.gene.clusters(joined.DEG.tacs.values,5)
ggsave(filename = paste0(sample,"GO-terms associated with all DEG from both TAC conditions_(k=5 clusters).jpeg"), width=5.8, height = 5)
ggsave(filename = paste0(sample,"GO-terms associated with all DEG from both TAC conditions_(k=5 clusters).svg"), width=5.8, height = 5)

#--------plot after migration scoring---------
mig_geneset <- read_delim("~/Documents/FP_scRNA/R stuff/references/positive regulation of cell migration GO_term_summary_20211130_085115.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
mig_geneset <-unlist(mig_geneset$Symbol,use.names = F)
Samples.combined<-AddModuleScore(Samples.combined,features = list(mig_geneset),ctrl = "35",name = "pos_migration")

FeaturePlot(Samples.combined, features = "pos_migration1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col) + ggtitle("positive regulation of cell migration")
ggsave(filename = paste0(sample,"_pos_migration_f.jpeg"), width=10 , height = 10) 

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "pos_migration1",pt.size = 0,cols = stim.col,split.by = "stim")+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust = 1)) + NoLegend()+ ggtitle("positive regulation of cell migration")
ggsave(filename = paste0(sample,"pos_migration_vln.jpeg"), width=8, height = 5)
ggsave(filename = paste0(sample,"pos_migration_vln.svg"), width=8, height = 5)

#statistic test
migration.score.data <- FetchData(Samples.combined,vars = c("pos_migration1","celltypes.stim"))
#test overall
k.result <- kruskal.test(migration.score.data$pos_migration1~migration.score.data$celltypes.stim)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(migration.score.data$pos_migration1,migration.score.data$celltypes.stim,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result pos_migration1score.xlsx"))

#--------plot after EndMA scoring--------------
Samples.combined <- scoreEndMA(Samples.combined)

FeaturePlot(Samples.combined, features = "EndMAscore1",pt.size = 0.5,order = T) + scale_colour_gradientn(colours = gradient.col) + ggtitle("EndMAscore1")
ggsave(filename = paste0(sample,"_EndMAscore_f.jpeg"), width=10 , height = 10) 

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "EndMAscore1",pt.size = 0,cols = stim.col,split.by = "stim")+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust = 1)) + NoLegend()+ ggtitle("EndMAscore1")
ggsave(filename = paste0(sample,"EndMAscore1_vln.jpeg"), width=8, height = 5)


#-------plots for MMPs and TIMPs-----------------
#source: https://fibrogenesis.biomedcentral.com/articles/10.1186/1755-1536-5-15
Idents(Samples.combined)="celltypes.stim"
avgexp = AverageExpression(Samples.combined, return.seurat = T)
## TIMPs
TIMPs = c("Timp1","Timp2","Timp3","Timp4")
TIMPs.values <-FetchData(avgexp,vars = TIMPs,slot = "scale.data")
TIMPs.heatmap <- pheatmap(t(TIMPs.values),cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100))
jpeg(filename = paste0(sample,"_heatmap of TIMPs.jpeg"), width=1800 , height =500,quality = 100,res = 200)
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
jpeg(filename = paste0(sample,"_heatmap of MMPs.jpeg"), width=1800 , height =700,quality = 100,res = 200)
print(MMPs.heatmap)
dev.off()
svg(filename = paste0(sample,"_heatmap of MMPs.svg"), width=5 , height =3.5)
print(MMPs.heatmap)
dev.off()

VlnPlot(Samples.combined,features = c("Mmp13","Mmp2","Mmp9","Mmp28","Mmp14"),stack = T,group.by = "celltypes.short",split.by = "stim",flip = T,cols = stim.col)
ggsave(filename = paste0(sample,"mmp_stack_vln.jpeg"), width=7, height = 6)

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
ggsave(filename = paste0(sample,"met.endo.genes_split_vln.jpeg"), width=10, height = 5)



