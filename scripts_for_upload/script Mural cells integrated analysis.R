#script for analysis of Mural Sham and TAC integrated
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

setwd("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Mural_integrated_all/")
sample="Mural_integrated_"


#colorblind friendly color panels for all figures
#for conditions/stim #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
sham.col = Okabe_Ito[1]
tac.col = Okabe_Ito[2]
tac28.col = Okabe_Ito[3]
stim.col=c(sham.col,tac.col)
stim.col.light=c("#E69F0080","#56B4E980")
#for genotypes #From Paul Tol: https://personal.sron.nl/~pault/ '#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'
col.pdgfrb="#EE6677"
col.myh11="#CCBB44"
col.ng2="#66CCEE"
geno.cl = c(col.pdgfrb,col.myh11,col.ng2)
geno.cl.light = c("#EE667780","#CCBB4480","#66CCEE80")
#for subclustering #RColorBrewer
display.brewer.pal(name = "Paired",n = 12)
cluster.col=brewer.pal(name = "Paired",n=11)
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)

################################################----start of sample processing before analysis----##################################################
Samples.combined = readRDS(file = "~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering_Sham&TAC_integrated/Mural_datasets_integrated_after_bacis_filters.rds")
#add stim based on sha, tac or tac28d
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP22_NG2_TAC','TAC',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP14_PDGFRb_TAC','TAC',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP57_Myh11_TAC','TAC',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP23_NG2_Sham','Sham',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP3_PDGFRb_Sham','Sham',Samples.combined@meta.data$stim)
Samples.combined[['stim']] <- ifelse(Samples.combined@meta.data$stim == 'FP56_Myh11_Sham','Sham',Samples.combined@meta.data$stim)
Samples.combined$stim <- factor(Samples.combined$stim,levels=c("Sham","TAC"))
#add genotype to metadata
# add colum for only the genotype
Samples.combined$orig.geno=Samples.combined$orig.ident
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("FP22_NG2_TAC","FP23_NG2_Sham")]="Ng2"
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("FP14_PDGFRb_TAC","FP3_PDGFRb_Sham")]="Pdgfrb"
Samples.combined$orig.geno[Samples.combined@meta.data$orig.ident %in% c("FP57_Myh11_TAC","FP56_Myh11_Sham")]="Myh11"
Samples.combined$orig.geno=factor(Samples.combined$orig.geno,levels = c("Pdgfrb","Myh11","Ng2"))
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
#add celltypes + stim with levels
lvl=NULL
for (levels in levels(Samples.combined$orig.geno)) {lvl=c(lvl,paste0(levels,"_",levels(Samples.combined$stim)))}
Samples.combined$celltypes.stim=factor(paste0(Samples.combined$orig.geno,"_",Samples.combined$stim),levels = lvl)
DimPlot(Samples.combined,split.by = "celltypes.stim", reduction = "umap", pt.size = 1,ncol = 2) + NoLegend()
ggsave(filename = paste0(sample,"_sample_split.jpeg"), width=10 , height = 15)

#-------------ploting only tdtom
#quality check on tdTom expression
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 0.5)
ggsave(filename = paste0(sample,"tdTom.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0)
ggsave(filename = paste0(sample,"tdTom_Vln_before.jpeg"), width=10 , height =10)
RidgePlot(Samples.combined,features = "tdTomato-Green",group.by = "orig.ident",slot = "counts",assay = "RNA")
ggsave(filename = paste0(sample,"tdTom_reads.jpeg"), width=10 , height = 10)

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
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_before.jpeg"), width=10 , height = 10)

#--------top10 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#get top10 markergenes per cluster
top10 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene))
#vertikal
DotPlot(Samples.combined, features = top10,dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"_top10_markers_unfiltered.jpeg"), width=7 , height = 25)

#----------subsetting to filter out nfeature low cluster and a remaining cluster of contaminating EC-----------
Samples.combined=subset(Samples.combined,idents = c("6","7"),invert=T)
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined=subset(Samples.combined,cells = WhichCells(Samples.combined,slot = "counts",expression = `tdTomato-Green` == 0),invert=T)
Samples.combined = recluster(Samples.combined,0.7)

#------reordering clusternumbers for a more comprehensible
types <- list("0"="P2","1"="P1","2"="sP","3"="V1","4"="V2","5"="SW","6"="sV","7"="iP")
Samples.combined <- RenameIdents(Samples.combined, types)
types <- list("P1"="0","P2"="1","V2"="2","V1"="3","sV"="4","sP"="5","iP"="6","SW"="7")
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$seurat_clusters <- Idents(Samples.combined)
Samples.combined$seurat_clusters = factor(Samples.combined$seurat_clusters,levels = c(0:7))
Idents(Samples.combined)="seurat_clusters"

#plot after tdTom filtering step and before nfeature filter
p1 <- DimPlot(Samples.combined, reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = stim.col.light) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"_filtered_clusters_res0.7.jpeg"), width=20 , height = 10)

#------------basic nfeature ncount plots and tdtom------------------
FeaturePlot(Samples.combined, features = c('nFeature_RNA', 'nCount_RNA'))
ggsave(filename = paste0(sample,"_nfeature_ncount.jpeg"), width=10 , height = 5)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA', 'nCount_RNA'), pt.size = 0.1) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_ncount_vln.jpeg"), width=10 , height =5)
VlnPlot(Samples.combined, features = c("tdTomato-Green"), pt.size = 0,cols = cluster.col)
ggsave(filename = paste0(sample,"tdTom_Vln_after.jpeg"), width=10 , height =10)
VlnPlot(Samples.combined,group.by = "seurat_clusters", features =c('nFeature_RNA'), pt.size = 0,cols = cluster.col) +NoLegend()
ggsave(filename = paste0(sample,"_nfeature_vln_after.jpeg"), width=10 , height =10)

#plot after processing  
DimPlot(Samples.combined, reduction = "umap",cols = cluster.col, label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"after_processing_clusters_res0.7.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.7.svg"), width=10 , height = 10)
#no label version
DimPlot(Samples.combined, reduction = "umap",cols = cluster.col,  pt.size = 1) + NoLegend()
ggsave(filename = paste0(sample,"after_processing_clusters_res0.7_noLabel.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.7_noLabel.svg"), width=10 , height = 10)

#-------normalize and scale RNA data------------------------
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

#--------save filtered object-----------
saveRDS(Samples.combined,file = paste0(sample,"_filtered_processed.rds"))
################################################----end of sample processing before analysis----##################################################
########################################################----start of sample analysis----###########################################################
#----------load filtered object---------
Samples.combined <- readRDS(file ="~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Mural_integrated_all/Mural_integrated__filtered_processed.rds")

#clustertree for optimazation of cluster resolution
sc.int=Samples.combined
DefaultAssay(sc.int)<-"integrated"
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)

#--------annotation for Fib clusters according to marker genes, GO terms and ECM scoring
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Pericytes 1","1"="Pericytes 2","2"="VSMC 2","3"="VSMC 1","4"="Stressed VSMC","5"="Stessed Pericytes","6"="Interferon Pericytes","7"="Schwann cells")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes <- Idents(Samples.combined)
cell.levels = c("Pericytes 1","Pericytes 2","VSMC 2","VSMC 1","Stressed VSMC","Stessed Pericytes","Interferon Pericytes","Schwann cells")
Samples.combined$celltypes = factor(Samples.combined$celltypes,levels = cell.levels)
#short abroviation of annotation
Idents(Samples.combined)="seurat_clusters"
types <- list("0"="Peri1","1"="Peri2","2"="VSMC2","3"="VSMC1","4"="StrVSMC","5"="StrPeri","6"="IntPeri","7"="Sw")
#Rename identities
Samples.combined <- RenameIdents(Samples.combined, types)
Samples.combined$celltypes.short <- Idents(Samples.combined)
cell.levels = c("Peri1","Peri2","VSMC2","VSMC1","StrVSMC","StrPeri","IntPeri","Sw")
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
ggsave(filename = paste0(sample,"after_processing_clusters_res0.3_labeled_short.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"after_processing_clusters_res0.3_labeled_short.svg"), width=10 , height = 10)

#--------top10 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
#horizontal
Idents(Samples.combined) = factor(Samples.combined$celltypes.short,levels = rev(levels(Samples.combined$celltypes.short))) #see lower section for annotation
DotPlot(Samples.combined,features = top5,dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1)) + scale_colour_gradientn(colours = gradient.col) +theme(axis.text.x = element_text(face = "italic"))
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_h.jpeg"), width=15, height = 3.7)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_h.svg"), width=15 , height = 3.7)

#feature plots of pericytes vs VSMC markers for the figure
fp.list=NULL
for (gene in c("Myh11","Acta2","Notch3","Kcnj8","Abcc9","Colec11")) {
  fp.list[[gene]]<-FeaturePlot(Samples.combined,features = gene,order = T,pt.size = 0.5)+ scale_colour_gradientn(colours =gradient.col)+theme(line = element_blank(),axis.title = element_blank(),axis.text = element_blank())+ggtitle(gene)
}
ggarrange(plotlist = fp.list,common.legend = T,legend = "right")
ggsave(filename = paste(sample,"_pericyte_VSMC_marker.jpeg"), width=15 , height = 10)
ggsave(filename = paste(sample,"_pericyte_VSMC_marker.svg"))

#----barplot of genotype contribution per cluster (normalized to even input per genotype)-----
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$orig.geno) 
total.myh11 = sum(breakdown[,2])
total.ng2 = sum(breakdown[,3])
total.pdgfrb = sum(breakdown[,1])
ratio.pdgfrb.ng2 = total.pdgfrb/total.ng2
ratio.pdgfrb.myh11 = total.pdgfrb/total.myh11
breakdown[,2] = round(breakdown[,2]*ratio.pdgfrb.myh11)
breakdown[,3] = round(breakdown[,3]*ratio.pdgfrb.ng2)
breakdown=t(breakdown)
breakdown <- round(apply(breakdown, 2, function(x){x*100/sum(x)}),2)
breakdown.df = as.data.frame(breakdown)
breakdown.df = melt(t(breakdown.df))
breakdown.df$Var2 = factor(breakdown.df$Var2,levels = levels(Samples.combined$orig.geno))
ggplot(data=breakdown.df, aes(x=Var1, y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = levels(factor(breakdown.df$Var1)),name="Cluster")+
  scale_y_continuous(name="Normalizes percentage of contribution [%]") +
  scale_fill_manual(values = c("Myh11" = col.myh11,  "Pdgfrb" = col.pdgfrb,"Ng2" = col.ng2))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.y = element_text(size=20,face = "bold"),axis.title.y = element_text(size = 15),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15),legend.title = element_blank(),legend.position = c(0.92,0.92),legend.key.size = unit(0.5, "cm"),legend.box.background = element_rect(colour = "black"))
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster_normalized.jpeg"), width=10 , height =10)
ggsave(filename = paste0(sample,"_Sample_contribution_per_custom_cluster_normalized.svg"), width=10 , height =10)

#-------barplot of condition contribution per cluster (input normalized)-----
breakdown<-table(Samples.combined@meta.data$seurat_clusters, Samples.combined@meta.data$stim) 
breakdown[,1]=100*(breakdown[,1]/sum(breakdown[,1]))
breakdown[,2]=100*(breakdown[,2]/sum(breakdown[,2]))
breakdown=t(breakdown)
breakdown <- apply(breakdown, 2, function(x){x*100/sum(x)})
breakdown = round(as.data.frame(breakdown),digits = 2)
breakdown.df = as.data.frame(breakdown)
breakdown.df = melt(t(breakdown.df))
breakdown.df$Var2 = factor(breakdown.df$Var2,levels = c("TAC","Sham"))
celltypes.short =   c("Peri1","Peri2","VSMC2","VSMC1","StrVSMC","StrPeri","IntPeri","Sw")

ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=c(celltypes.short,celltypes.short),y=1),hjust = 0, color="black", size=10)+
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
## save summary
write_xlsx(results.sham.vs.tac,path = paste0(sample,"_proportionTest_results.xlsx"))


##------------get GO-terms of markergenes per cluster-------------
types <- list("0"="Peri1","1"="Peri2","2"="VSMC2","3"="VSMC1","4"="StrVSMC","5"="StrPeri","6"="IntPeri","7"="Sw")
go.of.marker.genes <- get.markergenes.GOs(Samples.combined,sample,create.folder = F,all.markers = all.markers,marker.genes.nr = 50)
plot.markergenes.GOs(sample,go.of.marker.genes,provide.label = types)
plot.markergenes.GOs(sample,go.of.marker.genes,topGOs.number = 2,provide.label = types,custom.h = 3.2,custom.w = 7)

#find unique GO-terms of mural subtypes
Idents(Samples.combined)="celltypes.short"
mural.sub = subset(Samples.combined,ident=c("Peri1","Peri2","VSMC1","VSMC2"))
Idents(mural.sub)="seurat_clusters"
all.markers = FindAllMarkers(mural.sub, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#result of findallmarkers
saveRDS(all.markers,file = "mural.sub_all.markers_markergenes_per_cluster_result.RDS")
go.of.marker.genes <- get.markergenes.GOs(sc.object = mural.sub,sample = paste0(sample,"mural_sub_"),create.folder = F,all.markers = all.markers)
plot.markergenes.GOs(sample = paste0(sample,"mural_sub_"),go.of.marker.genes,provide.label = list("0"="Peri1","1"="Peri2","2"="VSMC2","3"="VSMC1"))

plot.markergenes.GOs(paste0(sample,"mural_sub_"),go.of.marker.genes,provide.label = types,uniqueGO.filter=T)
plot.markergenes.GOs(paste0(sample,"mural_sub_"),go.of.marker.genes,topGOs.number = 5,provide.label = types,uniqueGO.filter=T,custom.h = 4,custom.w = 5.4)

#---------------DEG testing-------------------------------------------------
Idents(Samples.combined)="stim"
types <- list("0"="Peri1","1"="Peri2","2"="VSMC2","3"="VSMC1","4"="StrVSMC","5"="StrPeri","6"="IntPeri","7"="Sw")
#Sham vs TAC_14
sample="Mural_integrated_shamvstac14_"
sham.vs.tac14.de <- get.DEG.between.conditions.per.cluster(Samples.combined,sample,control_sample = "Sham",treatment_sample = "TAC")
sham.vs.tac14.de <- check.for.unique.DE(sham.vs.tac14.de)

#test for gene set enrichment in DEG after grouping them as cluster specific or shared
GOs.of.unique_shared.DEG <- get.GOs.based.unique_shared.DEG(sham.vs.tac14.de)
plot.GOs.based.unique_shared.DEG(go.results = GOs.of.unique_shared.DEG,sample,treatment_sample = "TAC",provide.label = types,custom.h = 4,custom.w = 5.7)


#volcano plots for Peri1 & Peri2 of Sham vs. TAC14 as most interesting cluster of mural cells
sham.vs.tac14.de <- readRDS("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/Mural_integrated_all/DE_test_conditions_per_cluster/Mural_integrated_shamvstac14__test:MAST_FindMarkers_Sham_vs_TAC.rds")
sham.vs.tac14.de <- check.for.unique.DE(sham.vs.tac14.de)
#pericyte1
sham.vs.tac14.de.cl4 <- sham.vs.tac14.de[sham.vs.tac14.de$cluster==0,] #pick cluster 0 - peri1
sham.vs.tac14.de.cl4$high <-""
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,p_val_adj)
sham.vs.tac14.de.cl4$high[c(1:20)] <-sham.vs.tac14.de.cl4$gene_name[c(1:20)]
sham.vs.tac14.de.cl4$unique<-factor(sham.vs.tac14.de.cl4$unique,levels = c("shared","unique"))
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,unique)

ggplot(sham.vs.tac14.de.cl4) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj),colour=unique),size=2) +
  ggtitle(paste0("Peri1 DEG Sham vs TAC14")) +
  xlab("avg_logFC") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept=(mean(-log10(sham.vs.tac14.de.cl4$p_val_adj))+1), linetype="dashed", color = "grey")+
  geom_vline(xintercept=0.5, linetype="dashed", color = "grey")+
  geom_vline(xintercept=-0.5, linetype="dashed", color = "grey")+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = high),size = 5,min.segment.length = 0) +
  scale_color_manual(values = c("shared"="blue","unique"="red")) +
  theme(legend.title = element_blank(),legend.position = c(0.95,0.1),plot.title = element_text(size = 20, hjust = 0.5),axis.title = element_text(size = 20),axis.text = element_text(size=15))  +
  NoLegend()
ggsave(filename = paste0(sample,"volcano_plot_Peri1 DEG Sham vs TAC14.jpeg"), width=8, height =8)
ggsave(filename = paste0(sample,"volcano_plot_Peri1 DEG Sham vs TAC14.svg"), width=8, height =8)

#pericyte2
sham.vs.tac14.de.cl4 <- sham.vs.tac14.de[sham.vs.tac14.de$cluster==1,] #pick cluster 0 - peri2
sham.vs.tac14.de.cl4$high <-""
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,p_val_adj)
sham.vs.tac14.de.cl4$high[c(1:20)] <-sham.vs.tac14.de.cl4$gene_name[c(1:20)]
sham.vs.tac14.de.cl4$unique<-factor(sham.vs.tac14.de.cl4$unique,levels = c("shared","unique"))
sham.vs.tac14.de.cl4 <- arrange(sham.vs.tac14.de.cl4,unique)

ggplot(sham.vs.tac14.de.cl4) +
  geom_point(aes(x=avg_logFC, y=-log10(p_val_adj),colour=unique),size=2) +
  ggtitle(paste0("Peri2 DEG Sham vs TAC14")) +
  xlab("avg_logFC") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept=(mean(-log10(sham.vs.tac14.de.cl4$p_val_adj))+1), linetype="dashed", color = "grey")+
  geom_vline(xintercept=0.5, linetype="dashed", color = "grey")+
  geom_vline(xintercept=-0.5, linetype="dashed", color = "grey")+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = high),size = 5,min.segment.length = 0) +
  scale_color_manual(values = c("shared"="blue","unique"="red")) +
  theme(legend.title = element_blank(),legend.position = c(0.95,0.1),plot.title = element_text(size = 20, hjust = 0.5),axis.title = element_text(size = 20),axis.text = element_text(size=15))  +
  NoLegend()
ggsave(filename = paste0(sample,"volcano_plot_Peri2 DEG Sham vs TAC14.jpeg"), width=8, height =8)
ggsave(filename = paste0(sample,"volcano_plot_Peri2 DEG Sham vs TAC14.svg"), width=8, height =8)

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
number.of.DE.df<-arrange(number.of.DE.df,cluster)
number.of.DE.df$cluster=factor(levels(Samples.combined$celltypes.short),levels=levels(Samples.combined$celltypes.short))

ggplot(number.of.DE.df, aes(x = cluster)) + 
  geom_bar(aes(y = DE.up,fill="up"), stat = "identity") + 
  geom_bar(aes(y = -DE.down,fill="down"), stat = "identity") +
  scale_fill_manual(values =c("#3086e1","#e14f30")) +
  scale_y_continuous(labels = c(rev(seq(0,max(number.of.DE.df$DE.down)+5,10)),seq(0,max(number.of.DE.df$DE.up)+5,10)),
                     breaks = c(rev(seq(0,max(number.of.DE.df$DE.down)+5,10))*-1,seq(0,max(number.of.DE.df$DE.up)+5,10)))+
  theme(legend.title = element_blank(),axis.title.x = element_blank(),axis.text = element_text(size=16),axis.text.x = element_text(angle = 45,hjust = 1),axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ylab("Number of DEG")
ggsave(filename = paste0(sample,"number of Up&Down DE genes per cluster.jpeg"), width=5, height = 5)
ggsave(filename = paste0(sample,"number of Up&Down DE genes per cluster.svg"), width=5, height = 5)


#----score for ECM gene sets-----------
Samples.combined <- scoreECM(Samples.combined)

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "Core_matrisome1",pt.size = 0,cols = stim.col,split.by = "stim")+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust = 1))
ggsave(filename = paste0(sample,"core_matrisome_split_vln.jpeg"), width=10, height = 10)

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "Collagens1",pt.size = 0,cols = stim.col,split.by = "stim")+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust = 1))
ggsave(filename = paste0(sample,"collagens_split_vln.jpeg"), width=10, height = 10)

FeaturePlot(Samples.combined, features = "Collagens1",pt.size = 2,order = T,label = FALSE)+ scale_colour_gradientn(colours =gradient.col[2:11]) +theme(legend.position = c(0.9,0.2))
ggsave(filename = paste0(sample,"collagens.jpeg"), width=10, height = 10)
ggsave(filename = paste0(sample,"collagens.svg"), width=10, height = 10)

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


#----score for Collagen gene sets-----------
Samples.combined <- scoreCollagens(Samples.combined)
lvl=NULL
for (levels in levels(Samples.combined$celltypes)) {lvl=c(lvl,paste0(levels,"_",levels(Samples.combined$stim)))}
Idents(Samples.combined)=factor(paste0(Samples.combined$celltypes,"_",Samples.combined$stim),levels = lvl)
DotPlot(Samples.combined,dot.scale = 10,features =  c("fibrillar.col1","network.forming.col1","FACIT.col1","transmembrane.col1","multiplexin.col1"))+coord_flip()+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(size=8),axis.text = element_text(size=12),legend.key.size = unit(0.2, 'cm'))+ scale_colour_gradientn(colours =gradient.col)
ggsave(filename = paste0(sample,"Collagen-subgroup_scoring_dp.jpeg"), width=8, height = 3.5)
ggsave(filename = paste0(sample,"Collagen-subgroup_scoring_dp.svg"), width=8, height = 3.5)

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
#reduce to relevant expressed for heatmap
col.genes.rel = c("Col12a1","Col14a1","Col15a1","Col16a1","Col1a1","Col1a2","Col27a1","Col3a1","Col4a1","Col4a2","Col4a3","Col4a4","Col4a5","Col5a1","Col5a2","Col5a3","Col6a1","Col6a2","Col6a3","Col6a6","Col8a1")
col.genes = col.genes[col.genes %in% col.genes.rel]
Idents(Samples.combined)="celltypes.short"
mural.sub = subset(Samples.combined,ident=c("Peri1","Peri2","VSMC1","VSMC2"))
Idents(mural.sub)=paste0(mural.sub$celltypes.short,"_",mural.sub$stim)
avgexp = AverageExpression(mural.sub, return.seurat = T)
expression.values <- FetchData(avgexp,vars = col.genes,slot = "scale.data" )
expression.values<-as.data.frame(t(expression.values))
expression.values<-expression.values[,c("Peri1_Sham","Peri1_TAC","Peri2_Sham","Peri2_TAC","VSMC2_Sham","VSMC2_TAC","VSMC1_Sham","VSMC1_TAC")]

expression.values.heatmap <- pheatmap(expression.values,cluster_cols = F,cluster_rows = T,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_expressed_collagens.jpeg"), width=700 , height =800,quality = 100,res = 200)
print(expression.values.heatmap)
dev.off()
svg(filename = paste0(sample,"_expressed_collagens.svg"))
print(expression.values.heatmap)
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

VlnPlot(Samples.combined,features = c("Mmp3","Mmp13","Mmp2","Mmp9","Mmp28","Mmp14"),stack = T,group.by = "celltypes.short",split.by = "stim",flip = T,cols = stim.col)
ggsave(filename = paste0(sample,"mmp_stack_vln.jpeg"), width=7, height = 7)

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
