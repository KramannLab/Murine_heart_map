##------integration of fibroblasts from Pairwise filered Ng2, Myh11, Pdgfrb------------------------
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
setwd("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/cardiomyo_EC_excluded/")
sample="VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated_"
dir.create("VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated")

#colors for genotypes
col.pdgfrb="#EC407A"
col.myh11="#FFEB3B"
col.ng2="#795548"
color.cl = c(col.pdgfrb,col.myh11 ,col.ng2)

#load processed datafiles from pairwise
myh11.sham = readRDS(file = "~/Documents/FP_scRNA/R stuff/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/RDS/Myh11_integrated__Sham_indiv_after_pairwise_integration.rds")
myh11.tac = readRDS(file = "~/Documents/FP_scRNA/R stuff/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/RDS/Myh11_integrated__TAC_indiv_after_pairwise_integration.rds")

ng2.sham = readRDS(file = "~/Documents/FP_scRNA/R stuff/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/RDS/Ng2_integrated__Sham_indiv_after_pairwise_integration.rds")
ng2.tac = readRDS(file = "~/Documents/FP_scRNA/R stuff/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/RDS/Ng2_integrated__TAC_indiv_after_pairwise_integration.rds")

pdgfrb.sham = readRDS(file = "~/Documents/FP_scRNA/R stuff/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/RDS/PDGFRb_integrated__Sham_indiv_after_pairwise_integration_ONLY_VSMC&Pericytes.rds")
pdgfrb.tac = readRDS(file = "~/Documents/FP_scRNA/R stuff/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/RDS/PDGFRb_integrated__TAC_indiv_after_pairwise_integration_ONLY_VSMC&Pericytes.rds")

#----function to perfrom clustering after every filterstep------
recluster <- function(object,res){
  DefaultAssay(object) <- "integrated"
  # Run the standard workflow for visualization and clustering
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, verbose = FALSE)
  # t-SNE and Clustering
  object <- RunUMAP(object, reduction = "pca", dims = 1:20)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object, resolution = res)
  return(object)
}
####################################################---integration--#################################################################
Samples.anchors <- FindIntegrationAnchors(object.list = list(myh11.sham,myh11.tac,ng2.sham,ng2.tac,pdgfrb.tac,pdgfrb.sham), dims = 1:20)

#--------------------exclude tdTom from the anchorlist-------------------------
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list

#only keep variable features after integration
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined,file = paste0(sample,"after_integration.rds"))
Samples.combined <- readRDS("~/Documents/FP_scRNA/R stuff/final_seurat_analysis/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/RDS/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated_after_integration.rds") 

################################################----start of sample processing before analysis----##################################################

# subset the individual stim with just TAC or Sham
Samples.filtered=Samples.combined
TACs = Samples.filtered@meta.data$stim %in% c("FP22_NG2_TAC","FP14_PDGFRb_TAC","FP57_Myh11_TAC")
Samples.filtered@meta.data$stim[!TACs]="Sham"
Samples.filtered@meta.data$stim[TACs]="TAC" 
RidgePlot(Samples.filtered, assay = "RNA", features = "tdTomato-Green",group.by = "stim",log = T) + NoLegend()
Samples.combined=Samples.filtered
remove(Samples.filtered)

# add colum for only the genotype
Samples.filtered=Samples.combined
PDGFRb.ls = Samples.filtered@meta.data$orig.ident %in% c("FP14_PDGFRb_TAC","FP3_PDGFRb_Sham")
NG2.ls = Samples.filtered@meta.data$orig.ident %in% c("FP22_NG2_TAC","FP23_NG2_Sham")
Myh11.ls = Samples.filtered@meta.data$orig.ident %in% c("FP56_Myh11_Sham","FP57_Myh11_TAC")
Samples.filtered@meta.data$orig.geno[PDGFRb.ls]="PDGFRb"
Samples.filtered@meta.data$orig.geno[NG2.ls]="NG2"
Samples.filtered@meta.data$orig.geno[Myh11.ls]="Myh11"
Samples.combined=Samples.filtered
remove(Samples.filtered)
#settle the order
Samples.combined$orig.ident = factor(Samples.combined$orig.ident,levels = c("FP3_PDGFRb_Sham","FP14_PDGFRb_TAC","FP56_Myh11_Sham","FP57_Myh11_TAC","FP22_NG2_TAC","FP23_NG2_Sham"))
Samples.combined$orig.geno = factor(Samples.combined$orig.geno, levels = c("PDGFRb","Myh11","NG2"))

##-------first clustering after integration-------------------------------------------
Samples.combined = recluster(Samples.combined,0.6)

#plot for quality check
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 2, cols = c('#FA756C80','#00C5C080')) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"quality_check_clusters_res0.6.jpeg"), width=20 , height = 10)

#-------------ploting only tdtom
#quality check on tdTom expression
FeaturePlot(Samples.combined, features = c("tdTomato-Green"),pt.size = 1,order = T)
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

#--------after analysing the clusters based on steps until here cl9 (res0.6) was found to be cardiomyocytes, which is very likly to be a contamination
#--------also cl7 (res0.6) was found to be endothelial cell, which is very likly to be a contamination
#therefore cl9 and 7 were removed and the dataset recluster + analysed
Samples.combined=subset(Samples.combined,idents = c("7","9"),invert=T)
Samples.combined = recluster(Samples.combined,0.6) #-increased resolution after checking clustertree

#---------------clustered and analysed the data further downstream. in cluster 4 "neuronal" a larger pool of very likely contamination from doublets detected
#---------------subset cluster4 "neuronal" to check for subtypes
Idents(Samples.combined)
neuronal.subset = subset(Samples.combined,idents = 4)
neuronal.subset<-recluster(neuronal.subset,0.4)#checked on clustree
DimPlot(neuronal.subset, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"subset_cl4_reclustered.jpeg"), width=5 , height = 5)
DefaultAssay(neuronal.subset)="RNA"
VlnPlot(neuronal.subset,features = c("Plp1","Kcna1"),ncol = 1)
ggsave(filename = paste0(sample,"subset_cl4_reclustered_schwannCell_markers.jpeg"), width=10 , height = 10)
#--------top10 markers per cluster----------------------
all.markers = FindAllMarkers(neuronal.subset, test.use = "MAST",min.pct = 0.3,assay = "RNA")
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene)
top10 <-unique(top10)
#vertikal
DotPlot(neuronal.subset, features = top10,dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"subset_cl4_reclustered_markers.jpeg"), width=5 , height = 12)
#very likely doublets detected in the neuronal cluster -> isolating those cells for subseting and reanalysis overall
cont.cells <- WhichCells(neuronal.subset,idents = c(0,1,4,5))
Samples.combined <- subset(Samples.combined,cells = cont.cells,invert=T)
Samples.combined <- recluster(Samples.combined,0.6)

#----renaming cl 3 and 4 for better visualisation
Samples.combined <- RenameIdents(object = Samples.combined, `2` = "sP")
Samples.combined <- RenameIdents(object = Samples.combined, `4` = "2")
Samples.combined <- RenameIdents(object = Samples.combined, `sP` = "4")
Samples.combined$seurat_clusters = Idents(Samples.combined)
Samples.combined$seurat_clusters = factor(Samples.combined$seurat_clusters,levels=c(0,1,2,3,4,5,6,7))
Idents(Samples.combined) = Samples.combined$seurat_clusters
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined$orig.geno = factor(Samples.combined$orig.geno,levels = c("PDGFRb","Myh11","NG2"))
#reploted quality plots after optimization and number changes

#-------normalize and scale RNA data------------------------
DefaultAssay(Samples.combined) <- "RNA"
Samples.combined <- NormalizeData(Samples.combined, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE)

################################################----end of sample processing before analysis----##################################################

########################################################----start of sample analysis----###########################################################

#clustertree for optimazation of cluster resolution
sc.int=Samples.combined
DefaultAssay(sc.int)<-"integrated"
sc.int = FindClusters(sc.int, resolution = seq(from=0.1, to=1.3, by=0.1), verbose = FALSE)
clustree(sc.int)
remove(sc.int)
ggsave(filename = paste0(sample,"_cluster_resolution_tree.jpeg"), width=10 , height = 10)

#cluster plots - seurat cluster
DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 1, label.size = 12) + NoLegend()
ggsave(filename = paste0(sample,"_only_clusters_res0.6.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_only_clusters_res0.6.svg"), width=10 , height = 10)
DimPlot(Samples.combined, reduction = "umap", label = F, pt.size = 1, label.size = 12) + NoLegend() + xlim(-12,6.5)
ggsave(filename = paste0(sample,"_only_clusters_res0.6_nolab.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_only_clusters_res0.6_nolab.svg"), width=10 , height = 10)
#cluster plots - sham vs. tac
DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 1, cols = c('#00C5C080','#FA756C80'))+NoLegend()
ggsave(filename = paste0(sample,"_shamvstac_clusters.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_shamvstac_clusters.svg"), width=10 , height = 10)
#cluster plots - genotype
DimPlot(Samples.combined, reduction = "umap", split.by = "orig.geno",group.by = "orig.geno", pt.size = 1,cols = color.cl) + NoLegend()
ggsave(filename = paste0(sample,"_genotype_clusters.jpeg"), width=30 , height = 10)
ggsave(filename = paste0(sample,"_genotype_clusters.svg"), width=30 , height = 10)
DimPlot(Samples.combined, reduction = "umap",group.by = "orig.geno", pt.size = 1,cols = color.cl) + theme(legend.position = c(0.8,0.1))
ggsave(filename = paste0(sample,"_genotype_clusters_res0.6.jpeg"), width=10 , height = 10)
ggsave(filename = paste0(sample,"_genotype_clusters_res0.6.svg"), width=10 , height = 10)
#clusters + conditions
p1 <- DimPlot(Samples.combined, reduction = "umap", label = TRUE, pt.size = 2, label.size = 12) + NoLegend()
p2 <- DimPlot(Samples.combined, reduction = "umap", group.by = "stim", pt.size = 2, cols = c('#00C5C080','#FA756C80')) +theme(legend.position = "bottom")
plot_grid(p1, p2,ncol = 2)
ggsave(filename = paste0(sample,"clusters&conditions_res0.6.jpeg"), width=20 , height = 10)
ggsave(filename = paste0(sample,"clusters&conditions_res0.6.svg"), width=20 , height = 10)

#--------top5 markers per cluster----------------------
all.markers = FindAllMarkers(Samples.combined, test.use = "MAST",min.pct = 0.3,assay = "RNA")
#result of findallmarkers
saveRDS(all.markers,file = "all.markers_markergenes_per_cluster_result.RDS")
all.markers <-readRDS(file = "~/Documents/FP_scRNA/R stuff/final_seurat_analysis/VSMC&Pericytes_Myh11_Ng2_Pdgfrb_integrated/cardiomyo_EC_excluded/all.markers_markergenes_per_cluster_result.RDS")
top5 <- unique(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene))
#vertikal
DotPlot(Samples.combined, features = top5,dot.scale = 10,assay = "RNA")+ coord_flip()
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.jpeg"), width=5.7 , height = 12)
ggsave(filename = paste0(sample,"_top5_markers_per_cluster_v.svg"), width=5,7 , height = 12)
#horizontal
Idents(Samples.combined)= factor(Samples.combined@meta.data$seurat_clusters,levels(Idents(Samples.combined))[as.numeric(rev(levels(Idents(Samples.combined))))+1])
Samples.combined$sub_cell_type_short = factor(Samples.combined$sub_cell_type_short,levels = rev(sub_cell_type_short)) #see lower section for annotation
DotPlot(Samples.combined,group.by = "sub_cell_type_short",features = top5,dot.scale = 10,assay = "RNA")+theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title = element_blank())
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.jpeg"), width=14 , height = 3.4)
ggsave(filename = paste0(sample,"_top5_featureDOTplot_custom_cl_h.svg"), width=14 , height = 3.4)
Idents(Samples.combined)=Samples.combined@meta.data$seurat_clusters

#--------ploting of multiple featues and save individual------------------------------
features_VSMC = c("Myh11","Acta2","Tagln","Cnn1","Myocd","Vcl","Smtn","Ly6a","Myh10")
features_peri = c("Cspg4","Pdgfrb","Pdgfra","Colec11","Kcnj8","Abcc9","Vtn","Notch3")
features_fibro = c("Dcn","Lum","Cxcl1","Col1a1","Apoe","Ogn","Postn")
features_endo = c("Cdh5","Pecam1","Kdr","Fabp4","Cd36","Srgn")
features_all = c(features_endo,features_fibro,features_peri,features_VSMC,"end")
allfp = vector("list",1)
dir.create("multiple_marker_plots")
for (f in features_all) {
  if (f=="end") {
    print(ggarrange(plotlist = allfp, widths = 30, heights = 30))
    ggsave(filename = paste0("multiple_marker_plots/",sample,"_all_features.jpeg"), width=30 , height =30)
  }
  FeaturePlot(Samples.combined, features = f, min.cutoff = "q9",pt.size = 1) 
  ggsave(filename = paste0("multiple_marker_plots/",sample,"_",f,".jpeg"), width=10 , height = 10)
  allfp[[f]]<- FeaturePlot(Samples.combined, features = f, min.cutoff = "q9",pt.size = 1)
  
  VlnPlot(Samples.combined, features = f)
  ggsave(filename = paste0("multiple_marker_plots/",sample,"_",f,"_Vln.jpeg"), width=10 , height = 10)
}

#-----bar graph per genotype
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
ggplot(data=breakdown.df, aes(x=Var1, y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = levels(factor(breakdown.df$Var1)),name="Cluster")+
  scale_y_continuous(name="Normalizes percentage of contribution [%]") +
  scale_fill_manual( values = c("Myh11" = col.myh11, "NG2" = col.ng2, "PDGFRb" = col.pdgfrb))+
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
breakdown.df$Var2 = factor(breakdown.df$Var2, levels = c("TAC","Sham"))
sub_cell_type_short = c("Peri1","Peri2","VSMC1","VSMC2","StrPeri","StrVSMC","Neu","IntPeri")

ggplot(data=breakdown.df, aes(x=rev(Var1), y=value,fill=Var2)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=c(sub_cell_type_short,sub_cell_type_short) ,y=1),hjust = 0, color="black", size=10)+
  coord_flip() + 
  scale_fill_manual(values = c("Sham" = "#00C5C0", "TAC" = "#FA756C"))+
  scale_x_continuous(breaks = as.numeric(levels(factor(breakdown.df$Var1))),labels = rev(levels(factor(breakdown.df$Var1))),name="Cluster")+
  scale_y_continuous(name="Percentage of total cells [%]") +
  theme(axis.text.y = element_blank(),axis.ticks.y =element_blank() ,axis.title.y = element_blank(),axis.text.x = element_text(size=15),axis.title.x = element_text(size = 15)) + 
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_hline(yintercept=50, linetype="dashed", color = "black")+
  NoLegend()
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.jpeg"), width=5 , height = 10)
ggsave(filename = paste0(sample,"_barplot_percentage_per_stim and cl.svg"), width=5 , height = 10)


#---------------- add sub_cell_type according to annotation ---------------------------------
sub_cell_type = c("Pericyte 1", "Pericyte 2", "VSMC 1", "VSMC 2", "Stressed pericyte", "Stressed VSMC","Neuronal", "Interferon pericyte")
sub_cell_type_short = c("Peri1","Peri2","VSMC1","VSMC2","StrPeri","StrVSMC","Neu","IntPeri")
Samples.combined@meta.data$sub_cell_type=NULL
Samples.combined@meta.data$sub_cell_type_short=NULL
Samples.combined$sub_cell_type = Samples.combined$orig.ident
Samples.combined$sub_cell_type_short = Samples.combined$orig.ident

for (cluster.nr in levels(factor(Samples.combined$seurat_clusters))) {
  Samples.combined$sub_cell_type[Samples.combined$seurat_clusters==cluster.nr] <- sub_cell_type[as.numeric(cluster.nr)+1]
  Samples.combined$sub_cell_type_short[Samples.combined$seurat_clusters==cluster.nr] <- sub_cell_type_short[as.numeric(cluster.nr)+1]
}
Samples.combined$sub_cell_type = factor(Samples.combined$sub_cell_type,levels = sub_cell_type)
Samples.combined$sub_cell_type_short = factor(Samples.combined$sub_cell_type_short,levels = sub_cell_type_short)


###------sample specific plot for GOs-------------------------------
#selected two most interessting and representative GOs per cluster
go.picked = c("GO:0062023","REAC:R-MMU-216083","GO:0050839","GO:0016722","GO:0043292",	"GO:0005901",	"GO:0030029",	"REAC:R-MMU-445355",	"GO:0035976",	"GO:0010942",	"GO:0042552"	,"GO:0042063"	,"GO:0045087",	"GO:0035456" )

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
ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top2 representative GOs per cluster.jpeg"), width=5.2, height = 3.3)
ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/svg/",sample,"top2 representative GOs per cluster.svg"), width=5.2, height = 3.3)

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
  geom_text(aes(label=Samples.combined$Core_matrisome1_delta, y=c(rep(0.3,times=length(Samples.combined$orig.ident)))),hjust = -0.3, color="black", size=4)+
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