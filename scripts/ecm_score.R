#script for performing ECM scoring based on a script my mahmoud
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(readxl)
library(dplyr)
#################################################################################################################
#############---code in the following lines need to be adjusted before running---################################
#################################################################################################################
#the input seurat object needs to be an integrated dataset with clustering done.
#the experimental conditions need to be saved unter $stim

#set working directory to the output folder of your data
setwd("~/Documents/FP_scRNA/R stuff/")

#enter name of the integrated dataset (e.g. integration_gli1_Sham_TAC) 
sample="Fibroblast_Col1a1_Gli1_Pdgfrb_integrated_"

#in case of rename of stim is requiered (e.g. here FP21_Cdh5_Sham is changed to just Sham)
Samples.combined$stim[Samples.combined$stim=="FP21_Cdh5_Sham"] <- "Sham"
Samples.combined$stim[Samples.combined$stim=="FP20_Cdh5_TAC"] <- "TAC"

#load the integrated seurat object to the variable subset.cluster
Samples.combined <- readRDS(file = "C:/Users/fpeis/Desktop/R Stuff/maurice data rat/ XXXXX.rds") #if its an rds file

#################################################################################################################
#############---------------------------------end--------------------------------################################
######################----run everything below here in one go--##################################################

# Process mouse matrisome genes from #http://matrisomeproject.mit.edu/
matrisome_mm_masterlist <- as.data.frame(read_excel("~/Documents/FP_scRNA/R stuff/references/matrisome_mm_masterlist.xls"))
matrisome_mm_masterlist<-matrisome_mm_masterlist[c("Division","Category","Gene Symbol")]
matrisome_mm_masterlist$Division=gsub(pattern = " ",replacement = "_",x = matrisome_mm_masterlist$Division)
matrisome_mm_masterlist$Division=gsub(pattern = "-",replacement = "_",x = matrisome_mm_masterlist$Division)
matrisome_mm_masterlist$Category=gsub(pattern = " ",replacement = "_",x = matrisome_mm_masterlist$Category)
matrisome_mm_masterlist$Category=gsub(pattern = "-",replacement = "_",x = matrisome_mm_masterlist$Category)
matrisome_mm_masterlist.2<-matrisome_mm_masterlist
matrisome_mm_masterlist.2$Category[matrisome_mm_masterlist$Division=="Matrisome_associated"] = "Matrisome_associated"
matrisome_mm_masterlist.2$Category[matrisome_mm_masterlist$Division=="Core_matrisome"] = "Core_matrisome"
matrisome_mm_masterlist<-rbind(matrisome_mm_masterlist,matrisome_mm_masterlist.2)
matrisome_mm_masterlist<-matrisome_mm_masterlist[matrisome_mm_masterlist$Division!="Retired",]
rm(matrisome_mm_masterlist.2)
matrisome_mm_genesetlist = list()
for (geneset in unique(matrisome_mm_masterlist$Category)) {
  matrisome_mm_genesetlist[[geneset]] = matrisome_mm_masterlist$`Gene Symbol`[matrisome_mm_masterlist$Category==geneset]
}

# ECM scoring for the seurat object
DefaultAssay(Samples.combined)="RNA"
ctrl_genes = 35 #important
dir.create("ECM scoring output")
dir.create("ECM scoring output/svg")
pdf(paste0("ECM scoring output/",sample,"_ECM-scoring summary.pdf"))
for (gset in names(matrisome_mm_genesetlist)){
  features = matrisome_mm_genesetlist[gset]
  message(gset)
    # Add average expression of genes in gset minus the average expression of
    # 35 random genes 
    Samples.combined = AddModuleScore(object = Samples.combined, features = features, name = gset, ctrl = ctrl_genes)
    print(FeaturePlot(Samples.combined, features = paste0(gset, '1'),pt.size = 1.5,order = T,label = FALSE) + 
      scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', 
      midpoint = 0, limits = c(-1,1)))
    ggsave(filename = paste0("ECM scoring output/",sample,"_",gset,"_f.jpeg"), width=10 , height = 10)
    ggsave(filename = paste0("ECM scoring output/svg/",sample,"_",gset,"_f.svg"), width=10 , height = 10)
    
    print(VlnPlot(Samples.combined, features = paste0(gset, '1'), pt.size = 0,split.by = "stim",split.plot= T,cols = c('#00C5C0','#FA756C')) +
      geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
      theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust = 0.5)))
    ggsave(filename = paste0("ECM scoring output/",sample,"_",gset,"_vln.jpeg"), width=10 , height = 10)
    ggsave(filename = paste0("ECM scoring output/svg/",sample,"_",gset,"_vln.svg"), width=10 , height = 10)
}

#plot ECM between subclusters
print(VlnPlot(Samples.combined, features = "Core_matrisome1", pt.size = 0,group.by = "seurat_clusters") +
  geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust = 0.5))+ NoLegend())
ggsave(filename = paste0("ECM scoring output/",sample,"_ECMscore_seurat_clusters_vln.jpeg"), width=10 , height = 10) 
ggsave(filename = paste0("ECM scoring output/svg/",sample,"_ECMscore_seurat_clusters_vln.svg"), width=10 , height = 10) 

dev.off()