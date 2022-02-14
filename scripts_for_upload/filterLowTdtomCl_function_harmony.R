library(writexl)
library(Seurat)
library(cowplot)
library(ggplot2)
library(plyr)
library(genesorteR, quietly = TRUE)
library(tidyverse)

filterLowTdtomCl.harmony <- function(object,samplename){
  
  #---exclude cells from a genotype if they contribute to a cluster with less than 20 cells
  object@meta.data$joinedident = paste0(object@meta.data$orig.ident,object@meta.data$seurat_clusters)
  cluster.genotype = as.data.frame(table(object@meta.data$joinedident))
  low20 = cluster.genotype$Freq<20
  excludeL20=cluster.genotype$Var1[low20]
  write_xlsx(as.data.frame(excludeL20),path = paste0(samplename, "_clustercontribution less than 20 cells.xlsx"))
  Idents(object)=object@meta.data$joinedident
  object=subset(object,idents = excludeL20,invert=T)
  
  #----filter out cluster contributions of genotypes lower than geneProb of 80%------
  pdf(file = paste0(samplename, '_tdTomato_filtering_before-after.pdf'))
  #plot unfilered
  Idents(object)=object@meta.data$seurat_clusters
  print(VlnPlot(object, features = c("tdTomato-Green"),assay = "RNA", pt.size = 0))
 # ggsave(filename = paste0(samplename,"_tdTom_unfiltered_Vln.jpeg"), width=10 , height =10,limitsize = F)
  Idents(object)=object@meta.data$orig.ident
  print(RidgePlot(object,features = "tdTomato-Green",slot = "counts",assay = "RNA",log = T))
#  ggsave(filename = paste0(samplename,"_before_0tdTom_exclude.jpeg"), width=10 , height = 10)
  dev.off()
  
  #filter out
  all.geno=levels(factor(object@meta.data$orig.ident))
  for (geno in all.geno) {
    message(paste0("subsetting for tdTomlow filtering: ",geno))
    Idents(object)=object@meta.data$orig.ident
    object.subset = subset(object, idents = geno)
    Idents(object.subset)=object.subset@meta.data$seurat_clusters
    sg = sortGenes(object.subset@assays$RNA@data, Idents(object.subset))
    colnames(sg$condGeneProb) = paste0(levels(Idents(object.subset)))
    object.subset_geneProb = as.data.frame(sg$condGeneProb)["tdTomato-Green",]
    cl.to.exclude = ncol(object.subset_geneProb %>% select(where(~ .x < 0.8)))
    if (cl.to.exclude != 0) {
      object.subset_delet = as.vector(paste0(geno,colnames(object.subset_geneProb %>% select(where(~ .x < 0.8)))))
      message(paste0(" to delet: ",object.subset_delet))
      Idents(object)=object@meta.data$joinedident
      object = subset(object, idents = object.subset_delet , invert=T)  
      }
  }
  
  #-----exclude tdtom 0 for final reclustering--------------------------------------------
  DefaultAssay(object) <- "RNA"
  object=subset(object,cells = WhichCells(object,slot = "counts",expression = `tdTomato-Green` == 0),invert=T)
  
  Idents(object)=object@meta.data$seurat_clusters
  return(object)
}
