library(dorothea)
library(viper)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(readr)

#################################################################################################################
#############---code in the following lines need to be adjusted before running---################################
#################################################################################################################
#the input seurat object needs to be an integrated dataset with clustering done.
#the experimental conditions need to be saved unter $stim

#set working directory to the output folder of your data
setwd("C:/Users/fpeis/Desktop/R Stuff/")

#enter name of the integrated dataset (e.g. integration_gli1_Sham_TAC) 
sample ="Integration_rat_control_vs_stent"
control_sample ="Sham" #same as saved in $stim
treatment_sample ="TAC" #same as saved in $stim

#in case of rename of stim is requiered (e.g. here FP21_Cdh5_Sham is changed to just Sham)
Samples.combined$stim[Samples.combined$stim=="FP21_Cdh5_Sham"] <- "Sham"
Samples.combined$stim[Samples.combined$stim=="FP20_Cdh5_TAC"] <- "TAC"

#load the integrated seurat object to the variable subset.cluster
Samples.combined <- readRDS(file = "C:/Users/fpeis/Desktop/R Stuff/maurice data rat/ XXXXX.rds") #if its an rds file

#number of cored for viper (on windows error when >1)
cores=1 # for windows
cores=10 # for workstation

#in case some cluster have an overwhelming singal, subset
Idents(Samples.combined)="seurat_clusters"
Samples.combined = subset(Samples.combined,idents = c("3","5","7","8","9","10"),invert=T)
#################################################################################################################
#############---------------------------------end--------------------------------################################
######################----run everything below here in one go--##################################################
## We read Dorothea Regulons for mouse:
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores 
Samples.combined.vip <- run_viper(Samples.combined, regulon, assay_key = "RNA",
                  options = list(method = "scale", minsize = 4,
                                 eset.filter = FALSE, cores = cores,
                                 verbose = FALSE))

## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = Samples.combined.vip) <- "dorothea"
Samples.combined.vip <- ScaleData(Samples.combined.vip)
Samples.combined.vip <- RunPCA(Samples.combined.vip, features = rownames(Samples.combined.vip), verbose = FALSE)
Samples.combined.vip <- FindNeighbors(Samples.combined.vip, dims = 1:10, verbose = FALSE)
Samples.combined.vip <- FindClusters(Samples.combined.vip, resolution = 0.4, verbose = FALSE)

Samples.combined.vip <- RunUMAP(Samples.combined.vip, dims = 1:10, umap.method = "uwot", metric = "cosine")

DimPlot(Samples.combined.vip, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
ggsave(filename = paste0(sample,"_Dorothea_clustering_res0.4.jpeg"), width=10 , height =10)

#---------------------------# We transform Viper scores, scaled by seurat, into a data frame to better handling the results--------------
viper_scores_df <- GetAssayData(Samples.combined.vip, slot = "scale.data",
                                assay = "dorothea") %>%
  data.frame() %>%
  t()

## We create a data frame containing the cells and their clusters
Idents(Samples.combined)=Samples.combined@meta.data$stim
CellsClusters <- data.frame(cell = gsub("-",".",names(Idents(Samples.combined))),
                            cell_type = paste0(Samples.combined@meta.data$seurat_clusters,"_",Samples.combined@meta.data$stim),
                            stringsAsFactors = FALSE)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>%
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>%
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n((20*(2*length(levels(factor(Samples.combined$seurat_clusters))))), var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length,
                   max(summarized_viper_scores_df),
                   length.out=floor(palette_length/2)))

#flip and reorder dataframe for visualization
summarized_viper_scores_df=as.data.frame(t(summarized_viper_scores_df[,-1]))
o1 <- order(factor(names(summarized_viper_scores_df),levels=(c(paste0(levels(factor(Samples.combined$seurat_clusters)),"_",control_sample),paste0(levels(factor(Samples.combined$seurat_clusters)),"_",treatment_sample)))))
summarized_viper_scores_df <- summarized_viper_scores_df[ , o1]
summarized_viper_scores_df <- summarized_viper_scores_df[ , order(extract_numeric(names(summarized_viper_scores_df)))]


viper_hmap <- pheatmap(summarized_viper_scores_df,fontsize=14,
                       fontsize_row = 10,
                       color=my_color, breaks = my_breaks,
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA,cluster_cols = F)
jpeg(filename = paste0(sample,"_Dorothea_TF_heatmap_specific.jpeg"), width=1600 , height =1200,quality = 100,res = 200)
viper_hmap
dev.off()

##############################################################################################################################################
#loop for dorothea run sperate for each cluster
dir.create("dorothea_per_cluster")
for (cluster.nr in levels(factor(Samples.combined$seurat_clusters))) {
  Idents(Samples.combined)="seurat_clusters"
  #subest in for-loop for each cluster
  Samples.combined.subet = subset(Samples.combined,ident = cluster.nr)
  Samples.combined.subet$seurat_clusters = factor(Samples.combined.subet$seurat_clusters)
  
  ## We read Dorothea Regulons for mouse:
  dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))
  
  ## We obtain the regulons based on interactions with confidence level A, B and C
  regulon <- dorothea_regulon_mouse %>%
    dplyr::filter(confidence %in% c("A","B","C"))
  
  ## We compute Viper Scores 
  Samples.combined.vip <- run_viper(Samples.combined.subet, regulon,assay_key = "RNA",
                                    options = list(method = "scale", minsize = 4,
                                                   eset.filter = FALSE, cores = cores,
                                                   verbose = FALSE))
  
  ## We create the assay and scale the data
  DefaultAssay(object = Samples.combined.vip) <- "dorothea"
  Samples.combined.vip <- ScaleData(Samples.combined.vip)
  
  #---------------------------# We transform Viper scores, scaled by seurat, into a data frame to better handling the results--------------
  viper_scores_df <- GetAssayData(Samples.combined.vip, slot = "scale.data",
                                  assay = "dorothea") %>%
    data.frame() %>%
    t()
  
  ## We create a data frame containing the cells and their clusters
  Idents(Samples.combined)=Samples.combined.subet@meta.data$stim
  CellsClusters <- data.frame(cell = gsub("-",".",names(Idents(Samples.combined.subet))),
                              cell_type = paste0(Samples.combined.subet@meta.data$seurat_clusters,"_",Samples.combined.subet@meta.data$stim),
                              stringsAsFactors = FALSE)
  
  ## We create a data frame with the Viper score per cell and its clusters
  viper_scores_clusters <- viper_scores_df  %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)
  
  ## We summarize the Viper scores by cellpopulation
  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))
  
  ## We select the 20 most variable TFs. (20*9 populations = 180)
  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(2*20, var) %>%
    distinct(tf)
  
  ## We prepare the data for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "tf") %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(min(summarized_viper_scores_df), 0,
                     length.out=ceiling(palette_length/2) + 1),
                 seq(max(summarized_viper_scores_df)/palette_length,
                     max(summarized_viper_scores_df),
                     length.out=floor(palette_length/2)))
  
  #flip and reorder dataframe for visualization
  summarized_viper_scores_df=as.data.frame(t(summarized_viper_scores_df[,-1]))
  o1 <- order(factor(names(summarized_viper_scores_df),levels=(c(paste0(levels(factor(Samples.combined.subet$seurat_clusters)),"_",control_sample),paste0(levels(factor(Samples.combined.subet$seurat_clusters)),"_",treatment_sample)))))
  summarized_viper_scores_df <- summarized_viper_scores_df[ , o1]
  summarized_viper_scores_df <- summarized_viper_scores_df[ , order(extract_numeric(names(summarized_viper_scores_df)))]
  
  
  viper_hmap <- pheatmap(summarized_viper_scores_df,fontsize=14,
                         fontsize_row = 10,
                         color=my_color, breaks = my_breaks,
                         main = "DoRothEA (ABC)", angle_col = 0,
                         treeheight_col = 0,  border_color = NA,cluster_cols = F)
  jpeg(filename = paste0("dorothea_per_cluster/",sample,"_Dorothea_TF_heatmap_cluster_",cluster.nr,".jpeg"), width=800 , height =1200,quality = 100,res = 200)
  print(viper_hmap)
  dev.off()
  
}
############################################################################################################################
#------------ploting Dorothea between Sham and TAC
Idents(Samples.combined) = "stim"

average_stim= AverageExpression(Samples.combined,return.seurat = T,assays = "RNA",slot="data")

#regulon load above
average_stim <- run_viper(average_stim[["RNA"]]@data, regulon, assay_key = "RNA",
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = cores, 
                                           verbose = FALSE))
#get top10 TFs that are up and down
average_stim <- as.data.frame(average_stim)
average_stim$difference = ((sqrt(average_stim$Sham^2)+sqrt(average_stim$TAC^2))*average_stim$TAC)/sqrt(average_stim$TAC^2)
average_stim <- arrange(average_stim,-difference)
top10_tf_up_down <-rbind(head(average_stim,n=10),tail(average_stim,n=10))
top10_tf_up_down$difference=NULL

#ploting heatmap
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(top10_tf_up_down), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(top10_tf_up_down)/palette_length,
                   max(top10_tf_up_down),
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(top10_tf_up_down,fontsize=14,
                       fontsize_row = 10,
                       color=my_color, breaks = my_breaks,
                       main = "DoRothEA (ABC)", angle_col = 0,
                       treeheight_col = 0,  border_color = NA,cluster_cols = F)

jpeg(filename = paste0("dorothea_per_cluster/",sample,"_Dorothea_TF_heatmap_cluster_ShamVsTAC.jpeg"), width=800 , height =1200,quality = 100,res = 200)
print(viper_hmap)
dev.off()

# check for high TFs if they are subcluster specific by scoring there regulom from dorothea
regulon.df = as.data.frame(regulon) 
DefaultAssay(Samples.combined)="RNA"
Idents(Samples.combined)="seurat_clusters"
target.list=NULL
for (tf in rownames(top10_tf_up_down)) {
  target.list[[tf]] = regulon.df$target[regulon.df$tf==tf]
  Samples.combined<- AddModuleScore(Samples.combined,features = target.list[tf],name = tf, ctrl = 35)
  
  print(VlnPlot(Samples.combined, features = paste0(tf, '1'), pt.size = 0,split.by = "stim",split.plot= T,cols = c('#00C5C0','#FA756C')) +
          geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
          theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust = 0.5)))
  ggsave(filename = paste0("dorothea_per_cluster/",sample,"_TF_",tf,"_module_score.jpeg"), width=10 , height = 10)
}