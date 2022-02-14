#script for transcription factor activity prediction
#source for parts of the code: https://saezlab.github.io/dorothea/articles/single_cell_vignette.html
library(dorothea)
library(viper)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(readr)


## We read Dorothea Regulons for mouse:
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores 
Samples.combined <- run_viper(Samples.combined, regulon, assay_key = "RNA",
                  options = list(method = "scale", minsize = 4,
                                 eset.filter = FALSE, cores = 15,
                                 verbose = FALSE))

## We compute the Nearest Neighbours to perform cluster
DefaultAssay(Samples.combined) <- "dorothea"
Samples.combined <- ScaleData(Samples.combined)


#---------------------------# We transform Viper scores, scaled by seurat, into a data frame to better handling the results--------------
viper_scores_df <- GetAssayData(Samples.combined, slot = "scale.data",
                                assay = "dorothea") %>%
  data.frame() %>%
  t()

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = gsub("-",".",names(Idents(Samples.combined))),
                            cell_type = paste0(Samples.combined$celltypes.short,"_",Samples.combined$stim),
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

#optional step to exclude clusters with extreme score
  #summarized_viper_scores = summarized_viper_scores[!(summarized_viper_scores$cell_type %in% c("RepEC_TAC","RepEC_Sham","RepEC_TAC_28","CyclEC_TAC","CyclEC_Sham","CyclEC_TAC_28","EndoEC_TAC","EndoEC_Sham","EndoEC_TAC_28","IntEC_TAC","IntEC_Sham","IntEC_TAC_28","StrEC_TAC","StrEC_Sham","StrEC_TAC_28","LymEC_TAC","LymEC_Sham","LymEC_TAC_28")),]
  #summarized_viper_scores = summarized_viper_scores[!(summarized_viper_scores$cell_type %in% c("Interferon Fibroblast_Sham","Interferon Fibroblast_TAC","Interferon Fibroblast_TAC_28","Stressed Fibroblast_Sham","Stressed Fibroblast_TAC","Stressed Fibroblast_TAC_28")),]

## We select highly variable TFs
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n((30*(2*length(levels(factor(Samples.combined$seurat_clusters))))), var) %>%
  distinct(tf)

#optimal: select the top10
#highly_variable_tfs <- highly_variable_tfs[c(1:10),]

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)


palette_length = 100
my_breaks <- c(seq(min(summarized_viper_scores_df), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length,
                   max(summarized_viper_scores_df),
                   length.out=floor(palette_length/2)))

#reorder for fibroblasts subclusters
#summarized_viper_scores_df <- summarized_viper_scores_df[c(4,5,6,16,17,18,7,8,9,10,11,12,1,2,3,13,14,15),]
#reorder for endothelial subclusters and #top5 up and down for EC
#summarized_viper_scores_df <- summarized_viper_scores_df[c(10,11,12,7,8,9,13,14,15,1,2,3,4,5,6),]
#summarized_viper_scores_df <- summarized_viper_scores_df[,c("Snai2","Bhlhe40","Mnt","Zbtb7a","Zfp263","Smad3","Hif1a","Klf5","Klf9","Rest","Myc")]

#for EC overall plot
# summarized_viper_scores_df <- summarized_viper_scores_df[c(10,11,12,
#                                                            7,8,9,
#                                                            13,14,15,
#                                                            31,32,33,
#                                                            22,23,24,
#                                                            1,2,3,
#                                                            4,5,6,
#                                                            25,26,27,
#                                                            16,17,18,
#                                                            28,29,30,
#                                                            19,20,21),]



viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14,
                       fontsize_row = 10,
                       color=col.ramp(100), breaks = my_breaks,
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA,cluster_cols = F,cluster_rows = F)
jpeg(filename = paste0(sample,"_Dorothea_TF_heatmap_full.jpeg"), width=2000 , height =1000,quality = 100,res = 200)
viper_hmap
dev.off()
svg(filename = paste0(sample,"_Dorothea_TF_heatmap_full.svg"),width=10 , height =5)
viper_hmap
dev.off()

#--------ploting interessting TF activity scores for fibroblasts------------------
FeaturePlot(Samples.combined,features = "Tead1",order = T,pt.size = 0.5) + scale_colour_gradientn(colours = gradient.col) + ggtitle("Tead1 activity score")
ggsave(filename = paste0(sample,"_Tead1_TF_activity.jpeg"), width=10, height = 10)

VlnPlot(Samples.combined,group.by = "celltypes.short",features = "Tead1",pt.size = 0,cols = cluster.col)+
  geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.shape = NA,coef=0,lwd=0.3) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45,hjust =1,size = 15),axis.text.y = element_text(size = 10))+NoLegend()+ ggtitle("Tead1 activity score") + ylab("Activity Score")
ggsave(filename = paste0(sample,"_Tead1_TF_activity_vln.jpeg"), width=5, height = 5)
ggsave(filename = paste0(sample,"_Tead1_TF_activity_vln.svg"), width=5, height = 5) 

tead1.score.data <- FetchData(Samples.combined,vars = c("Tead1","celltypes"))
#test overall
k.result <- kruskal.test(tead1.score.data$Tead1~tead1.score.data$celltypes)
#post-hoc test individual comparisons
pair.result <- pairwise.wilcox.test(tead1.score.data$Tead1,tead1.score.data$celltypes,paired = F,p.adjust.method = "bonferroni")
pair.result.df <- as.data.frame(pair.result$p.value)
write_xlsx(pair.result.df,path = paste0(sample,"_test result Tead1score.xlsx"))

