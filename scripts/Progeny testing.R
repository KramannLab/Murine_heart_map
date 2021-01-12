#http://bioconductor.org/packages/release/bioc/vignettes/progeny/inst/doc/ProgenySingleCell.html

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
  
#################################################################################################################
#############---code in the following lines need to be adjusted before running---################################
#################################################################################################################
#the input seurat object needs to be an integrated dataset with clustering done.
#the experimental conditions need to be saved unter $stim

#set working directory to the output folder of your data
setwd("C:/Users/fpeis/Desktop/R Stuff/annotation for heart perivascular map/analysis on omen/cdh5 pairwise/gProfileR/")

#enter name of the integrated dataset (e.g. integration_gli1_Sham_TAC) 
sample ="Integration_rat_control_vs_stent"
control_sample ="Sham" #same as saved in $stim
treatment_sample ="TAC" #same as saved in $stim

#in case of rename of stim is requiered (e.g. here FP21_Cdh5_Sham is changed to just Sham)
Samples.combined$stim[Samples.combined$stim=="FP21_Cdh5_Sham"] <- "Sham"
Samples.combined$stim[Samples.combined$stim=="FP20_Cdh5_TAC"] <- "TAC"

#load the integrated seurat object to the variable subset.cluster
Samples.combined <- readRDS(file = "C:/Users/fpeis/Desktop/R Stuff/maurice data rat/ XXXXX.rds") #if its an rds file


#################################################################################################################
#############---------------------------------end--------------------------------################################
######################----run everything below here in one go--##################################################
## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
Idents(Samples.combined)<-Samples.combined$seurat_clusters
Samples.combined <- progeny(Samples.combined, scale=FALSE, organism="Mouse", top=500, perm=1,
                return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
Samples.combined <- Seurat::ScaleData(Samples.combined, assay = "progeny")

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <-
  as.data.frame(t(GetAssayData(Samples.combined, slot = "scale.data",
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

#here we use seurat_clusters combined with the condition as the cell tyes
CellsClusters <- data.frame(Cell = names(Idents(Samples.combined)), 
                            CellType = as.character(paste0(Idents(Samples.combined),"_",Samples.combined$stim)),
                            stringsAsFactors = FALSE)

## We match Progeny scores with the cell clusters.
#progeny_scores_df$Cell=gsub("-",".",progeny_scores_df$Cell)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

####################We plot the different pathway activities for the different cell populations

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))

#flip and reorder dataframe for visualization
summarized_progeny_scores_df=as.data.frame(t(summarized_progeny_scores_df[,-1]))
o1 <- order(factor(names(summarized_progeny_scores_df),levels=(c(paste0(levels(factor(Samples.combined$seurat_clusters)),"_",control_sample),paste0(levels(factor(Samples.combined$seurat_clusters)),"_",treatment_sample)))))
summarized_progeny_scores_df <- summarized_progeny_scores_df[ , o1]
summarized_progeny_scores_df <- summarized_progeny_scores_df[ , order(extract_numeric(names(summarized_progeny_scores_df)))]


progeny_hmap = pheatmap(summarized_progeny_scores_df,fontsize=14,
                        fontsize_row = 10,
                        color=myColor, breaks = progenyBreaks,
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA,cluster_cols = F)
jpeg(filename = paste0(sample,"_TF_PROGENy_heatmap.jpeg"), width=1600 , height =1200,quality = 100,res = 200)
progeny_hmap
dev.off()
svg(filename = paste0(sample,"_TF_PROGENy_heatmap.svg"))
progeny_hmap
dev.off()

