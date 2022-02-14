#script for pathway activity prediction based on progeny
#the following tutorial was adopted: http://bioconductor.org/packages/release/bioc/vignettes/progeny/inst/doc/ProgenySingleCell.html

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
  
## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
Idents(Samples.combined)<-Samples.combined$celltypes.short
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
                            CellType = as.character(paste0(Samples.combined$celltypes.stim)),
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
progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))

#flip and reorder dataframe for visualization
summarized_progeny_scores_df <-t(summarized_progeny_scores_df)[,levels(Samples.combined$celltypes.stim)]

#for specific heatmap plots, part of the matrix were subseted

progeny_hmap = pheatmap(summarized_progeny_scores_df,fontsize=14,
                        fontsize_row = 10,
                        color=col.ramp(100), breaks = progenyBreaks,
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA,cluster_cols = F,cluster_rows = T)
jpeg(filename = paste0(sample,"_TF_PROGENy_heatmap.jpeg"), width=1400 , height =1200,quality = 100,res = 200)
progeny_hmap
dev.off()
svg(filename = paste0(sample,"_TF_PROGENy_heatmap.svg"),width = 7,height = 6)
progeny_hmap
dev.off()

