#correlation analysis based on the core matrisome score from the matrisome database
#parts of the analysis strategy and the python script were adopted from Mayr et. al. 2021
#https://www.embopress.org/doi/full/10.15252/emmm.202012871
library(reticulate)
library(viper)
library(dorothea)

# python code for correlation analysis
source_python("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/Python-Correlation.py")

# Step 1: Calculate ECM Scores
Samples.combined <- scoreECM(Samples.combined)
ECM_Genes <- unlist(matrisome_mm_genesetlist) #load matrisome from the ECM score function

# Perform Correlation Analysis to find Top Correlating Genes with core matrisome Score

# Cluster at a high resolution for later Pseudobulking
#also split the high res cluster into 3 parts for the 3 different conditions (sham cells will always score lower then TAC and should not be combined)
DefaultAssay(Samples.combined) <- "integrated"
Samples.combined <- FindClusters(Samples.combined, resolution = 10)
DefaultAssay(Samples.combined) <- "RNA"
Idents(Samples.combined)<-paste0(Samples.combined$seurat_clusters,"_",Samples.combined$stim)
Samples.combined$cluster.stim = Idents(Samples.combined)
#calculate mean gene expression of high resolution clusters as well as mean ECM score. 
avgexp = AverageExpression(Samples.combined, return.seurat = T,assays = "RNA")
#remove ECM Genes for Correlation Analysis and get average gene expression per "mini"cluster
MeanExpr <-FetchData(avgexp,vars = rownames(Samples.combined@assays$RNA)[!(rownames(Samples.combined@assays$RNA)%in%ECM_Genes)],slot = "scale.data")
#get average ECM score per "mini" cluster 
ExprECM <-FetchData(Samples.combined,vars =c("Core_matrisome1","cluster.stim"))
MeanExprECM <- ExprECM %>%  group_by(cluster.stim) %>%  summarise(avg = mean(Core_matrisome1),std = sd(Core_matrisome1))
MeanExprECM <- as.data.frame(MeanExprECM)
rownames(MeanExprECM)<-MeanExprECM$cluster.stim
MeanExprECM$cluster.stim=NULL
MeanExprECM$std=NULL
colnames(MeanExprECM)<-"Core_matrisome_score"
#combine into one dataframe
CorMatrix=cbind(MeanExprECM, MeanExpr)
#correlation
cors=as.data.frame(py$correlate_means_to_gene(CorMatrix, corr_gene = "Core_matrisome_score"))

#plot heatmap of highest correlating genes sorted by ecm score TOP100
expression.values <-FetchData(avgexp,vars = rownames(cors)[2:100],slot = "scale.data")
expression.values=cbind(MeanExprECM, expression.values,)
expression.values=arrange(expression.values,Core_matrisome_score)
expression.values.heatmap <- pheatmap(t(expression.values),main = "Top100 genes correlating with core_matrisome_score",cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),labels_col = F)
jpeg(filename = paste0(sample,"_Top100 genes correlating with core_matrisome_score_res10_split.jpeg"), width=1200 , height =1200,quality = 100,res = 200)
print(expression.values.heatmap)
dev.off()

#test for GOs associated the top100 correlating genes
go.results <- gost(query = rownames(cors[cors$pvalue<0.01 &cors$spearman_corr>0,])[1:100], 
                   organism = "mmusculus", ordered_query = T, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                   measure_underrepresentation = FALSE, evcodes = FALSE, 
                   user_threshold = 0.05, correction_method = "g_SCS", 
                   domain_scope = "annotated", custom_bg = NULL, 
                   numeric_ns = "", sources = NULL, as_short_link = FALSE) 
go.results <- go.results$result
go.results <- go.results[go.results$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
go.results <- go.results[go.results$term_size<1000,] #filter out large and often very general GOs
go.results <- arrange(go.results,p_value)
go.results$term_name = factor(go.results$term_name,levels = rev(go.results$term_name))
ggplot(go.results[c(1:20),], aes(x=term_name, y="",size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
  geom_point() + 
  scale_size(range = c(2,6)) +
  coord_flip() +
  ylab("")+
  xlab("")+
  ggtitle(paste0("Top 20 GOs of top100 genes correlating to ECMscore"))+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle=45,hjust=1) ,legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5,size = 10),legend.key.size = unit(0.3, "cm"))
ggsave(filename = paste0(sample,"Top 20 GOs of top100 genes correlating to ECMscore.jpeg"), width=4, height = 4)
ggsave(filename = paste0(sample,"Top 20 GOs of top100 genes correlating to ECMscore.svg"), width=4, height = 4)



#--------correlate ECMscore to Dorothea TF acitivity scores----------------

### !!!!!!!!!!!!!!!!
## first perform ECMscoring and high resolution clustering as before (avg gene expression is not needed here)
### !!!!!!!!!!!!!!!!

## then run dorothea analysis with viper
#We read Dorothea Regulons for mouse:
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))
#We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))
#We compute Viper Scores 
Samples.combined <- run_viper(Samples.combined, regulon, assay_key = "RNA",
                                  options = list(method = "scale", minsize = 4,
                                                 eset.filter = FALSE, cores = 10,
                                                 verbose = FALSE))
#We add an assay and scale the scores
DefaultAssay(object = Samples.combined) <- "dorothea"
Samples.combined <- ScaleData(Samples.combined)
#We transform Viper scores, scaled by seurat, into a data frame to better handling the results
viper_scores_df <- GetAssayData(Samples.combined, slot = "scale.data",
                                assay = "dorothea") %>%
  data.frame() %>%
  t()
saveRDS(viper_scores_df,file = paste0(sample,"_viper_scores_df.RDS"))
#We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = gsub("-",".",names(Idents(Samples.combined))),
                            cell_type = Idents(Samples.combined),
                            stringsAsFactors = FALSE)
#We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>%
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)
#We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>%
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
#table of all tfs to include for correlation
table.tf = as_tibble(unique(summarized_viper_scores$tf))
names(table.tf)="tf"
## We prepare the data of tf acivity
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(table.tf, by = "tf") %>%
  dplyr::select(-std) %>%
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
#join ECM score and tf acivity score for correlation
CorMatrix=cbind(MeanExprECM, summarized_viper_scores_df)

#correlation
cors=as.data.frame(py$correlate_means_to_gene(CorMatrix, corr_gene = "Core_matrisome_score"))
top30.cor.TF <-rownames(cors)[2:31]
CorMatrix=arrange(CorMatrix,Core_matrisome_score)

#perform clustering of the acitiy scores per row
CorMatrix.clustered<-t(CorMatrix[c("Core_matrisome_score",top30.cor.TF)])
d = dist(CorMatrix.clustered, method = "euclidean") # get unique.combined.values from the analysis above
clust <- hclust(d, method = "complete")
CorMatrix.clustered<-CorMatrix.clustered[rev(clust$order),]
CorMatrix.clustered <- CorMatrix.clustered[c(19,1:10),]

EMC.tf.acivity.heatmap <- pheatmap(CorMatrix.clustered,main = "top10 TF activity correlating with Core_Matrisome_score",cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),fontsize_row = 10,show_colnames = F,fontsize = 10)
jpeg(filename = paste0(sample,"_top10 TF activity correlating with Core_Matrisome_score_res10_split.jpeg"), width=800 , height =800,quality = 100,res = 200)
print(EMC.tf.acivity.heatmap)
dev.off()

svg(filename = paste0(sample,"_top10 TF activity correlating with Core_Matrisome_score_res10_split.svg"))
print(EMC.tf.acivity.heatmap)
dev.off()
