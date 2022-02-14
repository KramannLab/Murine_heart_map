#script for runing monocle3 from seurat page vignette https://satijalab.org/signac/articles/monocle.html
library(monocle3)
library(SeuratWrappers)
library(gprofiler2)
library(viridis)
library(scales)
library(enrichR)

#------------convert seurat to cds for monocle3
monocle3.cds <- as.cell_data_set(Samples.combined)

#Calculate size factors using built-in function in monocle3
monocle3.cds <- estimate_size_factors(monocle3.cds)
#Add gene names into CDS
monocle3.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(sc[["RNA"]])
#preprocess
monocle3.cds<-preprocess_cds(monocle3.cds)

#process the cds object
monocle3.cds <- cluster_cells(cds = monocle3.cds, reduction_method = "UMAP",resolution = 0.2)

###instead of clustering from monocle add the seurat clusters here
monocle3.cds@clusters$UMAP$clusters <- Samples.combined$seurat_clusters

monocle3.cds <- learn_graph(monocle3.cds)

plot_cells(monocle3.cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_roots = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,cell_size = 1) + scale_color_manual(values = cluster.col)
ggsave(filename = paste0(sample,"monolce3-clusters.jpeg"), width=10 , height = 10)

# order cells
monocle3.cds <- order_cells(monocle3.cds, reduction_method = "UMAP")

plot_cells(monocle3.cds,cell_size = 1,graph_label_size = 5)
ggsave(filename = paste0(sample,"monolce3-pseudotime_root.jpeg"), width=10 , height = 10)

plot_cells(monocle3.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           trajectory_graph_segment_size = 2,
           graph_label_size = 5,cell_size = 1)
ggsave(filename = paste0(sample,"monolce3-pseudotime_root_dependent.jpeg"), width=10 , height = 10)


#---------select branch from Fib1 to ECM-Fib to test for genes over pseudotime----------------
cds_sub <- choose_graph_segments(monocle3.cds,clear_cds = F)

# plot
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           trajectory_graph_segment_size = 2,
           graph_label_size = 5,cell_size = 1)
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib.jpeg"), width=10 , height = 10)
# extract pseudotime from monocle3 and add to seurat, also plot
pseudotime <- as.data.frame(cds_sub@principal_graph_aux@listData$UMAP$pseudotime[colnames(cds_sub)])
colnames(pseudotime) = "pseudotime_value"
Samples.combined<-AddMetaData(Samples.combined,pseudotime)
FeaturePlot(Samples.combined,features = "pseudotime_value",pt.size = 1)+ scale_colour_gradientn(colours =inferno(20,begin = 0.1,end = 0.9))
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_peudotime.jpeg"), width=10 , height = 10)

FeaturePlot(Samples.combined,features = "pseudotime_value",split.by = "stim",pt.size = 1)
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_peudotime_split.jpeg"), width=30 , height = 10)

#split trajectory analysis into the three conditions and perform analysis for pseudotime associated genes separately
cds_sub.sham = cds_sub[,row.names(subset(pData(cds_sub),stim=="Sham"))]
cds_sub.tac = cds_sub[,row.names(subset(pData(cds_sub),stim=="TAC"))]
cds_sub.tac28 = cds_sub[,row.names(subset(pData(cds_sub),stim=="TAC_28"))]

# test for DEG over pseudotime seperate for each sample (Sham, TAC and TAC28)
## sham
cds_sub.sham_linDEG <- graph_test(cds_sub.sham, neighbor_graph="principal_graph", cores=16)
cds_sub.sham_linDEG <- arrange(cds_sub.sham_linDEG,-morans_test_statistic)
cds_sub.sham_linDEG <- arrange(cds_sub.sham_linDEG,q_value)
top.linDEG.sham <- cds_sub.sham_linDEG$gene_short_name[c(1:500)]
## tac
cds_sub.tac_linDEG <- graph_test(cds_sub.tac, neighbor_graph="principal_graph", cores=16)
cds_sub.tac_linDEG <- arrange(cds_sub.tac_linDEG,-morans_test_statistic)
cds_sub.tac_linDEG <- arrange(cds_sub.tac_linDEG,q_value)
top.linDEG.tac <- cds_sub.tac_linDEG$gene_short_name[c(1:500)]
## tac28
cds_sub.tac28_linDEG <- graph_test(cds_sub.tac28, neighbor_graph="principal_graph", cores=16)
cds_sub.tac28_linDEG <- arrange(cds_sub.tac28_linDEG,-morans_test_statistic)
cds_sub.tac28_linDEG <- arrange(cds_sub.tac28_linDEG,q_value)
top.linDEG.tac28 <- cds_sub.tac28_linDEG$gene_short_name[c(1:500)]

# combined list
overlapp = unique(c(top.linDEG.sham,top.linDEG.tac,top.linDEG.tac28))

# reduce seurat object to cells in the trajectory for ploting
cells.in.path = subset(Samples.combined,cells = rownames(Samples.combined[[]])[(!(is.na(Samples.combined$pseudotime_value)))])

## overall
# generate gene modules based on the hclust for the heatmap
cells.in.path$pseudotime_value_round = round(cells.in.path$pseudotime_value,digits = 1)
cells.in.path$pseudotime_value_round = factor(cells.in.path$pseudotime_value_round,levels = unique(cells.in.path$pseudotime_value_round)[order(unique(cells.in.path$pseudotime_value_round))])
Idents(cells.in.path)="pseudotime_value_round"
avgexp = AverageExpression(cells.in.path, return.seurat = T)
expression.values <- FetchData(avgexp,vars = overlapp,slot = "scale.data")
expression.values <- cluster.pseudotime.genes(expression.values)
genes.to.highlight <- c("Col4a1","Col4a2","Col6a1","Col1a1","Col3a1","Smoc2","Adamts5","Meox1","Mical2","Palld","Cald1","Tgfb1","Thbs4","Postn","Comp","Cilp","Runx1", "Sdc4")
newnames=NULL
for (row.nr in 1:nrow(expression.values)){
  if(row.names(expression.values[row.nr,]) %in% genes.to.highlight) {newnames=c(newnames,bquote(italic(.(row.names(expression.values[row.nr,])))))}else{
    newnames=c(newnames,"")}}
expression.heatmap <- pheatmap(expression.values,main = paste0(sample,"Pseudotime_DEG_overlapp"),labels_row = as.expression(newnames),cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_overall_pseudotime_associated_genes_heatmap.jpeg"), width=1000 , height =1800,quality = 100,res = 200)
print(expression.heatmap)
dev.off()
#save gene as list
write_xlsx(data.frame("gene_name"=rownames(expression.values)),path = paste0(sample,"_list of all pseudotime asociated genes.xlsx"))

#GO terms of gene clusters
expression.values <-FetchData(avgexp,vars = overlapp,slot = "scale.data")
expression.values <- cluster.pseudotime.genes(expression.values,keep.hclust = T)
GO.gene.clusters(expression.values,topGOs.number = 5,only.unique.GOs = F)
ggsave(filename = paste0(sample,"GO-terms associated with genes that correlate with pseudotime_Top5.jpeg"), width=4.5, height = 5)
ggsave(filename = paste0(sample,"GO-terms associated with genes that correlate with pseudotime_Top5.svg"), width=4.5, height = 5)


#load the GO-term extracellular matrix as this was found multiple times to specifically plot the associtated genes
GO_ECM <- read_excel("~/Documents/FP_scRNA/R stuff/references/GO_term_GO_0031012_extracellular_matrix.xlsx")
GO_ECM <- unique(GO_ECM[GO_ECM$`Annotated Term`=="extracellular matrix",]$Symbol)
expression.values <-FetchData(avgexp,vars = overlapp[overlapp %in% GO_ECM],slot = "scale.data")
expression.values <- cluster.pseudotime.genes(expression.values,keep.hclust = T)
rownames(expression.values)<- gsub(pattern = "_","",rownames(expression.values))
expression.values$hcluster = factor(expression.values$hcluster,levels = c(5,2,4,1,3))
expression.values<-arrange(expression.values,hcluster)
expression.values$hcluster=NULL
expression.heatmap <- pheatmap(expression.values,main = "ECM-GO-term genes over pseudotime",labels_row = italic_heatmap_labels(expression.values),cluster_cols = F,cluster_rows = F,angle_col = 45,border_color = 0,color = col.ramp(100),treeheight_row = 20)
jpeg(filename = paste0(sample,"_overall_ECM-GO-term genes over pseudotime.jpeg"), width=1000 , height =2000,quality = 100,res = 200)
print(expression.heatmap)
dev.off()

#test GSEA with enrichr in the pathway database bioplanet

dbs <- c("BioPlanet_2019")
expression.values <-FetchData(avgexp,vars = overlapp,slot = "scale.data")
expression.values <- cluster.pseudotime.genes(expression.values,keep.hclust = T)
enriched.all=NULL
for (k in unique(expression.values$hcluster)) {
  enriched <- enrichr(rownames(expression.values[expression.values$hcluster==k,]), dbs)
  enriched<-enriched$BioPlanet_2019
  enriched$cluster=k
  enriched.all=rbind(enriched.all,enriched)
}
enriched.all.top10 <- enriched.all %>% group_by(cluster) %>% top_n(n = -10,wt=Adjusted.P.value)
enriched.all.top10 <- as.data.frame(enriched.all.top10)
enriched.all.top10$count = parse_number(enriched.all.top10$Overlap)
enriched.all.top10$Term = paste0(rownames(enriched.all.top10),"_",enriched.all.top10$cluster,enriched.all.top10$Term)
enriched.all.top10$Term=factor(enriched.all.top10$Term,levels = rev(unique(enriched.all.top10$Term)))
#selected meaningfull top5 enriched pathways per pseudotime genecluster
enriched.all.selected <- enriched.all.top10[c(1:5,11,12,14,15,16,21,22,23,27,28,33,35,37,38,41,41,44,45,47,48,52),]
#plot
ggplot(enriched.all.selected, aes(x=Term, y="",size=count,color=-log10(Adjusted.P.value))) + 
geom_point() + 
  scale_size(range = c(2,6)) +
  coord_flip() +
  ylab("")+
  xlab("")+
  scale_x_discrete(position = "top")+
  ggtitle(paste0("BioPlanet pathways (EnrichR)"))+
  theme(axis.text.y = ,axis.text.x = element_text(angle=45,hjust=1) ,legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5,size = 10),legend.key.size = unit(0.3, "cm"))
ggsave(filename = paste0(sample,"BioPlanet pathways (EnrichR)_pseudotime_gene_clusters.jpeg"), width=5.8, height = 5)
ggsave(filename = paste0(sample,"BioPlanet pathways (EnrichR)_pseudotime_gene_clusters.svg"), width=5.8, height = 5)



#pseudotime pattern 
p1 <- plot.expression.over.pseudotime(cells.in.path,"Cald1",only.curves = T)
p1
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_Cald1.jpeg"), width=6, height = 5)

p2 <- plot.expression.over.pseudotime(cells.in.path,"Mical2",only.curves = T)
p2
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_Mical2.jpeg"), width=6, height = 5)

p3 <- plot.expression.over.pseudotime(cells.in.path,"Palld",only.curves = T)
p3
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_Palld.jpeg"), width=6, height = 5)

p4 <- plot.expression.over.pseudotime(cells.in.path,"Col4a1",only.curves = T)
p4
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_Col4a1.jpeg"), width=6, height = 5)

p5 <- plot.expression.over.pseudotime(cells.in.path,"Col3a1",only.curves = T)
p5
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_Col3a1.jpeg"), width=6, height = 5)


plotlist = list(p1,p2,p3,p4)
for (plot.nr in c(1:4)) {
  plotlist[[plot.nr]] <- plotlist[[plot.nr]] + theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank())
}
plotlist[[5]] <-p5 + theme(axis.title.y = element_blank())

ggarrange(plotlist = plotlist,common.legend = T,legend = "bottom",ncol = 1,nrow = 5)
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_5GOI.jpeg"), width=5, height = 25)
ggsave(filename = paste0(sample,"Fib1 to ECM-Fib_5GOI.svg"), width=5, height = 25)


# plot the enrichment of celltypes over the selected pseudotime
Idents(cells.in.path)="seurat_clusters"
data.sc <- subset(cells.in.path,ident=c(0,4)) 
R1 <- RidgePlot(data.sc,features = "pseudotime_value",group.by = "stim",cols = stim.col)+ylab("Cell density")+ylab("Cell density")+xlab("Pseudotime")+theme(axis.title.y = element_text(hjust = 0.5),plot.title = element_blank(),axis.title.x = element_text(hjust = 0.5))+NoLegend()
R2 <- RidgePlot(data.sc,features = "pseudotime_value",group.by = "celltypes.short",cols = cluster.col[c(1,5)])+theme(axis.title.y = element_text(hjust = 0.5),plot.title = element_blank(),axis.text.x = element_blank())+NoLegend()
ggarrange(plotlist = list(R2,R1),ncol = 1,nrow = 2)
ggsave(filename = paste0(sample,"Cell density over pseusotime_overall and per Fib1+ECMFib.jpeg"), width=10, height = 5)
ggsave(filename = paste0(sample,"Cell density over pseusotime_overall and per Fib1+ECMFib.svg"), width=10, height = 5)



# plot ECM and collagen subtype scores over pseudotime
# run ECM and collagen scoring in the full fibroblast dataset beforehand
#rescale the scoring between 0 and 100 to allow curve fitting
data.sc_df <- FetchData(data.sc,vars = c("Core_matrisome1","fibrillar.col1","network.forming.col1",))
data.sc_df$fibrillar.col1 <- rescale(data.sc_df$fibrillar.col1,to = c(0,100))
data.sc_df$network.forming.col1 <- rescale(data.sc_df$network.forming.col1,to = c(0,100))
data.sc_df$Core_matrisome1 <- rescale(data.sc_df$Core_matrisome1,to = c(0,100))
data.sc<-AddMetaData(data.sc,metadata = data.sc_df)

#plot the cruves, split by condition
p1 <- plot.expression.over.pseudotime(data.sc,"fibrillar.col1",only.curves = T) + ylab("relativ score") + ggtitle("Fibrillar col. score") +theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p2 <- plot.expression.over.pseudotime(data.sc,"network.forming.col1",only.curves = T) + ylab("relativ score") + ggtitle("Network forming col. score")+theme(axis.text.x = element_blank(),axis.title.x = element_blank())
p3 <- plot.expression.over.pseudotime(data.sc,"Core_matrisome1",only.curves = T) + ylab("relativ score") + ggtitle("ECM score")

ggarrange(plotlist = list(p1,p2,p3),common.legend = T,legend = "right",ncol = 1,nrow = 3)
ggsave(filename = paste0(sample,"ECM_fibri_net_score over pseudotime_split.jpeg"), width=10, height = 10)
ggsave(filename = paste0(sample,"ECM_fibri_net_score over pseudotime_split.svg"), width=10, height = 10)



