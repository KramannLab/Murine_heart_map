#script for functions
#libraries
library(writexl)
library(Seurat)
library(ggplot2)
library(gprofiler2)
library(dplyr)
library(ggrepel)
library(stringr)
require(svglite)
library(readxl)
library(KEGGREST)

#----function to perform clustering based on seurat cca integartion------
recluster <- function(object,res,rep.str=5,loc.con=1){
  DefaultAssay(object) <- "integrated"
  # Run the standard workflow for visualization and clustering
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, verbose = FALSE)
  # t-SNE and Clustering
  object <- RunUMAP(object, reduction = "pca", dims = 1:30,repulsion.strength = rep.str,local.connectivity = loc.con)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object, resolution = res)
  print(DimPlot(object))
  message("plot not saved")
  DefaultAssay(object) <- "RNA"
  return(object)
}

#----function to perform clustering based on RNA-----
recluster_RNA <- function(object,res,rep.str=5,loc.con=1){
  DefaultAssay(object)="RNA"
  # Run the standard workflow for visualization and clustering
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, verbose = FALSE)
  # t-SNE and Clustering
  object <- RunUMAP(object, reduction = "pca", dims = 1:30,repulsion.strength = rep.str,local.connectivity = loc.con)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object, resolution = res)
  print(DimPlot(object))
  message("plot not saved")
  return(object)
}

#----function to get the standard colours of seurat-----
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#----function to perform gene set enrichment analysis based marker genes, calculated by findmarkers-----
get.markergenes.GOs <-function(sc.object,sample,create.folder=T,test="MAST", marker.genes.nr="ALL_up",all.markers=NULL){
  #only for overrepresented markergenes with positiv logFC
  #result can be used with ->plot.markergenes.GOs
  #sc.object-> seurat object with numbered clustering
  #sample-> sample names for saving files
  #create.folder-> in case outputfolder should be created
  #test -> test to use for markergene calculation
  #marker.genes.nr -> number of marker genes to be used for the enrichment test. sometimes using only top50-100 is better for clearer picture
  #all.markers-> in case markergenes were already caculated previously they can be provided to safe time
  
  #check if folder needs to be created
  if (create.folder) {
    dir.create("gprofiler2_out_put")
    dir.create("gprofiler2_out_put/markergenes_GOs")}
  
  #caculate markergenes per seurat_cluster
  Idents(sc.object)<-sc.object$seurat_clusters
  if (is.null(all.markers)) {#in case the markergenes need to calculated first
  all.markers = FindAllMarkers(sc.object, test.use = test,min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA")
  saveRDS(all.markers,file = paste0("gprofiler2_out_put/markergenes_GOs/",sample,"_all.markers_markergenes_per_cluster_result.RDS"))  }
  #sort markergene results by logFC and clusternumber
  all.markers <- arrange(all.markers,avg_logFC)
  all.markers <- arrange(all.markers,cluster)

  for (cluster.nr in levels(sc.object$seurat_clusters)) {
    #GO of marker genes that are significant up
    genes_up = all.markers$gene[all.markers$avg_logFC>0 & all.markers$cluster==cluster.nr & all.markers$p_val_adj < 0.01]
    all.markers <- arrange(all.markers,p_val)

    #optional filter for only topXX markergenes to apply GO enrichment 
    if (marker.genes.nr != "ALL_up") {genes_up <- genes_up[1:marker.genes.nr]}
    
    #exclude Ribosomal genes form GO analysis
    genes_up = setdiff(genes_up,genes_up %>%  str_subset(pattern = "^Rp"))
    
    if (length(genes_up)<5) {
      message(paste0("cluster: ",cluster.nr," has less then 5 markergenes, therefore no GO testing"))
      go_up=NULL
    } else {
      #run gProfiler for genes that are UP in treatment_condition per cluster
      go_up <- gost(query = genes_up, 
                    organism = "mmusculus", ordered_query = T, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE) 
      }

    #transfer to temporary dataframe and add column to save cluster number and UP or Down regulation
    if (is.null(go_up)) {
      message(paste0("for cluster: ", cluster.nr," no GOs for markergenes found"))
      go_up.new.results=NULL
    } else {
      go_up.new.results = go_up$result
      go_up.new.results$cluster = as.numeric(cluster.nr)
      }
    
    #add results to overall dataframe
    if (cluster.nr == 0) {
      message(paste0("cluster: ",cluster.nr," now saved"))
      go.results <- go_up.new.results
      
    } else {
      message(paste0("cluster: ",cluster.nr," else for rbind"))
      go.results <- rbind(go.results,go_up.new.results)
    }
    
  }#end of for
  #save overall result
  saveRDS(go.results,file = paste0("gprofiler2_out_put/markergenes_GOs/",sample,"GOs per cluster based on_",marker.genes.nr,"_markergenes.rds"))
  write_xlsx(go.results,path = paste0("gprofiler2_out_put/markergenes_GOs/",sample,"GOs per cluster based on_",marker.genes.nr,"__markergenes.xlsx"), col_names = TRUE)
  return(go.results)
}#end of get.markergenes.GOs

#----two functions that help with visualization to find nice width and height by calculation------
GO.plot.width.heigt <- function(){
  perfect.width <<- function(gene.list = data.frame(),coord_fliped = logical()){
    if (coord_fliped == TRUE) {
      longest.term = unique(nchar(gene.list$term_name[nchar(gene.list$term_name)==max(nchar(gene.list$term_name))]))
      number.of.cl = length(levels(factor(gene.list$cluster)))
      width = round(number.of.cl*0.22 + longest.term*0.06,digits = 1)+1
      message(paste0("width is: ",width))
      return(width)
    } else {
      number.of.gos = length((unique(gene.list$term_name)))
      width = (number.of.gos*1.3/10)+2
      message(paste0("width is: ",width))
      return(width)
    }
  }
  perfect.height <<- function(gene.list = data.frame(),coord_fliped = logical()){
    if (coord_fliped == TRUE) {
      number.of.gos = length((unique(gene.list$term_name)))
      if (number.of.gos > 60) {
        message(paste0("high number of GOs: ",number.of.gos))
        height = (number.of.gos*1.3/10)
      } else {
        height = (number.of.gos*1.8/10)
      }
      message(paste0("height is: ",height))
      return(height)
    } else {
      longest.term = unique(nchar(gene.list$term_name[nchar(gene.list$term_name)==max(nchar(gene.list$term_name))]))
      number.of.cl = length(levels(factor(gene.list$cluster)))
      height = round(number.of.cl/4 + longest.term*0.05,digits = 1)
      message(paste0("height is: ",height))
      return(height)
    }
  }
}

#----function to make a dotplot of gene sets found enriched in the marker genes------
plot.markergenes.GOs <- function(sample,go.results,topGOs.number = 10,marker.genes.nr="ALL_up",custom.h=0,custom.w=0,provide.label=NULL,uniqueGO.filter=F){
  #only for overrepresented markergenes with positiv logFC
  #sc.object-> seurat object with numbered clustering
  #sample-> sample names for saving files
  #go.results-> calculated by get.markergenes.GOs function
  #topGOs.number-> how many top GOs should be included the final plot (ploting is more optimized for 10)
  #custom.h&w-> in case the perfect.width and height functions fail, this can be used to customize the graph dimensions for saving
  #provide.label-> in case x-axis label should be based on annotation provide a list of clusters+annotation e.g: "0"="Cap"
  #uniqueGO.filter -> in case the provided GOterms were filtered for unique ones per cluster
  GO.plot.width.heigt()
  if(uniqueGO.filter){
    per.cluster.unique.GOs<-NULL
    for (cluster in unique(go.results$cluster)) {
      all.GOs.found <- unique(go.results$term_id[go.results$cluster!=cluster])
      per.cluster.unique.GOs<-rbind(per.cluster.unique.GOs,go.results[go.results$cluster==cluster & !(go.results$term_id %in% all.GOs.found),])
    }#end of for
    go.results <- per.cluster.unique.GOs
    u="unique"
  }else{u=""}
  
  #----------------generate a dataframe of only the topXX GOs per cluster---------------
  go.results <- go.results[go.results$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results <- go.results[go.results$term_size<1000,] #filter out large and often very general GOs
  go.results <- arrange(go.results,p_value)
  go.results <- arrange(go.results,cluster)

  #------extract and plot the top10 by p_value per cluster
  go.results.top10 <- go.results %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  go.results.top10<-arrange(go.results.top10,cluster)
  go.results.top10$cluster=factor(go.results.top10$cluster)
  order_desc_up <- unique(go.results.top10$term_name)
  go.results.top10.all <- go.results[go.results$term_name %in% order_desc_up,]
  x.angle=0
  x.hjust=0.5
  #in case labels are provided
  if (!(is.null(provide.label))) {
    go.results.top10.all$cluster=as.character(go.results.top10.all$cluster)
    for (cluster.nr in names(provide.label)) {go.results.top10.all[["cluster"]]<-ifelse(go.results.top10.all$cluster==cluster.nr,provide.label[[cluster.nr]],go.results.top10.all$cluster)}
    provide.label <- provide.label[unlist(provide.label)%in%go.results.top10.all$cluster]
    go.results.top10.all$cluster = factor(go.results.top10.all$cluster,levels = unlist(provide.label))
    x.angle=45
    x.hjust=1
  }#end of if
  
  
  #-plot up vertikal
  print(ggplot(go.results.top10.all, aes(x=factor(term_name, level = order_desc_up), y=factor(cluster),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster"))+
    theme(axis.text.x = element_text(angle=x.angle,hjust=x.hjust) ,legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5,size = 10),legend.key.size = unit(0.3, "cm"))
  )#end of print
  
  if (custom.h==0 & custom.w==0) {
    ggsave(filename = paste0("gprofiler2_out_put/markergenes_GOs/",sample,"Top",topGOs.number," ",u," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.jpeg"), width=perfect.width(go.results.top10.all,T) , height = perfect.height(go.results.top10.all,T))
    ggsave(filename = paste0("gprofiler2_out_put/markergenes_GOs/",sample,"Top",topGOs.number," ",u," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.svg"), width=perfect.width(go.results.top10.all,T) , height = perfect.height(go.results.top10.all,T))
  }else{
    ggsave(filename = paste0("gprofiler2_out_put/markergenes_GOs/",sample,"Top",topGOs.number," ",u," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.jpeg"), width=custom.w , height = custom.h)
    ggsave(filename = paste0("gprofiler2_out_put/markergenes_GOs/",sample,"Top",topGOs.number," ",u," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.svg"), width=custom.w , height = custom.h)
  }#end of if else
}#end of plot.markergenes.GOs

#----calculated DEG between two condition in a seurat object------
get.DEG.between.conditions.per.cluster <-function(sc.object,sample,test="MAST",control_sample,treatment_sample){
  #test DE genes for all clusters between two conditions (in stim) and save summarized
  #sc.object-> seurat object with numbered clustering
  #sample-> sample names for saving files
  #test -> test to use for markergene calculation
  #control&treatment_sample-> names of the two conditions in stim
  
  dir.create("DE_test_conditions_per_cluster")
  #create a new column for joining cluster and condition in order to compare the conditions later within one cluster
  sc.object$cluster.stim <- paste(sc.object$seurat_clusters, sc.object$stim, sep = "_")
  Idents(sc.object) <- "cluster.stim"
  for (cluster.nr in levels(sc.object$seurat_clusters)) {
    #make sure there are enough cells in the cluster to do DE
    cellnumber.control = length(sc.object$cluster.stim[sc.object$cluster.stim==paste0(cluster.nr,"_",control_sample)])
    cellnumber.treatment = length(sc.object$cluster.stim[sc.object$cluster.stim==paste0(cluster.nr,"_",treatment_sample)])
    
    if (cellnumber.control > 10 & cellnumber.treatment > 10) {
      stim.response.new <- FindMarkers(sc.object,assay = "RNA",test.use = test,min.pct=0.25,logfc.threshold = 0.3 ,ident.1 = paste0(cluster.nr,"_",treatment_sample), ident.2 = paste0(cluster.nr,"_",control_sample), verbose = FALSE)
      stim.response.new$gene_name=row.names(stim.response.new)
      stim.response.new <- arrange(stim.response.new,-avg_logFC)
      stim.response.new$cluster = cluster.nr
      stim.response.new <- stim.response.new[stim.response.new$p_val_adj < 0.01,] #only keep significant genes
      
      #collect all results in a dataframe of stim.response
      if (cluster.nr == min(levels(sc.object$seurat_clusters))) {
        message(paste0("cluster: ",cluster.nr," DE result now saved"))
        stim.response <- stim.response.new
      } else {
        message(paste0("cluster: ",cluster.nr," DE result for rbind"))
        stim.response <- rbind(stim.response,stim.response.new)
      }
      
      #save the results of the DE testing itself
      write_xlsx(stim.response.new,path = paste0("DE_test_conditions_per_cluster/",sample,"_DE_analysis_with_",test,"_for_cluster_",cluster.nr,".xlsx"), col_names = TRUE)
    }else{
      message(paste("Custer_",cluster.nr,"has only ",cellnumber.control," cell for control and ",cellnumber.treatment," cells for treatment"))  
    }
  }
  #save overall DE result
  saveRDS(stim.response,file = paste0("DE_test_conditions_per_cluster/",sample,"_test:",test,"_FindMarkers_",control_sample,"_vs_",treatment_sample,".rds"))
  write_xlsx(stim.response,path = paste0("DE_test_conditions_per_cluster/",sample,"_DE_analysis_with_",test,"_all_DE_genes.xlsx"), col_names = TRUE)
  stim.response$cluster=factor(stim.response$cluster,levels = levels(sc.object$seurat_clusters))
  return(stim.response)
}#end of get.DEG.between.conditions.per.cluster

#----check results from get.DEG.between.conditions.per.cluster for genes that were only DEG in one cluster (or might be rather batch related)------
check.for.unique.DE <- function(stim.response,overall.DE=NULL){
  #overall.DE -> can be provided if it should be tested which genes are also overall DE and which only in a subcluster
  if (is.null(overall.DE)) {
  #check all genes if they are unique DE for a cluster
  for (cluster.nr in levels(factor(stim.response$cluster))) {
    DE.genes.up = unique(stim.response$gene_name[stim.response$avg_logFC>0& stim.response$p_val_adj < 0.01&stim.response$cluster!=cluster.nr]) #genes of all other cluster
    DE.genes.up.cl = stim.response$gene_name[stim.response$avg_logFC>0& stim.response$p_val_adj < 0.01&stim.response$cluster==cluster.nr] #genes of the cluster to look at
    
    DE.genes.down = unique(stim.response$gene_name[stim.response$avg_logFC<0& stim.response$p_val_adj < 0.01&stim.response$cluster!=cluster.nr]) #see above
    DE.genes.down.cl = stim.response$gene_name[stim.response$avg_logFC<0& stim.response$p_val_adj < 0.01&stim.response$cluster==cluster.nr] #see above
    
    DE.genes.all = unique(stim.response$gene_name[stim.response$cluster!=cluster.nr& stim.response$p_val_adj < 0.01]) #see above
    DE.genes.all.cl = stim.response$gene_name[stim.response$cluster==cluster.nr& stim.response$p_val_adj < 0.01] #see above
    
    if (cluster.nr == 0) { #create the result list starting with cluster 0
      genes.de.up.unique <- stim.response[stim.response$cluster==cluster.nr & stim.response$gene_name %in% setdiff(DE.genes.up.cl,DE.genes.up),] #compare DE genes of a cluster to the DE genes of all other clusters to identify unique DE genes
      genes.de.down.unique <- stim.response[stim.response$cluster==cluster.nr & stim.response$gene_name %in% setdiff(DE.genes.down.cl,DE.genes.down),] #see above
      genes.de.all.unique <- stim.response[stim.response$cluster==cluster.nr & stim.response$gene_name %in% setdiff(DE.genes.all.cl,DE.genes.all),] #see above
      
    } else { #add genes+cluster to the list
      genes.de.up.unique.new <- stim.response[stim.response$cluster==cluster.nr & stim.response$gene_name %in% setdiff(DE.genes.up.cl,DE.genes.up),] #see above
      genes.de.down.unique.new <- stim.response[stim.response$cluster==cluster.nr & stim.response$gene_name %in% setdiff(DE.genes.down.cl,DE.genes.down),] #see above
      genes.de.all.unique.new <- stim.response[stim.response$cluster==cluster.nr & stim.response$gene_name %in% setdiff(DE.genes.all.cl,DE.genes.all),] #see above
      
      genes.de.up.unique <- rbind(genes.de.up.unique,genes.de.up.unique.new)
      genes.de.down.unique <- rbind(genes.de.down.unique,genes.de.down.unique.new)
      genes.de.all.unique <- rbind(genes.de.all.unique,genes.de.all.unique.new)
    }
  }
  #adding to stim.response if the gene is unique for the cluster or not
  for (gene in unique(stim.response$gene_name)) {
    if (gene %in% c(genes.de.up.unique$gene_name,genes.de.down.unique$gene_name)) { #check if unique at all
      if (gene %in% genes.de.up.unique$gene_name) { #check if it is unique in UP
        stim.response$unique[stim.response$gene_name==gene & stim.response$avg_logFC>0] = paste0("unique")
      }else { #if its not in unique UP it must be in unique DOWN
        stim.response$unique[stim.response$gene_name==gene & stim.response$avg_logFC<0] = paste0("unique")
      }
    } else {
      stim.response$unique[stim.response$gene_name==gene] ="shared"
    }
  }
  stim.response$unique[is.na(stim.response$unique)] ="shared" #all NAs in the unique check must be shared in UP but unique in down and wise versa
  }else{#end of is.null(overall.DE)
    stim.response[['unique']] <- ifelse(stim.response$gene_name %in% rownames(overall.DE),'shared','unique')
  }#end of else
  return(stim.response)
}#end of check.for.unique.DE

#----function to perform gene set enrichment analysis based DEG (take unique/shared per cluster into account)-----
get.GOs.based.unique_shared.DEG <- function(stim.response,test="MAST",topGOs.number=10){
  #stim.response -> result from per cluster DE testing + check.for.unique.DE
  ####-----------GOs based only on the UNIQUE  genes-------------###########################################################################
  for (cluster.nr in levels(factor(stim.response$cluster))){
    if (cluster.nr == 0) {go.results.unique <- NULL}
    
    #testing for GOs that are unique per clster -> here are already just significant genes included
    genes_up.unique = stim.response$gene_name[stim.response$cluster==cluster.nr & stim.response$unique=="unique" & stim.response$avg_logFC>0]
    genes_down.unique = stim.response$gene_name[stim.response$cluster==cluster.nr & stim.response$unique=="unique"& stim.response$avg_logFC<0]
    stim.response <- arrange(stim.response,p_val)
    genes_all.unique = stim.response$gene_name[stim.response$cluster==cluster.nr & stim.response$unique=="unique"]
    
    #exclude Ribosomal genes form GO analysis
    genes_up.unique = setdiff(genes_up.unique,genes_up.unique %>%  str_subset(pattern = "^Rp"))
    genes_down.unique = setdiff(genes_down.unique,genes_down.unique %>%  str_subset(pattern = "^Rp"))
    genes_all.unique = setdiff(genes_all.unique,genes_all.unique %>%  str_subset(pattern = "^Rp"))
    
    #check if there are at least 5 unique up genes for GO
    if (length(genes_up.unique)<5) {
      message(paste0("cluster: ",cluster.nr," has less then 5 up unique DE genes, therefore no GO testing"))
      go_up.unique = NULL
    } else {
      #run gProfiler for genes that are UP in treatment_sample per cluster
      go_up.unique <- gost(query = genes_up.unique, 
                           organism = "mmusculus", ordered_query = T, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = FALSE, 
                           user_threshold = 0.05, correction_method = "g_SCS", 
                           domain_scope = "annotated", custom_bg = NULL, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE)}
    
    #check if there are at least 5 unique down genes for GO
    if (length(genes_down.unique)<5) {
      message(paste0("cluster: ",cluster.nr," has less then 5 down unique DE genes, therefore no GO testing"))
      go_down.unique = NULL
    } else {
      #run gProfiler for genes that are UP in treatment_condition per cluster
      go_down.unique <- gost(query = rev(genes_down.unique), 
                             organism = "mmusculus", ordered_query = T, 
                             multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                             measure_underrepresentation = FALSE, evcodes = FALSE, 
                             user_threshold = 0.05, correction_method = "g_SCS", 
                             domain_scope = "annotated", custom_bg = NULL, 
                             numeric_ns = "", sources = NULL, as_short_link = FALSE)}
    
    #check if there are at least 5 unique genes for GO
    if (length(genes_all.unique)<5) {
      message(paste0("cluster: ",cluster.nr," has less then 5 unique DE genes, therefore no GO testing"))
      go_all.unique=NULL
    } else {
      #run gProfiler for genes that are UP in treatment_condition per cluster
      go_all.unique <- gost(query = genes_all.unique, 
                            organism = "mmusculus", ordered_query = T, 
                            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                            measure_underrepresentation = FALSE, evcodes = FALSE, 
                            user_threshold = 0.05, correction_method = "g_SCS", 
                            domain_scope = "annotated", custom_bg = NULL, 
                            numeric_ns = "", sources = NULL, as_short_link = FALSE)}
    
    
    #transfer to temporary dataframe and add colum to save cluster number and UP or Down regulation
    if (is.null(go_up.unique)) {
      message(paste0("for cluster: ", cluster.nr," no unique up GOs found"))
      go_up.new.results.unique=NULL
    } else {
      go_up.new.results.unique = go_up.unique$result
      go_up.new.results.unique$cluster = as.numeric(cluster.nr)
      go_up.new.results.unique$up_or_down = "up"}
    if (is.null(go_down.unique)) {
      message(paste0("for cluster: ", cluster.nr," no unique down GOs found"))
      go_down.new.results.unique=NULL
    } else {
      go_down.new.results.unique = go_down.unique$result
      go_down.new.results.unique$cluster = as.numeric(cluster.nr)
      go_down.new.results.unique$up_or_down = "down"}
    if (is.null(go_all.unique)) {
      message(paste0("for cluster: ", cluster.nr," no unique GOs found"))
      go_all.new.results.unique=NULL
    } else {
      go_all.new.results.unique = go_all.unique$result
      go_all.new.results.unique$cluster = as.numeric(cluster.nr)
      go_all.new.results.unique$up_or_down = "all"}
    
    #add results to overall dataframe
    if (cluster.nr == 0) {
      message(paste0("cluster: ",cluster.nr," now saved"))
      go.results.unique <- rbind(go_up.new.results.unique,go_down.new.results.unique,go_all.new.results.unique)
      
    } else {
      message(paste0("cluster: ",cluster.nr," else for rbind"))
      go.results.unique <- rbind(go.results.unique,go_up.new.results.unique,go_down.new.results.unique,go_all.new.results.unique)
    }
  }#end of for loop
  
  #save GO based on unique DE list
  write_xlsx(go.results.unique,path = paste0(sample,"top",topGOs.number,"_GOs_unique_for_a_cluster.xlsx"), col_names = TRUE)
  saveRDS(go.results.unique,file = paste0(sample,"_test:",test,"_GOs_unique_for_a_cluster.rds"))
  
  ####-----------GOs based only on the SHARED genes and adding to common list with UNIQUE-------###########################################
  #generate a dataframe of DE genes per cluster keeping only the unique DE genes per cluster
  gene.list.up.shared= stim.response$gene_name[stim.response$unique=="shared" & stim.response$avg_logFC>0]
  gene.list.down.shared= stim.response$gene_name[stim.response$unique=="shared" & stim.response$avg_logFC<0]
  stim.response <- arrange(stim.response,p_val)
  gene.list.all.shared= stim.response$gene_name[stim.response$unique=="shared"]
  
  #exclude Ribosomal genes form GO analysis
  gene.list.up.shared = setdiff(gene.list.up.shared,gene.list.up.shared %>%  str_subset(pattern = "^Rp"))
  gene.list.down.shared = setdiff(gene.list.down.shared,gene.list.down.shared %>%  str_subset(pattern = "^Rp"))
  gene.list.all.shared = setdiff(gene.list.all.shared,gene.list.all.shared %>%  str_subset(pattern = "^Rp"))
  
  #test for GOs
  go_up.shared <- gost(query = gene.list.up.shared, 
                       organism = "mmusculus", ordered_query = F, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = FALSE, 
                       user_threshold = 0.05, correction_method = "g_SCS", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE)
  #test for GOs
  go_down.shared <- gost(query = rev(gene.list.down.shared), 
                         organism = "mmusculus", ordered_query = F, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "g_SCS", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE)
  #test for GOs
  go_all.shared <- gost(query = gene.list.all.shared, 
                        organism = "mmusculus", ordered_query = F, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "g_SCS", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
  #transfer to temporary dataframe and add colum to save cluster number and UP or Down regulation
  go_up.shared.results= go_up.shared$result
  go_up.shared.results$cluster = "shared"
  go_up.shared.results$up_or_down = "up"
  
  go_down.shared.results= go_down.shared$result
  go_down.shared.results$cluster = "shared"
  go_down.shared.results$up_or_down = "down"
  
  go_all.shared.results= go_all.shared$result
  go_all.shared.results$cluster = "shared"
  go_all.shared.results$up_or_down = "all"
  
  #join to one list
  go_up.shared.and.unique.results <- rbind(go.results.unique,go_up.shared.results,go_down.shared.results,go_all.shared.results)
  
  return(go_up.shared.and.unique.results)
}#end of function get.GOs.based.unique_shared.DEG

#----function to make a dotplot of gene sets found enriched in the DEG, split in unique and shared DEG per cluster------
plot.GOs.based.unique_shared.DEG <- function(go.results,treatment_sample,sample,topGOs.number=10,provide.label=NULL,custom.h=0,custom.w=0){
  #goresults -> output from get.GOs.based.unique_shared.DEG
  #provide.label-> in case x-axis label should be based on annotation provide a list of clusters+annotation e.g: "0"="Cap"
  
  GO.plot.width.heigt()
  #########################--GOs based on SHARED + UNIQUE DE genes per cluster-------########################################
  #####--------plotting GOs based on only UP unique genes-----------------
  go.results.unique.filtered <- go.results[go.results$up_or_down =="up",]
  go.results.unique.filtered <- go.results.unique.filtered[go.results.unique.filtered$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results.unique.filtered <- go.results.unique.filtered[go.results.unique.filtered$term_size<1000,] #filter out large and often very general GOs
  
  #soring was done when joining unique and shared gos
  #------extract and plot the topXX by p_value per cluster
  go.results.unique.filtered.top10 <- go.results.unique.filtered %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  go.results.unique.filtered.top10 <- arrange(go.results.unique.filtered.top10,cluster)
  order_desc_up <- unique(go.results.unique.filtered.top10$term_name)
  go.results.unique.filtered.top10.all <- go.results.unique.filtered[go.results.unique.filtered$term_name %in% order_desc_up,]
  go.results.unique.filtered.top10.all$cluster <- factor(go.results.unique.filtered.top10.all$cluster,levels = c(levels(factor(as.numeric(c(na.omit(parse_number(levels(factor(go.results.unique.filtered.top10.all$cluster)))))))),"shared"))
  go.results.unique.filtered.top10.all<-arrange(go.results.unique.filtered.top10.all,cluster)
  x.angle=0
  x.hjust=0.5
  #in case labels are provided
  if (!(is.null(provide.label))) {
    provide.label[["shared"]]="shared"
    go.results.unique.filtered.top10.all$cluster=as.character(go.results.unique.filtered.top10.all$cluster)
    for (cluster.nr in names(provide.label)) {go.results.unique.filtered.top10.all[["cluster"]]<-ifelse(go.results.unique.filtered.top10.all$cluster==cluster.nr,provide.label[[cluster.nr]],go.results.unique.filtered.top10.all$cluster)}
    provide.label <- provide.label[unlist(provide.label)%in%go.results.unique.filtered.top10.all$cluster]
    go.results.unique.filtered.top10.all$cluster = factor(go.results.unique.filtered.top10.all$cluster,levels = unlist(provide.label))
    x.angle=45
    x.hjust=1
    }#end of if
  
  ggplot(go.results.unique.filtered.top10.all, aes(x=factor(term_name, level = rev(order_desc_up)), y=cluster,size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs upregulated in ",treatment_sample," per cluster based unique + shared upregulated DE genes"))+
    theme(axis.text.x = element_text(angle = x.angle,hjust = x.hjust),legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0,size = 10),legend.key.size = unit(0.2, "cm"),legend.box.just = "left") + scale_x_discrete(position = "top") #+theme(legend.position = c(5.5, 0.28))
  if (custom.h==0 & custom.w==0) {
    ggsave(filename = paste0(sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.jpeg"), width=perfect.width(go.results.unique.filtered.top10.all,T)+1 , height = perfect.height(go.results.unique.filtered.top10.all,T))
    ggsave(filename = paste0(sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.svg"), width=perfect.width(go.results.unique.filtered.top10.all,T)+1 , height = perfect.height(go.results.unique.filtered.top10.all,T))
    }else{
    ggsave(filename = paste0(sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.jpeg"), width=custom.w , height =custom.h)
    ggsave(filename = paste0(sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.svg"), width=custom.w , height = custom.h)
  }#end of else
  message("GO plot saved")
}#end of function plot.GOs.based.unique_shared.DEG

#----function to compare cluster specific DEG and overall sample DEG-----
find.batch <- function(DE.results,condition.markers,control_sample,treatment_sample){
  #this function takes DE results of comparision within a cluster and checks if the genes found there are also found be just comparing conditions global
  #these genes have a potential to be rather unspecific DE genes (sample or batch effects)
  condition.markers$gene = rownames(condition.markers)
  DE.results$batch.eff = "specific"
  DE.results$batch.eff[DE.results$gene_name %in% condition.markers$gene]="batch"
  return(DE.results)
}#end of find batch

#----function to score for the expression of several ECM related gene sets-----
scoreECM <- function(sc.object){
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
  DefaultAssay(sc.object)="RNA"
  ctrl_genes = 35
  for (gset in names(matrisome_mm_genesetlist)){
    features = matrisome_mm_genesetlist[gset]
    sc.object = AddModuleScore(object = sc.object, features = features, name = gset, ctrl = ctrl_genes)
  }
  return(sc.object)
}#end of scoreECM

#----function to score for the expression of several functional collagen subgroups-----
scoreCollagens <- function(sc.object){
  #collagens: https://onlinelibrary.wiley.com/doi/10.1111/liv.14390
  fibrillar.col = c("Col1a1","Col1a2", "Col2a1","Col3a1","Col5a1","Col5a2","Col5a3","Col11a1","Col24a1","Col27a1")
  network.forming.col = c("Col4a1","Col4a2","Col4a3","Col4a4","Col4a5","Col4a6","Col6a1","Col6a2","Col6a3","Col6a4","Col6a5","Col6a6","Col7a1","Col8a1","Col8a2","Col10a1")
  FACIT.col = c( "Col9a1","Col9a2" , "Col9a3" ,"Col12a1","Col14a1","Col19a1","Col20a1","Col21a1","Col22a1") #fibril‐associated collagens with interrupted triple helices (FACIT) 
  transmembrane.col = c("Col25a1","Col13a1","Col17a1","Col23a1")
  multiplexin.col = c("Col15a1","Col16a1","Col18a1")#Two subfamilies of collagens, namely the plasma membrane-associated collagens with interrupted triple-helices (MACITs, including ColXIII, ColXXIII and ColXXV) and the basement membrane-associated collagens with multiple triple-helix domains with interruptions (multiplexins, including ColXV and ColXVIII)
  
  col.functional.groups = list(fibrillar.col,network.forming.col,FACIT.col,transmembrane.col,multiplexin.col)
  names(col.functional.groups) = c("fibrillar.col","network.forming.col","FACIT.col","transmembrane.col","multiplexin.col")
  #  scoring for the seurat object
  DefaultAssay(sc.object)="RNA"
  ctrl_genes = 35 #important
  
  for (gset in names(col.functional.groups)){
    features = col.functional.groups[gset]
    sc.object = AddModuleScore(object = sc.object, features = features, name = gset, ctrl = ctrl_genes)
  }
  return(sc.object)
}#end of score Collagens

#----function to score for the expression of stress related genes-----
scoreStress <- function(sc.object){
  #-----------------STRESS score
  #geneset from: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02048-6
  #Fosb, Fos, Jun, Junb, Jund, Atf3, Egr1, Hspa1a, Hspa1b, Hsp90ab1, Hspa8, Hspb1, Ier3, Ier2, Btg1, Btg2, Dusp1
  
  #geneset from: https://www.nature.com/articles/s41598-019-45842-4
  #Brd2, Dnaja1, Dnajb1, Egr1, Hsp90aa1, Hsp90ab1, Hspe1, Hspj1, Ier3, Nr4a1, Nr4a2, Scl38a2
  
  #geneset from:https://www.nature.com/articles/nmeth.4437/figures/1
  #Junb, Fosb, Fos, Jun, Zfp36, Egr1, Hspa1a, Hspa1b, Hspa8, Hspb1, Cebpd, Cebpb, Atf3, Socs3, Jund, Hspe1, Hsp90ab1
  
  stress.markers = c("Fosb", "Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1", "Hspa1a", "Hspa1b", "Hsp90ab1", "Hspa8", "Hspb1", "Ier3", "Ier2", "Btg1", "Btg2", "Dusp1","Brd2", "Dnaja1", "Dnajb1", "Egr1", "Hsp90aa1", "Hsp90ab1", "Hspe1",  "Ier3", "Nr4a1","Junb", "Fosb", "Fos", "Jun", "Zfp36", "Egr1", "Hspa1a", "Hspa1b", "Hspa8", "Hspb1", "Cebpd", "Cebpb", "Atf3", "Socs3", "Jund", "Hspe1", "Hsp90ab1")
  stress.markers = unique(stress.markers)
  features = list(stress.markers)
  sc.object = AddModuleScore(object = sc.object, features = features, name = "STRESSscore", ctrl = 35,search = T)
  return(sc.object)
}#end of scoreStress

#----function to score for the expression of EndMA related genes-----
scoreEndMA <- function(sc.object){
  #load EndMA genes from #https://www.nature.com/articles/s41467-021-20905-1
  EndMA_gene_set <- c("Tgbr2","Mmp14","Fn1","Postn","Eln","Vim","Col1a2","Col1a1","Mgp","Col3a1","Bgn","Dcn")
    sc.object = AddModuleScore(object = sc.object, features = list(EndMt_gene_set), name = "EndMAscore", ctrl = 35)
  return(sc.object)
}#end of scoreEndMA

#----function to get nice clustering and ordering of the genes associated with pseudotime-----
cluster.pseudotime.genes <- function(expression.values,k.treecut=5,keep.hclust=F){
  #expression.values: expression values from FetchData
  #k.treecut: cutoff value for the clustertree, based on this there will be more or less groups of genes
  expression.values <- as.data.frame(t(expression.values))
  d = dist(expression.values, method = "euclidean")
  clust <- hclust(d, method = "complete")
  plot(clust, labels = FALSE)
  clust <- cutree(clust, k = k.treecut) %>%  data.frame()
  names(clust)="hcluster"
  expression.values <- cbind(expression.values,clust)
  expression.values <- arrange(expression.values,hcluster)
  
  order.df=NULL
  for (k in levels(factor(expression.values$hcluster))) {
    k_sub = expression.values[expression.values$hcluster==k,]
    k_sub$hcluster=NULL
    k_sub.means = colMeans(k_sub)
    df = data.frame("ps"=names(k_sub.means[k_sub.means == max(k_sub.means)]),"k"=k)
    order.df=rbind(order.df,df)
  }
  order.df <- arrange(order.df,ps)
  expression.values$hcluster <- factor(expression.values$hcluster,levels =order.df$k )
  expression.values <- arrange(expression.values,hcluster)
  if (keep.hclust) {
    return(expression.values) 
  }else{
  expression.values$hcluster = sapply(expression.values$hcluster,function(k) strrep("_",k) )
  rownames(expression.values) = paste0(expression.values$hcluster,rownames(expression.values))
  expression.values$hcluster=NULL
  return(expression.values)}#end of else
}

#----function to perform gene set enrichment analysis based clustered genes and generate dotplot-----
GO.gene.clusters <-function(expression.values,topGOs.number=10,only.unique.GOs=F){
  #expression.values: from fetch data and run through the cluster.pseudotime.genes function with keep.hclust=T
  gene.cluster.gos = NULL
  global.lvls=NULL
  if(!(only.unique.GOs)){
  for (k in levels(expression.values$hcluster)) {
    genelist = rownames(expression.values[expression.values$hcluster==k,])
    go <- gost(query = genelist, 
               organism = "mmusculus", ordered_query = F, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "g_SCS", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE)
    go.results <- go$result
    if (!(is.null(go.results))) {
      go.results <- go.results[go.results$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
      go.results <- go.results[go.results$term_size<1000,] #filter out large and often very general GOs
      go.results <- arrange(go.results,p_value)
      if ((topGOs.number+1)<length(go.results$query)) { go.results <- go.results[c(1:topGOs.number),]}else{go.results <- go.results[c(1:length(go.results$query)),] }
      go.results$term_name <- paste0(k,"_",go.results$term_name)
      go.results$term_name <- factor(go.results$term_name,levels = go.results$term_name)
      global.lvls<-c(global.lvls,levels(go.results$term_name))
      gene.cluster.gos = rbind(gene.cluster.gos,go.results)
      gene.cluster.gos$term_name <- factor(gene.cluster.gos$term_name,levels = rev(global.lvls)) 
    }#end of if
  }#end of for loop
    u=""
    }#end of if only unique GOs=FALSE
  if(only.unique.GOs){
    for (k in levels(expression.values$hcluster)) {
      genelist = rownames(expression.values[expression.values$hcluster==k,])
      go <- gost(query = genelist, 
                 organism = "mmusculus", ordered_query = F, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = FALSE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL, as_short_link = FALSE)
      go.results <- go$result
      if (!(is.null(go.results))) {
        go.results <- go.results[go.results$source %in% c("KEGG"),] #only take most relevant database
        #go.results <- go.results[go.results$term_size<1000,] #filter out large and often very general GOs
        go.results <- arrange(go.results,p_value)
        go.results$cluster =  k
        go.results$term_name <- paste0(k,"_",go.results$term_name)
        go.results$term_name <- factor(go.results$term_name,levels = go.results$term_name)
        global.lvls<-c(global.lvls,levels(go.results$term_name))
        gene.cluster.gos = rbind(gene.cluster.gos,go.results)
        gene.cluster.gos$term_name <- factor(gene.cluster.gos$term_name,levels = rev(global.lvls))
      }#end of if
    }#end of for loop
      per.cluster.unique.GOs<-NULL
      for (cluster in unique(gene.cluster.gos$cluster)) {
        all.GOs.found <- unique(gene.cluster.gos$term_id[gene.cluster.gos$cluster!=cluster])
        per.cluster.unique.GOs<-rbind(per.cluster.unique.GOs,gene.cluster.gos[gene.cluster.gos$cluster==cluster & !(gene.cluster.gos$term_id %in% all.GOs.found),])
      }#end of for
      gene.cluster.gos <- per.cluster.unique.GOs
      u="unique"
      gene.cluster.gos <- gene.cluster.gos %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  }#end of if only unique GOs =TRUE
  ggplot(gene.cluster.gos, aes(x=term_name, y="",size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("")+
    xlab("")+
    scale_x_discrete(position = "top")+
    ggtitle(paste0("Top ", topGOs.number,u," GOs per gene cluster"))+
    theme(axis.text.y = ,axis.text.x = element_text(angle=45,hjust=1) ,legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5,size = 10),legend.key.size = unit(0.3, "cm"))
}

#----function to plot smoothener of gene expression over pseudotime-----
plot.expression.over.pseudotime <-function(cells.in.path,gene,only.curves=F){
  #' @param cells.in.path Seuratobject containing the pseudotime_values as metadata and reduced to cell that have acctual values, not NA
  #' @param gene Gene to be plotted over pseudotime. Uses nbGAM for the fitting
  #' @param stim whether or not to plot separate lines for the conditions in stim. only for 3 condition implemented
  #' @param only.curves whether ot not to only plot the lines or also the gene expression points per cell
  require(mgcv)
  #rescale gene expression without centering
  cells.in.path.rescaled <-ScaleData(cells.in.path,do.center = F,features = gene,verbose = F)
  #get pseudotime, gene expression and conditions
  data <-FetchData(cells.in.path.rescaled,vars = c("pseudotime_value","stim",gene),slot = "scale.data")
  names(data) <- c("Pseudotime","Condition","gene")
  
  if (only.curves) {
    ggplot(data, aes(x=Pseudotime, y= gene,color=Condition))+
      geom_smooth(method = "gam", formula = y ~ s(x), method.args = list(family = "nb"),size=2,fill = "#D4D4D4")+
      theme(strip.background = element_rect(colour = "white",fill = "white"),panel.border = element_blank(),axis.line.x = element_line(size = 0.25, color = "black"), 
            axis.line.y = element_line(size = 0.25, color = "black"),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.background = element_rect(fill = "white"),legend.key = element_blank())+ 
      scale_color_manual(values=c("Sham"=stim.col.light[1],"TAC"=stim.col.light[2],"TAC_28"=stim.col.light[3]))+
      ylab("Gene expression")+
      ggtitle(gene)
  }else{
    ggplot(data, aes(x=Pseudotime, y= gene,color=Condition))+
      geom_point(size=1, position = position_jitter(0.5))+
      geom_smooth(method = "gam", formula = y ~ s(x), method.args = list(family = "nb"),size=2,fill = "#D4D4D4")+
      theme(strip.background = element_rect(colour = "white",fill = "white"),panel.border = element_blank(),axis.line.x = element_line(size = 0.25, color = "black"), 
            axis.line.y = element_line(size = 0.25, color = "black"),panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.background = element_rect(fill = "white"),legend.key = element_blank())+ 
      scale_color_manual(values=c("Sham"=stim.col.light[1],"TAC"=stim.col.light[2],"TAC_28"=stim.col.light[3]))+
      ylab("Gene expression")+
      ggtitle(gene)
  }#end of else
}

#----function to change the row labels of pheatmap to italic style with labels_row = italic_heatmap_labels(t(expression.values))------
italic_heatmap_labels <- function(data){
  newnames <- lapply(
    rownames(data),
    function(x) bquote(italic(.(x))))
  return(as.expression(newnames))
}