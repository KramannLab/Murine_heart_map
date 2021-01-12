#source: https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
#source: https://f1000research.com/articles/9-709
#-----------DE testing coupled with GO enrichment testing----------------
library(writexl)
library(Seurat)
library(ggplot2)
library(gprofiler2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(Vennerable)
library(stringr)
require(svglite)

#################################################################################################################
#############---code in the following lines need to be adjusted before running---################################
#################################################################################################################
#the input seurat object needs to be an integrated dataset with clustering done.
#the experimental conditions need to be saved unter $stim

#set working directory to the output folder of your data
setwd("C:/Users/fpeis/Desktop/R Stuff/annotation for heart perivascular map/analysis on omen/cdh5 pairwise/gProfileR/")

#enter name of the integrated dataset (e.g. integration_gli1_Sham_TAC) 
sample="full_int_fibroblast_harmony"
control_sample ="Sham" #same as saved in $stim
treatment_sample ="TAC" #same as saved in $stim

#in case of rename of stim is requiered
table(Samples.combined$stim)
Samples.combined$stim[Samples.combined$stim=="FP21_Cdh5_Sham"] <- "Sham"
Samples.combined$stim[Samples.combined$stim=="FP20_Cdh5_TAC"] <- "TAC"

#load the integrated seurat object to the variable subset.cluster
subset.cluster <- Samples.combined #if its loaded already in Samples.combined
rm(Samples.combined)
subset.cluster <- readRDS(file = "C:/Users/fpeis/Desktop/R Stuff/ XXXXX.rds") #if its an rds file

#choose DE test (I recommend MAST, its faster)
test = "DESeq2"
test = "MAST"

#how many top GOs should be included the final plot (ploting is more optimized for 10)
topGOs.number = 10
#how many marker genes should be used for GO enrichment (lower number might give more specific GOs)
marker.genes.nr = 30 #number of genes
marker.genes.nr = "ALL" #incase that all genes should be used

#the lower part of the script tests GOs that are different between the clusters

#################################################################################################################
#############---------------------------------end--------------------------------################################
######################----run everything below here in one go--##################################################
  
  #--two functions that help with visualisation to find nice width and heigth by calculation
  perfect.width <- function(gene.list = data.frame(),coord_fliped = logical()){
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
  perfect.height <- function(gene.list = data.frame(),coord_fliped = logical()){
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
  
  ##############################################################################################################
  #----loop for DE analysis between control_condition and treatment_condition per cluster with MAST or DESeq2(which is slow af)  ----------
  ##############################################################################################################
  #create a new column for joining cluster and condition in order to compare the conditions later within one cluster
  subset.cluster$cluster.stim <- paste(subset.cluster$seurat_clusters, subset.cluster$stim, sep = "_")
  Idents(subset.cluster) <- "cluster.stim"
  dir.create("gprofiler2_out_put")
  dir.create(paste0("gprofiler2_out_put/DE_with_",test))
  dir.create(paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample))
  dir.create(paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/plot_from_gProfiler2"))
  dir.create(paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg"))
  for (i in levels(subset.cluster@meta.data$seurat_clusters)) {
    #DE testing + saving list and plot
    cluster.nr = i
    
    #make sure there are enough cells in the cluster to do DE
    cellnumber.control = length(subset.cluster$cluster.stim[subset.cluster$cluster.stim==paste0(cluster.nr,"_",control_sample)])
    cellnumber.treatment = length(subset.cluster$cluster.stim[subset.cluster$cluster.stim==paste0(cluster.nr,"_",treatment_sample)])
    
    if (cellnumber.control > 10 & cellnumber.treatment > 10) {
      stim.response.new <- FindMarkers(subset.cluster,assay = "RNA",test.use = test,min.pct=0.25,logfc.threshold = 0.3 ,ident.1 = paste0(cluster.nr,"_",treatment_sample), ident.2 = paste0(cluster.nr,"_",control_sample), verbose = FALSE)
      stim.response.new <- arrange(stim.response.new,-avg_logFC)
      stim.response.new$gene_name=rownames(stim.response.new)
      stim.response.new$cluster = cluster.nr
      
      #collect all results in a dataframe of stim.response
      if (cluster.nr == 0) {
        message(paste0("cluster: ",cluster.nr," DE result now saved"))
        stim.response <- stim.response.new
        
      } else {
        message(paste0("cluster: ",cluster.nr," DE result for rbind"))
        stim.response <- rbind(stim.response,stim.response.new)
      }
    
      #save the results of the DE testing itself
      write_xlsx(stim.response.new,path = paste0("gprofiler2_out_put/DE_with_",test,"/",sample,"_DE_analysis_with_",test,"_for_cluster_",cluster.nr,".xlsx"), col_names = TRUE)
      top10=row.names(head(stim.response.new, n = 10))
      top10_low10=c(top10,row.names(tail(stim.response.new, n = 10)))
      Idents(subset.cluster) <- "seurat_clusters"
      VlnPlot(subset(subset.cluster,idents = cluster.nr),assay = "RNA", features =top10_low10,split.by = "stim", pt.size = 0.1,col=c('#00C5C0','#FA756C')) +RestoreLegend()
      ggsave(filename = paste0("gprofiler2_out_put/DE_with_",test,"/",sample,"DE_with_",test,"_markers_for_cluster_",cluster.nr,".jpeg"), width=15 , height =15)
      Idents(subset.cluster) <- "cluster.stim"
      
      #GO of genes that are significant up and down and both combined
      genes_up = stim.response.new$gene_name[stim.response.new$avg_logFC > 0 & stim.response.new$p_val_adj < 0.01]
      genes_down = stim.response.new$gene_name[stim.response.new$avg_logFC < 0 & stim.response.new$p_val_adj < 0.01]
      stim.response.new <- arrange(stim.response.new,p_val)
      genes_all = stim.response.new$gene_name[stim.response.new$p_val_adj < 0.01]
      
      #exclude Ribosomal genes form GO analysis
      genes_up = setdiff(genes_up,genes_up %>%  str_subset(pattern = "^Rp"))
      genes_down = setdiff(genes_down,genes_down %>%  str_subset(pattern = "^Rp"))
      genes_all = setdiff(genes_all,genes_all %>%  str_subset(pattern = "^Rp"))
      
      if (length(genes_up)<5) {
        message(paste0("cluster: ",cluster.nr," has less then 5 up DE genes, therefore no GO testing"))
        go_up=NULL
      } else {
      #run gProfiler for genes that are UP in treatment_condition per cluster
      go_up <- gost(query = genes_up, 
                    organism = "mmusculus", ordered_query = T, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE) }
      
      if (length(genes_down)<5) {
        message(paste0("cluster: ",cluster.nr," has less then 5 down DE genes, therefore no GO testing"))
        go_down=NULL
      } else {
      #run gProfiler for genes that are DOWN in treatment_condition per cluster
      go_down <- gost(query = rev(genes_down), 
                  organism = "mmusculus", ordered_query = T, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL, 
                  numeric_ns = "", sources = NULL, as_short_link = FALSE)}
      
      if (length(genes_all)<5) {
        message(paste0("cluster: ",cluster.nr," has less then 5 DE genes, therefore no GO testing"))
        go_all=NULL
      } else {
      #run gProfiler for genes that are DE overall in treatment_condition per cluster
      go_all <- gost(query = genes_all, 
                      organism = "mmusculus", ordered_query = T, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = NULL, as_short_link = FALSE)}
      
    
      #transfer to temporary dataframe and add colum to save cluster number and UP / Down regulation or from all genes
      if (is.null(go_up)) {
        message(paste0("for cluster: ", cluster.nr," no up GOs found"))
        go_up.new.results=NULL
      } else {
      go_up.new.results = go_up$result
      go_up.new.results$cluster = as.numeric(cluster.nr)
      go_up.new.results$up_or_down = "up"
      #plot with gProfiler plot function the manhattan plot and the top10 GOs
      go_up$result <- arrange(go_up$result,p_value)
      p<-gostplot(go_up, capped = TRUE, interactive =F)
      publish_gostplot(p, highlight_terms = go_up$result$term_id[1:10] ,width = 10, height = 12, filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/","plot_from_gProfiler2/Top10 GOs from DE genes UPregulated in ",treatment_sample," in cluster_",cluster.nr,"_gProfileR2.jpeg"))}
      
      if (is.null(go_down)) {
        message(paste0("for cluster: ", cluster.nr," no down GOs found"))
        go_down.new.results=NULL
      } else {
      go_down.new.results = go_down$result
      go_down.new.results$cluster = as.numeric(cluster.nr)
      go_down.new.results$up_or_down = "down"
      #plot with gProfiler plot function the manhattan plot and the top10 GOs
      go_down$result <- arrange(go_down$result,p_value)
      p2<-gostplot(go_down, capped = TRUE, interactive =F)
      publish_gostplot(p2, highlight_terms = go_down$result$term_id[1:10] ,width = 10, height = 12, filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/","plot_from_gProfiler2/Top10 GOs from DE genes DOWNregulated in ",treatment_sample," in cluster_",cluster.nr,"_gProfileR2.jpeg"))}
    
      if (is.null(go_all)) {
        message(paste0("for cluster: ", cluster.nr," no all GOs found"))
        go_all.new.results=NULL
      } else {
      go_all.new.results = go_all$result
      go_all.new.results$cluster = as.numeric(cluster.nr)
      go_all.new.results$up_or_down = "all"
      #plot with gProfiler plot function the manhattan plot and the top10 GOs
      go_all$result <- arrange(go_all$result,p_value)
      p3<-gostplot(go_all, capped = TRUE, interactive =F)
      publish_gostplot(p2, highlight_terms = go_all$result$term_id[1:10] ,width = 10, height = 12, filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/","plot_from_gProfiler2/Top10 GOs from all DE genes regulated in ",treatment_sample," in cluster_",cluster.nr,"_gProfileR2.jpeg"))}
      
      #add results to overall dataframe
      if (cluster.nr == 0) {
        message(paste0("cluster: ",cluster.nr," now saved"))
        go.results <- rbind(go_up.new.results,go_down.new.results,go_all.new.results)
      
        } else {
        message(paste0("cluster: ",cluster.nr," else for rbind"))
        go.new.results <- rbind(go_up.new.results,go_down.new.results,go_all.new.results)
        go.results <- rbind(go.results,go.new.results)
        }
    }else{
    message(paste("Custer_",cluster.nr,"has only ",cellnumber.control," cell for control and ",cellnumber.treatment," cells for treatment"))  
    }
  }
  #save overall GO result
  saveRDS(go.results,file = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"_test:",test,"_goProfileR2_",control_sample,"_vs_",treatment_sample,".rds"))
  
  #save overall DE result
  saveRDS(stim.response,file = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"_test:",test,"_FindMarkers_",control_sample,"_vs_",treatment_sample,".rds"))
  write_xlsx(stim.response,path = paste0("gprofiler2_out_put/DE_with_",test,"/",sample,"_DE_analysis_with_",test,"_all_DE_genes.xlsx"), col_names = TRUE)
  
  
  ####-----------GOs based only on the UNIQUE  genes-------------###########################################################################
  #generate a dataframe of DE significant genes per cluster keeping only the significant unique DE genes per cluster
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
  
  for (cluster.nr in levels(factor(stim.response$cluster))){
    if (cluster.nr == 0) {go.results.unique <- NULL}
    
    #testing for GOs that are unique per clster -> here are already just significant genes included
    genes_up.unique = genes.de.up.unique$gene_name[genes.de.up.unique$cluster==cluster.nr]
    genes_down.unique = genes.de.down.unique$gene_name[genes.de.down.unique$cluster==cluster.nr]
    genes.de.all.unique <- arrange(genes.de.all.unique,p_val)
    genes_all.unique = genes.de.all.unique$gene_name[genes.de.all.unique$cluster==cluster.nr]
    
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
  
  #save overall unique DE list
  write_xlsx(genes.de.up.unique,path = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"_test:",test,"_DE_genes_unique_for_a_cluster.xlsx"), col_names = TRUE)
  saveRDS(genes.de.up.unique,file = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"_test:",test,"_DE_genes_unique_for_a_cluster.rds"))
  write_xlsx(go.results.unique,path = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number,"_GOs_unique_for_a_cluster.xlsx"), col_names = TRUE)
  saveRDS(go.results.unique,file = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"_test:",test,"_GOs_unique_for_a_cluster.rds"))
  
  ####-----------GOs based only on the SHARED genes and adding to common list with UNIQUE-------###########################################
  #generate a dataframe of DE genes per cluster keeping only the unique DE genes per cluster
  gene.list.up.shared= setdiff(unique(stim.response$gene_name[stim.response$avg_logFC>0& stim.response$p_val_adj < 0.01]),genes.de.up.unique$gene_name)
  gene.list.down.shared= setdiff(unique(stim.response$gene_name[stim.response$avg_logFC<0& stim.response$p_val_adj < 0.01]),genes.de.down.unique$gene_name)
  stim.response <- arrange(stim.response,p_val)
  gene.list.all.shared= setdiff(unique(stim.response$gene_name[stim.response$p_val_adj < 0.01]),genes.de.all.unique$gene_name)
  
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
  
  #########################################################################################################################
  #########################----------------start of plotting GOs from control vs. treatment DE genes---####################
  #########################################################################################################################
  
  #########################--GOs based on only UPregulated genes---####
  #----------------generate a dataframe of only the top10 GOs per cluster---------------
  go.results.up <- go.results[go.results$up_or_down =="up",]
  go.results.up <- go.results.up[go.results.up$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results.up <- go.results.up[go.results.up$term_size<1000,] #filter out large and often very general GOs
  
  go.results.up <- arrange(go.results.up,p_value)
  go.results.up$cluster <-as.numeric(go.results.up$cluster)
  go.results.up <- arrange(go.results.up,cluster)
  
  #save top 50 per cluster
  go.results.up.top50 <- go.results.up %>% group_by(cluster) %>% top_n(n = -50,wt=p_value)
  write_xlsx(go.results.up.top50,path = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number,"_top50_UP_GOs.xlsx"), col_names = TRUE)
  
  #------extract and plot the topXX by p_value per cluster
  go.results.up.top10 <- go.results.up %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  order_desc_up <- unique(go.results.up.top10$term_name)
  go.results.up.top10.all <- go.results.up[go.results.up$term_name %in% order_desc_up,]
  
  #-plot up
  ggplot(go.results.up.top10.all, aes(x=factor(term_name, level = order_desc_up), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs upregulated in ",treatment_sample," per cluster based upregulated DE genes"))+
    theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number," GOs UPregulated in ",treatment_sample," on UP DE.jpeg"), width=perfect.width(go.results.up.top10.all,T) , height = perfect.height(go.results.up.top10.all,T))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top",topGOs.number," GOs UPregulated in ",treatment_sample," on UP DE.svg"), width=perfect.width(go.results.up.top10.all,T) , height = perfect.height(go.results.up.top10.all,T))
  
  #------extract and plot the top5 (p_value) per scource per cluster
  for (cluster in levels(factor(go.results.up$cluster))) {
    if (cluster == 0) {
      go.results.up.top5.perSource <- go.results.up[go.results.up$cluster==cluster,]  %>% group_by(source) %>% top_n(n = -5,wt=p_value)
    } else{
      go.results.up.top5.perSource.new <- go.results.up[go.results.up$cluster==cluster,]  %>% group_by(source) %>% top_n(n = -5,wt=p_value)  
      go.results.up.top5.perSource <-rbind(go.results.up.top5.perSource,go.results.up.top5.perSource.new)
    }
  }
  #plot with ordering after source
  go.results.up.top5.perSource.all <- go.results.up[go.results.up$term_name %in% unique(go.results.up.top5.perSource$term_name),]
  go.results.up.top5.perSource.all <- arrange(go.results.up.top5.perSource.all,source)
  order_desc_up <- unique(go.results.up.top5.perSource.all$term_name)
  
  ggplot(go.results.up.top5.perSource.all, aes(x=factor(term_name, level = order_desc_up), y=factor(as.numeric(cluster)),size=-log10(p_value),fill=source,color=source)) + 
    geom_point() + 
    scale_size(range = c(1,5)) +
    ylab("Cluster") + xlab("") + ggtitle(paste0("GOs upregulated in ",treatment_sample," based on upregulated DE genes"))+
    theme(legend.title = element_text(size = 8),axis.text.x = element_text(angle = 90,hjust = 1),plot.title.position = "plot",plot.title = element_text(hjust = 0.5),legend.key.size = unit(0.2, "cm"))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample," top5 GOs per database based on UPregulated genes in ",treatment_sample,"_col.jpeg"), width=perfect.width(go.results.up.top5.perSource.all,F) , height = perfect.height(go.results.up.top5.perSource.all,coord_fliped = F))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample," top5 GOs per database based on UPregulated genes in ",treatment_sample,"_col.svg"), width=perfect.width(go.results.up.top5.perSource.all,F) , height = perfect.height(go.results.up.top5.perSource.all,coord_fliped = F))
  
  #plot with ordering on clusternumber
  go.results.up.top5.perSource.all <- go.results.up[go.results.up$term_name %in% unique(go.results.up.top5.perSource$term_name),]
  go.results.up.top5.perSource.all <- arrange(go.results.up.top5.perSource.all,source)
  go.results.up.top5.perSource.all$go_and_source=paste0(go.results.up.top5.perSource.all$term_name,"_",go.results.up.top5.perSource.all$source)
  order_desc_up <- unique(go.results.up.top5.perSource.all$go_and_source)
  ggplot(go.results.up.top5.perSource.all, aes(x=factor(go_and_source, level = order_desc_up), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(1,5)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("GOs upregulated in ",treatment_sample," based on upregulated DE genes"))+
    theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5),legend.key.size = unit(0.2, "cm"),legend.box.just = "left")
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample," top5 GOs per database based on UPregulated genes in ",treatment_sample,".jpeg"), width=perfect.width(go.results.up.top5.perSource.all,T) , height = perfect.height(go.results.up.top5.perSource.all,T))
  #ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample," top5 GOs per database based on UPregulated genes in ",treatment_sample,".jpeg"), width=6 , height = 3)
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample," top5 GOs per database based on UPregulated genes in ",treatment_sample,".svg"), width=perfect.width(go.results.up.top5.perSource.all,T) , height = perfect.height(go.results.up.top5.perSource.all,T))
  
  #########################--GOs based on only DOWNregulated genes---####
  #----------------generate a dataframe of only the top10 GOs per cluster---------------
  go.results.down <- go.results[go.results$up_or_down =="down",]
  go.results.down <- go.results.down[go.results.down$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results.down <- go.results.down[go.results.down$term_size<1000,] #filter out large and often very general GOs
  
  go.results.down <- arrange(go.results.down,p_value)
  go.results.down$cluster <- as.numeric(go.results.down$cluster)
  go.results.down <- arrange(go.results.down,cluster)
  
  #------extract and plot the topXX by p_value per cluster
  go.results.down.top10 <- go.results.down %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  order_desc_down <- unique(go.results.down.top10$term_name)
  go.results.down.top10.all <- go.results.down[go.results.down$term_name %in% order_desc_down,]
  
  #-plot down
  ggplot(go.results.down.top10.all, aes(x=factor(term_name, level = order_desc_down), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs downregulated in ",treatment_sample," per cluster based downregulated DE genes"))+
    theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number," GOs DOWNregulated in ",treatment_sample," on DOWN DE.jpeg"), width=perfect.width(go.results.down.top10,T) , height = perfect.height(go.results.down.top10,T))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top",topGOs.number," GOs DOWNregulated in ",treatment_sample," on DOWN DE.svg"), width=perfect.width(go.results.down.top10,T) , height = perfect.height(go.results.down.top10,T))
  
  #------extract and plot the top5 (p_value) per scource per cluster
  for (cluster in levels(factor(go.results.down$cluster))) {
    if (cluster == 0) {
      go.results.down.top5.perSource <- go.results.down[go.results.down$cluster==cluster,]  %>% group_by(source) %>% top_n(n = -5,wt=p_value)
    } else{
      go.results.down.top5.perSource.new <- go.results.down[go.results.down$cluster==cluster,]  %>% group_by(source) %>% top_n(n = -5,wt=p_value)  
      go.results.down.top5.perSource <-rbind(go.results.down.top5.perSource,go.results.down.top5.perSource.new)
    }
  }
  #plot with ordering after source
  go.results.down.top5.perSource.all <- go.results.down[go.results.down$term_name %in% unique(go.results.down.top5.perSource$term_name),]
  go.results.down.top5.perSource.all <- arrange(go.results.down.top5.perSource.all,source)
  order_desc_down <- unique(go.results.down.top5.perSource.all$term_name)
  
  ggplot(go.results.down.top5.perSource.all, aes(x=factor(term_name, level = order_desc_down), y=factor(as.numeric(cluster)),size=-log10(p_value),fill=source,color=source)) + 
    geom_point() + 
    scale_size(range = c(1,5)) +
    ylab("Cluster") + xlab("") + ggtitle(paste0("GOs downregulated in ",treatment_sample," based on downregulated DE genes"))+
    theme(legend.title = element_text(size = 8),axis.text.x = element_text(angle = 90,hjust = 1),plot.title.position = "plot",plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top5 GOs per database based on DOWNregulated genes in ",treatment_sample,"_col.jpeg"), width=perfect.width(go.results.down.top5.perSource,F) , height = perfect.height(go.results.down.top5.perSource,coord_fliped = F))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top5 GOs per database based on DOWNregulated genes in ",treatment_sample,"_col.svg"), width=perfect.width(go.results.down.top5.perSource,F) , height = perfect.height(go.results.down.top5.perSource,coord_fliped = F))
  
  #plot with ordering on clusternumber
  go.results.down.top5.perSource.all <- go.results.down[go.results.down$term_name %in% unique(go.results.down.top5.perSource$term_name),]
  go.results.down.top5.perSource.all <- arrange(go.results.down.top5.perSource.all,source)
  go.results.down.top5.perSource.all$go_and_source=paste0(go.results.down.top5.perSource.all$term_name,"_",go.results.down.top5.perSource.all$source)
  order_desc_down <- unique(go.results.down.top5.perSource.all$go_and_source)
  ggplot(go.results.down.top5.perSource.all, aes(x=factor(go_and_source, level = order_desc_down), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(1,5)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("GOs downregulated in ",treatment_sample," based on downregulated DE genes"))+
    theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top5 GOs per database based on DOWNregulated genes in ",treatment_sample,".jpeg"), width=perfect.width(go.results.down.top5.perSource,T) , height = perfect.height(go.results.down.top5.perSource,T))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top5 GOs per database based on DOWNregulated genes in ",treatment_sample,".svg"), width=perfect.width(go.results.down.top5.perSource,T) , height = perfect.height(go.results.down.top5.perSource,T))
  
  #########################--GOs based on SHARED + UNIQUE DE genes per cluster-------########################################
  #####--------plotting GOs based on only UP unique genes-----------------
  go.results.unique.filtered <- go_up.shared.and.unique.results[go_up.shared.and.unique.results$up_or_down =="up",]
  go.results.unique.filtered <- go.results.unique.filtered[go.results.unique.filtered$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results.unique.filtered <- go.results.unique.filtered[go.results.unique.filtered$term_size<1000,] #filter out large and often very general GOs
  
  #soring was done when joining unique and shared gos
  
  #------extract and plot the topXX by p_value per cluster
  go.results.unique.filtered.top10 <- go.results.unique.filtered %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  go.results.unique.filtered.top10$cluster <- as.numeric(go.results.unique.filtered.top10$cluster)
  go.results.unique.filtered.top10 <- arrange(go.results.unique.filtered.top10,cluster)
  order_desc_up <- unique(go.results.unique.filtered.top10$term_name)
  go.results.unique.filtered.top10.all <- go.results.unique.filtered[go.results.unique.filtered$term_name %in% order_desc_up,]
  go.results.unique.filtered.top10.all$cluster <- factor(go.results.unique.filtered.top10.all$cluster,levels = c(levels(factor(go.results.unique$cluster)),"shared"))
  
  ggplot(go.results.unique.filtered.top10.all, aes(x=factor(term_name, level = rev(order_desc_up)), y=cluster,size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs upregulated in ",treatment_sample," per cluster based unique + shared upregulated DE genes"))+
    theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0,size = 10),legend.key.size = unit(0.2, "cm"),legend.box.just = "left") + scale_x_discrete(position = "top") #+theme(legend.position = c(5.5, 0.28))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.jpeg"), width=perfect.width(go.results.unique.filtered.top10.all,T)+1 , height = perfect.height(go.results.unique.filtered.top10.all,T))
  #ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.jpeg"), width=6.5 , height = 6)
  
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.svg"), width=perfect.width(go.results.unique.filtered.top10.all,T)+1 , height = perfect.height(go.results.unique.filtered.top10.all,T))
  #ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top",topGOs.number,"GOs upregulated in ",treatment_sample," based on unique + shared UP DE genes.svg"), width=6.5 , height = 6)
  
  ######--------plotting GOs based on only ALL unique genes-----------------
  go.results.unique.filtered <- go_up.shared.and.unique.results[go_up.shared.and.unique.results$up_or_down =="all",]
  go.results.unique.filtered <- go.results.unique.filtered[go.results.unique.filtered$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results.unique.filtered <- go.results.unique.filtered[go.results.unique.filtered$term_size<1000,] #filter out large and often very general GOs
  
  #soring was done when joining unique and shared gos
  
  #------extract and plot the topXX by p_value per cluster
  go.results.unique.filtered.top10 <- go.results.unique.filtered %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  order_desc_up <- unique(go.results.unique.filtered.top10$term_name)
  go.results.unique.filtered.top10.all <- go.results.unique.filtered[go.results.unique.filtered$term_name %in% order_desc_up,]
  go.results.unique.filtered.top10.all$cluster <- factor(go.results.unique.filtered.top10.all$cluster,levels = c(levels(factor(go.results.unique$cluster)),"shared"))
  
  ggplot(go.results.unique.filtered.top10.all, aes(x=factor(term_name, level = rev(order_desc_up)), y=cluster,size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs upregulated in ",treatment_sample," based on unique + shared overall DE genes"))+
    theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0,size = 10),legend.key.size = unit(0.2, "cm"),legend.box.just = "left") + scale_x_discrete(position = "top") #+theme(legend.position = c(5.5, 0.28))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number," GOs upregulated in ",treatment_sample," based on unique + shared overall DE genes.jpeg"), width=perfect.width(go.results.unique.filtered.top10.all,T)+1 , height = perfect.height(go.results.unique.filtered.top10.all,T))
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/",sample,"top",topGOs.number," GOs upregulated in ",treatment_sample," based on unique + shared overall DE genes.jpeg"), width=6 , height = 3.3)
  
  ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top",topGOs.number," GOs upregulated in ",treatment_sample," based on unique + shared overall DE genes.svg"), width=perfect.width(go.results.unique.filtered.top10.all,T)+1 , height = perfect.height(go.results.unique.filtered.top10.all,T))
 # ggsave(filename = paste0("gprofiler2_out_put/",control_sample,"_vs_",treatment_sample,"/svg/",sample,"top",topGOs.number," GOs upregulated in ",treatment_sample," based on unique + shared overall DE genes.svg"), width=6.2 , height = 4.2)

  ############################################################################################################################################
  ###################---------test GOs per cluster based on FindAllMarkers--------------#####################################################
  ###########################################################################################################################################
  Idents(subset.cluster)<-subset.cluster$seurat_clusters
  all.markers = FindAllMarkers(subset.cluster, test.use = "MAST",min.pct = 0.3,logfc.threshold = 0.3,assay = "RNA")
  dir.create("gprofiler2_out_put/comparison_between_clusters")
  dir.create("gprofiler2_out_put/comparison_between_clusters/plot_from_gProfiler2")
  dir.create("gprofiler2_out_put/comparison_between_clusters/svg")
  saveRDS(all.markers,file = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"_all.markers_markergenes_per_cluster_result.RDS"))
  #all.markers <-readRDS(file = "~/Documents/FP_scRNA/R stuff/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated/res0.3/gprofiler2_out_put/comparison_between_clusters/Fibroblast_Col1a1_Gli1_Pdgfrb_integrated__all.markers_markergenes_per_cluster_result.RDS")
  all.markers <- arrange(all.markers,avg_logFC)
  all.markers <- arrange(all.markers,cluster)
  
  for (i in levels(subset.cluster@meta.data$seurat_clusters)) {
    cluster.nr = i
    #GO of marker genes that are up and down
    genes_up = all.markers$gene[all.markers$avg_logFC>0 & all.markers$cluster==cluster.nr & all.markers$p_val_adj < 0.01]
    genes_down = all.markers$gene[all.markers$avg_logFC<0 & all.markers$cluster==cluster.nr& all.markers$p_val_adj < 0.01]
    all.markers <- arrange(all.markers,p_val)
    genes_all = all.markers$gene[all.markers$cluster==cluster.nr& all.markers$p_val_adj < 0.01]
    
  
    #optional filter for only topXX markergenes to apply GO enrichment 
    if (marker.genes.nr != "ALL") {
      genes_up <- genes_up[1:marker.genes.nr]
      genes_down <- genes_down[1:marker.genes.nr]
      genes_all <- genes_all[1:marker.genes.nr]
    }
    #exclude Ribosomal genes form GO analysis
    genes_up = setdiff(genes_up,genes_up %>%  str_subset(pattern = "^Rp"))
    genes_down = setdiff(genes_down,genes_down %>%  str_subset(pattern = "^Rp"))
    genes_all = setdiff(genes_all,genes_all %>%  str_subset(pattern = "^Rp"))
    
    if (length(genes_up)<5) {
      message(paste0("cluster: ",cluster.nr," has less then 5 up DE genes, therefore no GO testing"))
      go_up=NULL
    } else {
      #run gProfiler for genes that are UP in treatment_condition per cluster
      go_up <- gost(query = genes_up, 
                    organism = "mmusculus", ordered_query = T, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE) }
    
    if (length(genes_down)<5) {
      message(paste0("cluster: ",cluster.nr," has less then 5 down DE genes, therefore no GO testing"))
      go_down=NULL
    } else {
      #run gProfiler for genes that are DOWN in treatment_condition per cluster
      go_down <- gost(query = rev(genes_down), 
                      organism = "mmusculus", ordered_query = T, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = NULL, as_short_link = FALSE)}
    
    if (length(genes_all)<5) {
      message(paste0("cluster: ",cluster.nr," has less then 5 DE genes, therefore no GO testing"))
      go_all=NULL
    } else {
      #run gProfiler for genes that are DE overall in treatment_condition per cluster
      go_all <- gost(query = genes_all, 
                     organism = "mmusculus", ordered_query = T, 
                     multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                     measure_underrepresentation = FALSE, evcodes = FALSE, 
                     user_threshold = 0.05, correction_method = "g_SCS", 
                     domain_scope = "annotated", custom_bg = NULL, 
                     numeric_ns = "", sources = NULL, as_short_link = FALSE)}
    
    #transfer to temporary dataframe and add column to save cluster number and UP or Down regulation
    if (is.null(go_up)) {
      message(paste0("for cluster: ", cluster.nr," no up GOs found"))
      go_up.new.results=NULL
    } else {
      go_up.new.results = go_up$result
      go_up.new.results$cluster = as.numeric(cluster.nr)
      go_up.new.results$up_or_down = "up"}
    if (is.null(go_down)) {
      message(paste0("for cluster: ", cluster.nr," no down GOs found"))
      go_down.new.results=NULL
    } else {
      go_down.new.results = go_down$result
      go_down.new.results$cluster = as.numeric(cluster.nr)
      go_down.new.results$up_or_down = "down"}
    if (is.null(go_all)) {
      message(paste0("for cluster: ", cluster.nr," no GOs found"))
      go_all.new.results=NULL
    } else {
      go_all.new.results = go_all$result
      go_all.new.results$cluster = as.numeric(cluster.nr)
      go_all.new.results$up_or_down = "all"}
    
    #add results to overall dataframe
    if (cluster.nr == 0) {
      message(paste0("cluster: ",cluster.nr," now saved"))
      go.results <- rbind(go_up.new.results,go_down.new.results,go_all.new.results)
      
    } else {
      message(paste0("cluster: ",cluster.nr," else for rbind"))
      go.new.results <- rbind(go_up.new.results,go_down.new.results,go_all.new.results)
      go.results <- rbind(go.results,go.new.results)
    }
    
    #plot with gProfiler plot function the manhattan plot and the top10 GOs 
    if (is.null(go_up)){}else{ go_up$result <- arrange(go_up$result,p_value)
    p<-gostplot(go_up, capped = TRUE, interactive =F)
    publish_gostplot(p, highlight_terms = go_up$result$term_id[1:10] ,width = 10, height = 12, filename = paste0("gprofiler2_out_put/comparison_between_clusters/plot_from_gProfiler2/Top10 GOs for UP markergenes of cluster_",cluster.nr,"_gProfileR2.jpeg"))}
    if (is.null(go_down)){}else{go_down$result <- arrange(go_down$result,p_value)
    p2<-gostplot(go_down, capped = TRUE, interactive =F)
    publish_gostplot(p2, highlight_terms = go_down$result$term_id[1:10] ,width = 10, height = 12, filename = paste0("gprofiler2_out_put/comparison_between_clusters/plot_from_gProfiler2/Top10 GOs for DOWN markergenes of cluster_",cluster.nr,"_gProfileR2.jpeg"))}
    if (is.null(go_all)){}else{go_all$result <- arrange(go_all$result,p_value)
    p3<-gostplot(go_all, capped = TRUE, interactive =F)
    publish_gostplot(p, highlight_terms = go_all$result$term_id[1:10] ,width = 10, height = 12, filename = paste0("gprofiler2_out_put/comparison_between_clusters/plot_from_gProfiler2/Top10 GOs for all markergenes of cluster_",cluster.nr,"_gProfileR2.jpeg"))}
    }
  #save overall result
  saveRDS(go.results,file = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"GOs per cluster based on_",marker.genes.nr,"_FindAllMarkers_result.rds"))
  write_xlsx(go.results,path = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"GOs per cluster based on_",marker.genes.nr,"_FindAllMarkers_result_all_GOs.xlsx"), col_names = TRUE)
  
  #########################--up marker genes plotting---########################################################################
  #----------------generate a dataframe of only the topXX GOs per cluster---------------
  go.results.up <- go.results[go.results$up_or_down =="up",]
  go.results.up <- go.results.up[go.results.up$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results.up <- go.results.up[go.results.up$term_size<1000,] #filter out large and often very general GOs
  
  go.results.up <- arrange(go.results.up,p_value)
  go.results.up <- arrange(go.results.up,cluster)
  
  #save top 50 and top 20 per cluster
  go.results.up.top50 <- go.results.up %>% group_by(cluster) %>% top_n(n = -50,wt=p_value)
  write_xlsx(go.results.up.top50,path = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"nMarker_",marker.genes.nr,"_top50_UP_GOs.xlsx"), col_names = TRUE)
  go.results.up.top20 <- go.results.up %>% group_by(cluster) %>% top_n(n = -20,wt=p_value)
  write_xlsx(go.results.up.top20,path = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"nMarker_",marker.genes.nr,"_top20_UP_GOs.xlsx"), col_names = TRUE)
   
  #------extract and plot the top10 by p_value per cluster
  go.results.up.top10 <- go.results.up %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  go.results.up.top10$cluster=as.numeric(go.results.up.top10$cluster)
  go.results.up.top10<-arrange(go.results.up.top10,cluster)
  order_desc_up <- unique(go.results.up.top10$term_name)
  #order_desc_up <- setdiff(order_desc_up,longest.term)
  go.results.up.top10.all <- go.results.up[go.results.up$term_name %in% order_desc_up,]
    
  #-plot up vertikal
  ggplot(go.results.up.top10.all, aes(x=factor(term_name, level = order_desc_up), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
    geom_point() + 
    scale_size(range = c(2,6)) +
    coord_flip() +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster"))+
    theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5,size = 10),legend.key.size = unit(0.3, "cm"))
  ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.jpeg"), width=perfect.width(go.results.up.top10.all,T) , height = perfect.height(go.results.up.top10.all,T))
  ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/svg/",sample,"Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.svg"), width=perfect.width(go.results.up.top10.all,T) , height = perfect.height(go.results.up.top10.all,T))
  #ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.jpeg"), width=7.7, height = 4.2)
  #ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/svg/",sample,"Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_v.svg"), width=7.7, height = 4.2)
  
  
  #-plot up horizontal
  ggplot(go.results.up.top10.all, aes(x=factor(term_name, level = order_desc_up), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) +
    geom_point() +
    scale_size(range = c(1,5)) +
    ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," marker genes per cluster")) +
    theme(legend.title = element_text(size = 8),axis.text.x = element_text(angle = 90,hjust = 1),plot.title.position = "plot",plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_h.jpeg"), width=perfect.width(go.results.up.top10.all,F) , height = perfect.height(go.results.up.top10.all,F))
  #ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top",topGOs.number,"GOs of ",marker.genes.nr,"UP marker genes per cluster_v.jpeg"), width=7.7, height = 4.2)
  ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/svg/",sample,"Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," upregulated marker genes per cluster_h.svg"), width=perfect.width(go.results.up.top10.all,F) , height = perfect.height(go.results.up.top10.all,F))
  
  
  #########################--down marker genes plotting---########################################################################
  #not implemented because of the way DE genes are detected this might be missleading
  #########################-----------------------------########################################################################
  
  # #########################--ALL marker genes plotting---########################################################################
  # #not recommended because very unspecific results
  # #----------------generate a dataframe of only the top10 GOs per cluster---------------
  # go.results.all <- go.results[go.results$up_or_down =="all",]
  # go.results.all <- go.results.all[go.results.all$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  # go.results.all <- go.results.all[go.results.all$term_size<1000,] #filter out large and often very general GOs
  # 
  # go.results.all <- arrange(go.results.all,p_value)
  # go.results.all <- arrange(go.results.all,cluster)
  # 
  # #save top 50 and top 20 per cluster
  # go.results.all.top50 <- go.results.all %>% group_by(cluster) %>% top_n(n = -50,wt=p_value)
  # write_xlsx(go.results.all.top50,path = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"nMarker",marker.genes.nr,"_top50_all_GOs.xlsx"), col_names = TRUE)
  # go.results.all.top20 <- go.results.all %>% group_by(cluster) %>% top_n(n = -20,wt=p_value)
  # write_xlsx(go.results.all.top20,path = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"nMarker",marker.genes.nr,"_top20_all_GOs.xlsx"), col_names = TRUE)
  # 
  # #------extract and plot the top10 by p_value per cluster
  # go.results.all.top10 <- go.results.all %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  # go.results.all.top10$cluster=as.numeric(go.results.all.top10$cluster)
  # go.results.all.top10<-arrange(go.results.all.top10,cluster)
  # order_desc_up <- unique(go.results.all.top10$term_name)
  # go.results.all.top10.all <- go.results.all[go.results.all$term_name %in% order_desc_up,]
  # 
  # #-plot all vertikal
  # ggplot(go.results.all.top10.all, aes(x=factor(term_name, level = order_desc_up), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
  #   geom_point() + 
  #   scale_size(range = c(2,6)) +
  #   coord_flip() +
  #   ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," marker genes per cluster"))+
  #   theme(legend.title = element_text(size = 8),plot.title.position = "plot",plot.title = element_text(hjust = 0.5),legend.key.size = unit(0.3, "cm"))
  # ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top",topGOs.number,"GOs of ",marker.genes.nr,"_ALL marker genes per cluster_v.jpeg"), width=perfect.width(go.results.all.top10.all,T) , height = perfect.height(go.results.all.top10.all,T))
  # #ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top",topGOs.number,"GOs of ",marker.genes.nr,"_all marker genes per cluster_v.jpeg"), width=7.7, height = 4.2)
  # ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top",topGOs.number,"GOs of ",marker.genes.nr,"_ALL marker genes per cluster_v.svg"), width=perfect.width(go.results.all.top10.all,T) , height = perfect.height(go.results.all.top10.all,T))
  # 
  # 
  # #-plot all horizontal
  # ggplot(go.results.all.top10.all, aes(x=factor(term_name, level = order_desc_up), y=factor(as.numeric(cluster)),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
  #   geom_point() + 
  #   scale_size(range = c(1,5)) +
  #   ylab("Cluster") + xlab("") + ggtitle(paste0("Top",topGOs.number," GOs per cluster based on ",marker.genes.nr," marker genes per cluster")) +
  #   theme(legend.title = element_text(size = 8),axis.text.x = element_text(angle = 90,hjust = 1),plot.title.position = "plot",plot.title = element_text(hjust = 0.5))
  # ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top",topGOs.number,"GOs of ",marker.genes.nr,"_all marker genes per cluster_h.jpeg"), width=perfect.width(go.results.all.top10.all,F) , height = perfect.height(go.results.all.top10.all,F))
  # ggsave(filename = paste0("gprofiler2_out_put/comparison_between_clusters/",sample,"top",topGOs.number,"GOs of ",marker.genes.nr,"_all marker genes per cluster_h.svg"), width=perfect.width(go.results.all.top10.all,F) , height = perfect.height(go.results.all.top10.all,F))
  
  
  #################################----------- Volcano and heatmaps plot DE genes-----------------##############################
  dir.create("gprofiler2_out_put/volc&heatmap_plots_of_DE")
  dir.create("gprofiler2_out_put/volc&heatmap_plots_of_DE/svg")
  #generate a list from stim.response with only significant DE genes
  stim.response.sig = stim.response[stim.response$p_val_adj < 0.01,]
  
  #---adding to stim.response.sig if the gene is unique for the cluster or not
  #adding to stim.response.sig if the gene is unique for the cluster or not
  for (gene in unique(stim.response.sig$gene_name)) {
    if (gene %in% c(genes.de.up.unique$gene_name,genes.de.down.unique$gene_name)) { #check if unique at all
      if (gene %in% genes.de.up.unique$gene_name) { #check if it is unique in UP
        stim.response.sig$unique[stim.response.sig$gene_name==gene & stim.response.sig$avg_logFC>0] = paste0("unique for cl ",genes.de.up.unique$cluster[genes.de.up.unique$gene_name==gene])
      }else { #if its not in unique UP it must be in unique DOWN
        stim.response.sig$unique[stim.response.sig$gene_name==gene & stim.response.sig$avg_logFC<0] = paste0("unique for cl ",genes.de.down.unique$cluster[genes.de.down.unique$gene_name==gene])
      }
    } else {
      stim.response.sig$unique[stim.response.sig$gene_name==gene] ="shared"
    }
  }
  stim.response.sig$unique[is.na(stim.response.sig$unique)] ="shared" #all NAs in the unique check must be shared in UP but unique in down and wise versa
  
  
  #------plot the volcano plot for all samples
  for (cluster.nr in levels(factor(stim.response.sig$cluster))) {
    cluster = cluster.nr
    library(ggrepel)
    DE.result.cl = stim.response.sig[stim.response.sig$cluster==cluster,]
    DE.result.cl$p_val_adj[DE.result.cl$p_val_adj==0] = 1.000000e-300 #genes with significe of 0 are set to 1.000000e-300
    #look for genes about threshhold (cutoff for avg_logFC is deactivated)(cutoff for the p_val_adj ist dynamically changed for each sample to the man p_val_adj of all samples)
    for (gene in DE.result.cl$gene_name) {
      #if (DE.result.cl$avg_logFC[DE.result.cl$gene_name == gene]>0.5 | DE.result.cl$avg_logFC[DE.result.cl$gene_name == gene]< (-0.5) ) {
        if (-log10(DE.result.cl$p_val_adj[DE.result.cl$gene_name == gene])>(mean(-log10(DE.result.cl$p_val_adj)))) { 
          DE.result.cl$high[DE.result.cl$gene_name == gene] = gene 
          DE.result.cl$overth[DE.result.cl$gene_name == gene] = TRUE
        } else {DE.result.cl$overth[DE.result.cl$gene_name == gene] = FALSE
        DE.result.cl$high[DE.result.cl$gene_name == gene] = ""}
      #} else {
      #  DE.result.cl$high[DE.result.cl$gene_name == gene] = ""
      #  DE.result.cl$overth[DE.result.cl$gene_name == gene] = FALSE
      #}
    }
    #add if alo unique or shared
    for (gene in DE.result.cl$gene_name) {
      if (DE.result.cl$overth[DE.result.cl$gene_name == gene] == FALSE) {
        DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "shared"
        if (DE.result.cl$unique[DE.result.cl$gene_name == gene] != "shared") {DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "unique"}
      } else {
        if (DE.result.cl$unique[DE.result.cl$gene_name == gene] == "shared") {
          DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "shared"
        } else {
          DE.result.cl$high.unique[DE.result.cl$gene_name == gene] = "unique"
        }
      }
    }
    #volcanoplot
    ggplot(DE.result.cl) +
      geom_point(aes(x=avg_logFC, y=-log10(p_val_adj),colour=high.unique),size=2.5) +
      ggtitle(paste0("Cluster ",cluster, " overexpression in ",treatment_sample)) +
      xlab("avg_logFC") + 
      ylab("-log10 adjusted p-value") +
      #scale_y_continuous(limits = c(0,50)) +
      geom_hline(yintercept=(mean(-log10(DE.result.cl$p_val_adj))+1), linetype="dashed", color = "grey")+
      geom_vline(xintercept=0.5, linetype="dashed", color = "grey")+
      geom_vline(xintercept=-0.5, linetype="dashed", color = "grey")+
      theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
      geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = high),size = 5) +
      scale_color_manual(values = c("shared"="blue","unique"="red")) +
      theme(legend.title = element_blank(),legend.position = c(0.95,0.1),plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))  
    ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"volcano_plot_cluster",cluster,".jpeg"), width=8, height =8)
    ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/svg/",sample,"volcano_plot_cluster",cluster,".svg"), width=8, height =8)
    }
  ###------------------generate a heatmap of all genes that were found DE in any cluster-----------------------------
  Idents(subset.cluster) <- "cluster.stim"
  subset.cluster$cluster.stim = factor(subset.cluster$cluster.stim,levels((Idents(subset.cluster)))[order(parse_character(levels(Idents(subset.cluster))))])
  if(length(levels(subset.cluster$seurat_clusters))>10){subset.cluster$cluster.stim = factor(subset.cluster$cluster.stim,levels((Idents(subset.cluster)))[order(parse_number(levels(Idents(subset.cluster))))])}
  Idents(subset.cluster) <- "cluster.stim"
  #generate average expression of genes per cluster and condition
  avgexp = AverageExpression(subset.cluster, return.seurat = T)
  
  stim.response.sig<-arrange(stim.response.sig,avg_logFC)
  stim.response.sig$cluster <-as.numeric(stim.response.sig$cluster)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
  stim.response.sig<-arrange(stim.response.sig,cluster)
  
  #plot the average only for the DE genes UP 
  DoHeatmap(avgexp,features = unique(stim.response.sig$gene_name[stim.response.sig$avg_logFC>0]),draw.lines = F) 
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"DE genes UP average_per_cl_cond.jpeg"), width=10, height = 20)
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/svg/",sample,"DE genes UP average_per_cl_cond.svg"), width=10, height = 20)
  
  #plot the average only for the DE genes DOWN 
  DoHeatmap(avgexp,features = unique(stim.response.sig$gene_name[stim.response.sig$avg_logFC<0]),draw.lines = F)
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"DE genes DOWN average_per_cluster_condition.jpeg"), width=10, height = 20)
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/svg/",sample,"DE genes DOWN average_per_cluster_condition.svg"), width=10, height = 20)
  
  #---------following section for UP DE genes
  #plot the average only for the DE genes UP + label of shared and unique genes
  DE.label.for.hm=NULL
  DE.label.for.hm$gene_name = gene.list.up.shared
  DE.label.for.hm = data.frame(DE.label.for.hm)
  DE.label.for.hm$from = "shared"
  DE.uni.df=NULL
  DE.uni.df$gene_name = genes.de.up.unique$gene_name
  DE.uni.df = data.frame(DE.uni.df)
  DE.uni.df$from = paste0("cluster",genes.de.up.unique$cluster)
  DE.uni.df<-DE.uni.df[order(DE.uni.df$from ),]
  DE.heatmap = rbind(DE.label.for.hm,DE.uni.df)
  #count how often a gene was found UP and add to shared 
  gene.counter = NULL
  for (gene in stim.response.sig$gene_name[stim.response.sig$avg_logFC>0]) {
    gene.counter[[gene]] = sum(stim.response.sig$gene_name[stim.response.sig$avg_logFC>0] == gene)
  }
  gene.counter.df = data.frame("genes"= names(gene.counter),"count"=unlist(gene.counter, use.names=FALSE))
  gene.counter.df = arrange(gene.counter.df,-count)
  for (gene in DE.heatmap$gene_name) {
    if (DE.heatmap$from[DE.heatmap$gene_name==gene] == "shared") {
      DE.heatmap$from[DE.heatmap$gene_name==gene] = paste0("shared",gene.counter.df$count[gene.counter.df$genes==gene])
    }
  }
  
  #--------following section for DOWN DE genes 
  #plot the average only for the DE genes DOWN + label of shared and unique genes
  DE.label.for.hm.down=NULL
  DE.label.for.hm.down$gene_name = gene.list.down.shared
  DE.label.for.hm.down = data.frame(DE.label.for.hm.down)
  DE.label.for.hm.down$from = "shared"
  #count how often a gene was found DOWN and add to shared 
  gene.counter.down = NULL
  for (gene in stim.response.sig$gene_name[stim.response.sig$avg_logFC<0]) {
    gene.counter.down[[gene]] = sum(stim.response.sig$gene_name[stim.response.sig$avg_logFC<0] == gene)
  }
  gene.counter.down.df = data.frame("genes"= names(gene.counter.down),"count"=unlist(gene.counter.down, use.names=FALSE))
  gene.counter.down.df = arrange(gene.counter.down.df,-count)
  gene.counter.down.df = gene.counter.down.df[!(gene.counter.down.df$genes %in% c(gene.counter.down.df$genes %>%  str_subset(pattern = "^Rp"))),] #filter out ribosomal genes
  for (gene in  DE.label.for.hm.down$gene_name) {
    DE.label.for.hm.down$from[DE.label.for.hm.down$gene_name==gene] = paste0("shared",gene.counter.down.df$count[gene.counter.down.df$genes==gene])
  }
  DE.label.for.hm.down<-DE.label.for.hm.down[order(extract_numeric(DE.label.for.hm.down$from ),decreasing = T),]
  #get unique down genes
  DE.uni.df.down=NULL
  DE.uni.df.down$gene_name = genes.de.down.unique$gene_name
  DE.uni.df.down = data.frame(DE.uni.df.down)
  DE.uni.df.down$from = paste0("cluster",genes.de.down.unique$cluster)
  DE.uni.df.down<-DE.uni.df.down[order(extract_numeric(DE.uni.df.down$from )),]
  DE.heatmap.down = rbind(DE.label.for.hm.down,DE.uni.df.down) # put shared and unique together

###### -------adjust the ordering of labels from most shared to unique UP per cluster
  DE.heatmap = arrange(DE.heatmap,desc(from))
  label.bottom = DE.heatmap$from[DE.heatmap$from %in% c(paste0("shared",2:length(levels(factor(stim.response.sig$cluster)))))]
  DE.heatmap = arrange(DE.heatmap,from)
  label.top = DE.heatmap$from[(DE.heatmap$from %in% c(paste0("cluster",levels(factor(stim.response.sig$cluster)))))]
  DE.heatmap$from = factor(DE.heatmap$from,levels = unique(c(label.bottom,label.top)))
  DE.heatmap = arrange(DE.heatmap,from)
  
#####-----heatmap of UP with labeling if shared and how often
  DoHeatmap(avgexp,features = DE.heatmap$gene_name,draw.lines = F) +
    scale_y_discrete(labels=rev(DE.heatmap$from))
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"DE genes UP average_per_cl_cond_labeled.jpeg"), width=10, height = 20)
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/svg/",sample,"DE genes UP average_per_cl_cond_labeled.svg"), width=10, height = 20)
  #save smaller
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"DE genes UP average_per_cl_cond_labeled_small.jpeg"), width=10, height = 10)
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/svg/",sample,"DE genes UP average_per_cl_cond_labeled_small.svg"), width=10, height = 10)
  #heatmap with labeling if shared and how often + the gene names UP
  DoHeatmap(avgexp,features = DE.heatmap$gene_name,draw.lines = F) +
    scale_y_discrete(labels=rev(paste0(DE.heatmap$from,"_",DE.heatmap$gene_name)))
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"DE genes UP average_per_cl_cond_labeled+gene.jpeg"), width=10, height = 20)
  
####--plot which genes are shared how often UP
  gene.counter.df$genes = factor(gene.counter.df$genes, levels = DE.heatmap$gene_name)
  ggplot(gene.counter.df, aes(x=genes, y=count, fill=count)) +
    geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"UP DE gene_count.jpeg"), width=length(gene.counter.df$genes)/5 , height = 5,limitsize = FALSE)
  
####--plot which genes are shared how often DOWN
  gene.counter.down.df$genes = factor(gene.counter.down.df$genes, levels = DE.heatmap.down$gene_name)
  ggplot(gene.counter.down.df, aes(x=genes, y=count, fill=count)) +
    geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + scale_fill_gradient(low = "black",high = "red")
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"DOWN DE gene_count.jpeg"), width=length(gene.counter.df$genes)/5 , height = 5,limitsize = FALSE)
  
  #-------------------venn diagram of shared + unique DE genes between the clusters-------------------------------------
  genes.per.cluster=NULL
  for (cluster in levels(factor(stim.response.sig$cluster))) {
    genes.per.cluster[[paste0("cluster",cluster)]]=stim.response.sig$gene_name[stim.response.sig$avg_logFC>0 & stim.response.sig$cluster==cluster]
  }
  
  venn.genes.per.cl <- Venn(genes.per.cluster,SetNames = paste0("cluster", levels(factor(stim.response.sig$cluster))))
  jpeg(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"_Venn_digramm_DE genes_UP.jpeg"), width=1600 , height =1200,quality = 100,res = 200)
  plot(venn.genes.per.cl, type = "battle")
  dev.off()
  
  #######################################################################################################################################
  ############------DE testing for overall changes between control and treatment independant of clusters------------#####################
  dir.create("gprofiler2_out_put/overall_DE_test_between_conditions")
  Idents(subset.cluster)="stim"
  stim.response.sig.global <- FindMarkers(subset.cluster,assay = "RNA",test.use = test,min.pct=0.25,logfc.threshold = 0.3 ,ident.1 = treatment_sample, ident.2 = control_sample, verbose = FALSE)
  stim.response.sig.global$gene_name=rownames(stim.response.sig.global)
  
  #---plot heatmap of top 30 DE genes
  stim.response.sig.global$p_val_adj[stim.response.sig.global$p_val_adj == 0] = 1.000000e-300 #genes with significe of 0 are set to 1.000000e-300
  stim.response.sig.global <- arrange(stim.response.sig.global,-avg_logFC)
  top30.up.and.down =  unique(c(stim.response.sig.global$gene_name[1:30],stim.response.sig.global$gene_name[(length(stim.response.sig.global$gene_name)-30):length(stim.response.sig.global$gene_name)]))
  avgexp = AverageExpression(subset.cluster, return.seurat = T)
  DoHeatmap(avgexp,features = top30.up.and.down,draw.lines = F,angle = 0) +ggtitle(paste0("number of total DE genes: ",length(stim.response.sig.global$gene_name))) 
  ggsave(filename = paste0("gprofiler2_out_put/overall_DE_test_between_conditions/",sample,"top_low30 DE genes.jpeg"), width=10, height = 10)
  
  #volcanoplot
  ggplot(stim.response.sig.global) +
    geom_point(aes(x=avg_logFC, y=-log10(p_val_adj)),size=2.5) +
    ggtitle(paste0(" overexpression in ",treatment_sample)) +
    xlab("avg_logFC") + 
    ylab("-log10 adjusted p-value") +
    #scale_y_continuous(limits = c(0,50)) +
    geom_hline(yintercept=(mean(-log10(stim.response.sig.global$p_val_adj))+1), linetype="dashed", color = "grey")+
    geom_vline(xintercept=0.5, linetype="dashed", color = "grey")+
    geom_vline(xintercept=-0.5, linetype="dashed", color = "grey")+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = gene_name),size = 5) +
    theme(legend.title = element_blank(),legend.position = c(0.95,0.1),plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))  
  ggsave(filename = paste0("gprofiler2_out_put/overall_DE_test_between_conditions/",sample,"top_low30 DE genes_volcano.jpeg"), width=8, height =8)
  
  ###################---------plot number of DE genes per cluster Up/Down--------------#####################################################
  ###########################################################################################################################################
  number.of.DE.up = NULL
  number.of.DE.down = NULL
  number.of.DE.df=NULL
  for (cluster.nr in levels(factor(stim.response$cluster))) {
    number.of.DE.up=length(stim.response$gene_name[stim.response$cluster==cluster.nr & stim.response$avg_logFC>0 & stim.response$p_val_adj<0.01])
    number.of.DE.down=length(stim.response$gene_name[stim.response$cluster==cluster.nr & stim.response$avg_logFC<0 & stim.response$p_val_adj<0.01])
    number.of.DE.df.new= data.frame("cluster"=cluster.nr, "DE.up"=number.of.DE.up,"DE.down"=number.of.DE.down)
    number.of.DE.df=rbind(number.of.DE.df.new,number.of.DE.df)
  }
  number.of.DE.df$cluster=factor(number.of.DE.df$cluster,levels = (0:(length(levels(factor(stim.response$cluster)))-1)))
  number.of.DE.df<-arrange(number.of.DE.df,cluster)
  
  ggplot(number.of.DE.df, aes(x = cluster)) + 
    geom_bar(aes(y = DE.up,fill="up"), stat = "identity") + 
    geom_bar(aes(y = -DE.down,fill="down"), stat = "identity") +
    scale_fill_manual(values =c("#3086e1","#e14f30")) +
    scale_y_continuous(labels = c(rev(seq(0,max(number.of.DE.df$DE.down)+5,5)),seq(5,max(number.of.DE.df$DE.up)+5,5)),
                       breaks = c(rev(seq(0,max(number.of.DE.df$DE.down)+5,5))*-1,seq(5,max(number.of.DE.df$DE.up)+5,5)))+
    theme(legend.title = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    xlab("Cluster") +ylab("number of DE genes")
  
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/",sample,"number of Up&Down DE genes per cluster.jpeg"), width=5, height = 5)
  ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/svg/",sample,"number of Up&Down DE genes per cluster.svg"), width=5, height = 5)
  
  ###################---------plot all shared and unique UP DE genes individually----------#####################################################
  ###########################################################################################################################################
  #vln plots for unique genes per cluster
  dir.create("gprofiler2_out_put/volc&heatmap_plots_of_DE/individual_DE_up_genes_plots")
  Idents(subset.cluster)="seurat_clusters"
  for (cluster.nr in unique(genes.de.up.unique$cluster)) {
    if (length(genes.de.up.unique$gene_name[genes.de.up.unique$cluster==cluster.nr])>1) {
      print(VlnPlot(subset.cluster,features = genes.de.up.unique$gene_name[genes.de.up.unique$cluster==cluster.nr],stack = T,flip = T,split.by = "stim",pt.size = 0,cols = c('#00C5C0','#FA756C')) +
        geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
        ggtitle(label = paste0("unique DE genes cluster ",cluster.nr))+
        theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0)) + NoLegend())
      ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/individual_DE_up_genes_plots/",sample,"_",cluster.nr,"_vln.jpeg"), width=length(unique(genes.de.up.unique$cluster)) , height = length(unique(genes.de.up.unique$gene_name[genes.de.up.unique$cluster==cluster.nr]))+1)
    } else { gene = genes.de.up.unique$gene_name[genes.de.up.unique$cluster==cluster.nr]
      print(VlnPlot(subset.cluster,features = gene,split.by = "stim",pt.size = 0,cols = c('#00C5C0','#FA756C')) +
              geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
              ggtitle(label = paste0("cluster ",cluster.nr," one unique gene: ",gene))+
              theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0)) + NoLegend())
      ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/individual_DE_up_genes_plots/",sample,"_",cluster.nr,"_vln.jpeg"), width=length(unique(genes.de.up.unique$cluster)) , height = 2)
    }
   }

  #vln plots for shared genes grouped by how often shared
  dir.create(paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/shared_DE_up_genes_plots")) 
  Idents(subset.cluster)="seurat_clusters"
  shared.df = gene.counter.df[gene.counter.df$count>1,]
  for (shared.nr in unique(shared.df$count)) {
    if (length(shared.df$genes[shared.df$count==shared.nr])>1) {
      print(VlnPlot(subset.cluster,features = shared.df$genes[shared.df$count==shared.nr],stack = T,flip = T,split.by = "stim",pt.size = 0,cols = c('#00C5C0','#FA756C')) +
              geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
              ggtitle(label = paste0("shared DE genes between ",shared.nr," clusters"))+
              theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.5)) + NoLegend())
      ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/shared_DE_up_genes_plots/",sample,"_shared_DE_genes_between_",shared.nr,"_clusters_vln.jpeg"), width=length(unique(genes.de.up.unique$cluster)) , height = length(shared.df$genes[shared.df$count==shared.nr])+1)
    } else {gene = shared.df$genes[shared.df$count==shared.nr]
    print(VlnPlot(subset.cluster,features = gene,split.by = "stim",pt.size = 0,cols = c('#00C5C0','#FA756C')) +
            geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) +
            ggtitle(label = paste0("shared DE genes between ",shared.nr," one shared gene: ",gene))+
            theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.5)) + NoLegend())
    ggsave(filename = paste0("gprofiler2_out_put/volc&heatmap_plots_of_DE/shared_DE_up_genes_plots/",sample,"_shared_DE_gene_between_",shared.nr,"_clusters_vln.jpeg"), width=length(unique(genes.de.up.unique$cluster)) , height = 2)
    }
  }
  
 
  