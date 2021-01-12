#https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
############-------------------Cell cycle analysis ---------------------------------------------------
  library(Seurat)
  library(stringr)
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- str_to_title(cc.genes$s.genes) #cc.genes is part of seurat
  g2m.genes <- str_to_title(cc.genes$g2m.genes)
  
  
  #perform scoring for each cell
  dir.create("cell_cycle_analysis")
  Samples.combined <- CellCycleScoring(Samples.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  FeaturePlot(Samples.combined,features = "S.Score",pt.size = 1,order = T)
  ggsave(filename = paste0("cell_cycle_analysis/",sample,"_S.score.jpeg"), width=10 , height = 10)
  VlnPlot(Samples.combined,features = "S.Score",group.by = "seurat_clusters",cols = c('#00C5C0','#FA756C'),split.by = "stim",split.plot = T,pt.size = 0) + theme(axis.text.x = element_text(hjust = 0.5,angle = 0))
  ggsave(filename = paste0("cell_cycle_analysis/",sample,"_S.score_vln.jpeg"), width=10 , height = 10)
  
  
  FeaturePlot(Samples.combined,features = "G2M.Score",pt.size = 1,order = T)
  ggsave(filename = paste0("cell_cycle_analysis/",sample,"_G2M.score.jpeg"), width=10 , height = 10)
  VlnPlot(Samples.combined,features = "G2M.Score",group.by = "seurat_clusters",cols = c('#00C5C0','#FA756C'),split.by = "stim",split.plot = T,pt.size = 0) + theme(axis.text.x = element_text(hjust = 0.5,angle = 0))
  ggsave(filename = paste0("cell_cycle_analysis/",sample,"_G2M.score_vln.jpeg"), width=10 , height = 10)
  
  
  DimPlot(Samples.combined,group.by = "Phase",pt.size = 1)
  ggsave(filename = paste0("cell_cycle_analysis/",sample,"_cc_Phase.jpeg"), width=10 , height = 10)
  