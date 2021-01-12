library(Seurat)
library(ggplot2)
setwd("~/Documents/FP_scRNA/R stuff")

# ------------------------------------------------------------------------------------------- #
# basic filter setting for each sample of the perivascular map (nfeature + mt-procentage cutoffs )
# filtering out cells with Hba-a1, Hba-a2 and Hbb-bs contamination
# integration of Sham and Tac per genotype (tdtomato excluded form anchor set)
# ------------------------------------------------------------------------------------------- #

#load for -------------------------------------------------------PDGFRb samples--------------------------------------------------
sample="PDGFRb_Sham&TAC_integrated"
ctrl_sample = 'FP3_PDGFRb_Sham'
stim_sample = 'FP14_PDGFRb_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP3_PDGFRb_Sham/outs/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP14_PDGFRb_TAC/outs/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(ctrl,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(ctrl_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a1` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a2` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hbb-bs` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
ctrl <- NormalizeData(ctrl, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(ctrl, file = paste0(ctrl_sample, "_seurat_object_after_bacis_filters.rds"))

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(stim,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(stim_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a1` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a2` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hbb-bs` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
stim <- NormalizeData(stim, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(stim, file = paste0(stim_sample, "_seurat_object_after_bacis_filters.rds"))

#------------perform pairwise integration
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Gli1 samples--------------------------------------------------
sample="Gli1_Sham&TAC_integrated"
ctrl_sample = 'FP13_Gli1_Sham'
stim_sample = 'FP12_Gli1_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP13_Gli1_Sham/outs/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP12_Gli1_TAC/outs/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(ctrl,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(ctrl_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a1` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a2` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hbb-bs` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
ctrl <- NormalizeData(ctrl, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(ctrl, file = paste0(ctrl_sample, "_seurat_object_after_bacis_filters.rds"))

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(stim,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(stim_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a1` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a2` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hbb-bs` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
stim <- NormalizeData(stim, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(stim, file = paste0(stim_sample, "_seurat_object_after_bacis_filters.rds"))

#------------perform pairwise integration
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Col1a1 samples--------------------------------------------------
sample="Col1a1_Sham&TAC_integrated"
ctrl_sample = 'FP19_Col1a1_Sham'
stim_sample = 'FP18_Col1a1_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP19_Col1a1_Sham/outs/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP18_Col1a1_TAC/outs/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(ctrl,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(ctrl_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a1` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a2` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hbb-bs` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
ctrl <- NormalizeData(ctrl, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(ctrl, file = paste0(ctrl_sample, "_seurat_object_after_bacis_filters.rds"))

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(stim,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(stim_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a1` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a2` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hbb-bs` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
stim <- NormalizeData(stim, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(stim, file = paste0(stim_sample, "_seurat_object_after_bacis_filters.rds"))

#------------perform pairwise integration
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Cdh5 samples--------------------------------------------------
sample="Cdh5_Sham&TAC_integrated"
ctrl_sample = 'FP21_Cdh5_Sham'
stim_sample = 'FP20_Cdh5_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP21_Cdh5_Sham/outs/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP20_Cdh5_TAC/outs/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(ctrl,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(ctrl_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a1` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a2` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hbb-bs` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
ctrl <- NormalizeData(ctrl, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(ctrl, file = paste0(ctrl_sample, "_seurat_object_after_bacis_filters.rds"))

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(stim,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(stim_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a1` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a2` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hbb-bs` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
stim <- NormalizeData(stim, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(stim, file = paste0(stim_sample, "_seurat_object_after_bacis_filters.rds"))

#------------perform pairwise integration
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------NG2 samples--------------------------------------------------
sample="NG2_Sham&TAC_integrated"
ctrl_sample = 'FP23_NG2_Sham'
stim_sample = 'FP22_NG2_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP23_NG2_Sham/outs/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP22_NG2_TAC/outs/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(ctrl,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(ctrl_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a1` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a2` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hbb-bs` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
ctrl <- NormalizeData(ctrl, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(ctrl, file = paste0(ctrl_sample, "_seurat_object_after_bacis_filters.rds"))

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(stim,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(stim_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a1` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a2` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hbb-bs` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
stim <- NormalizeData(stim, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(stim, file = paste0(stim_sample, "_seurat_object_after_bacis_filters.rds"))

#------------perform pairwise integration
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Myh11 samples--------------------------------------------------
sample="Myh11_Sham&TAC_integrated"
ctrl_sample = 'FP56_Myh11_Sham'
stim_sample = 'FP57_Myh11_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP56/outs/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/Seagate Expansion Drive/cellranger outs/FP57/outs/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(ctrl,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(ctrl_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a1` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hba-a2` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = `Hbb-bs` ==0))
ctrl=subset(ctrl,cells = WhichCells(ctrl,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
ctrl <- NormalizeData(ctrl, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(ctrl, file = paste0(ctrl_sample, "_seurat_object_after_bacis_filters.rds"))

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)

#filter low feature and high mt cells
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 6)

VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)

#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(stim,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(stim_sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a1` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hba-a2` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = `Hbb-bs` ==0))
stim=subset(stim,cells = WhichCells(stim,expression = Ptprc ==0))

#---normalize and findvariable features-----------------------
stim <- NormalizeData(stim, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(stim, file = paste0(stim_sample, "_seurat_object_after_bacis_filters.rds"))

#------------perform pairwise integration
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))
