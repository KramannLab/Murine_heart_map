library(Seurat)
library(ggplot2)
setwd("~/Documents/FP_scRNA/R stuff/PVM RDS/")

# ------------------------------------------------------------------------------------------- #
# basic filter setting for each sample of the perivascular map (nfeature + mt-procentage cutoffs )
# filtering out cells with Hba-a1, Hba-a2 and Hbb-bs contamination, as well as immune cells by Ptprc
# followed by integration of Sham and Tac per genotype (tdtomato excluded form anchor set)
# ------------------------------------------------------------------------------------------- #

#load for -------------------------------------------------------PDGFRb samples--------------------------------------------------
sample="PDGFRb_Sham&TAC_integrated"
ctrl_sample = 'FP3_PDGFRb_Sham'
stim_sample = 'FP14_PDGFRb_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP3_PDGFRb_Sham/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP14_PDGFRb_TAC/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')
VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
ctrl.pdgfrb = ctrl

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')
VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
stim.pdgfrb = stim

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Gli1 samples--------------------------------------------------
sample="Gli1_Sham&TAC_integrated"
ctrl_sample = 'FP13_Gli1_Sham'
stim_sample = 'FP12_Gli1_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP13_Gli1_Sham/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP12_Gli1_TAC/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')
VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
ctrl.gli1 = ctrl

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')
VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 6)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
stim.gli1 = stim

#######-Gli1 samples--28day time point-##############Set up stimulated object - standard filter + normalization
sample="FP71_Gli1_TAC_28d"

gli1.28d.tac.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/21Oct37-DL012-FP71/outs/filtered_feature_bc_matrix/")
# - standard filter + normalization
gli1.28d.tac <- CreateSeuratObject(counts = gli1.28d.tac.data, project = sample, min.cells = 3, min.features = 200)
gli1.28d.tac$stim <- "FP71_Gli1_TAC_28d"
gli1.28d.tac[['percent.mt']] = PercentageFeatureSet(gli1.28d.tac, pattern = '^mt-')
VlnPlot(gli1.28d.tac, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
gli1.28d.tac <- subset(gli1.28d.tac, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 6)
VlnPlot(gli1.28d.tac, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)
#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(gli1.28d.tac,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
gli1.28d.tac=subset(gli1.28d.tac,cells = WhichCells(gli1.28d.tac,expression = `Hba-a1` ==0))
gli1.28d.tac=subset(gli1.28d.tac,cells = WhichCells(gli1.28d.tac,expression = `Hba-a2` ==0))
gli1.28d.tac=subset(gli1.28d.tac,cells = WhichCells(gli1.28d.tac,expression = `Hbb-bs` ==0))
gli1.28d.tac=subset(gli1.28d.tac,cells = WhichCells(gli1.28d.tac,expression = Ptprc ==0))
#---normalize and findvariable features-----------------------
gli1.28d.tac <- NormalizeData(gli1.28d.tac, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
gli1.28d.tac <- FindVariableFeatures(gli1.28d.tac, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(gli1.28d.tac, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Col1a1 samples--------------------------------------------------
sample="Col1a1_Sham&TAC_integrated"
ctrl_sample = 'FP19_Col1a1_Sham'
stim_sample = 'FP18_Col1a1_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP19_Col1a1_Sham/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP18_Col1a1_TAC/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')
VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 6)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
ctrl.col1a1 = ctrl

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')
VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 6)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
stim.col1a1= stim

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Cdh5 samples--14day time point --------------------------------
sample="Cdh5_Sham&TAC_integrated"
ctrl_sample = 'FP21_Cdh5_Sham'
stim_sample = 'FP20_Cdh5_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP21_Cdh5_Sham/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP20_Cdh5_TAC/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')
VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
ctrl.cdh5 = ctrl

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')
VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
stim.cdh5 = stim

#######-Cdh5 samples--28day time point-##############Set up stimulated object - standard filter + normalization
sample="FP72_Cdh5_TAC_28d"

chd5.28d.tac.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/21Oct37-DL013-FP72/outs/filtered_feature_bc_matrix/")
# - standard filter + normalization
chd5.28d.tac <- CreateSeuratObject(counts = chd5.28d.tac.data, project = sample, min.cells = 3, min.features = 200)
chd5.28d.tac$stim <- "FP72_Cdh5_TAC_28d"
chd5.28d.tac[['percent.mt']] = PercentageFeatureSet(chd5.28d.tac, pattern = '^mt-')
VlnPlot(chd5.28d.tac, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
chd5.28d.tac <- subset(chd5.28d.tac, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 6)
VlnPlot(chd5.28d.tac, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(sample,"_Gene&RNA_count_filtered.jpeg"), width=10 , height = 10)
#---filter out contaminations by erythrocytes and immune cells (cd45)-----------------
VlnPlot(chd5.28d.tac,features = c("Hba-a1","Hba-a2","Hbb-bs","Ptprc"),pt.size = 0.1,ncol = 4)
ggsave(filename = paste0(sample,"_hemoglobin&immune_contamination.jpeg"), width=10 , height = 10)
chd5.28d.tac=subset(chd5.28d.tac,cells = WhichCells(chd5.28d.tac,expression = `Hba-a1` ==0))
chd5.28d.tac=subset(chd5.28d.tac,cells = WhichCells(chd5.28d.tac,expression = `Hba-a2` ==0))
chd5.28d.tac=subset(chd5.28d.tac,cells = WhichCells(chd5.28d.tac,expression = `Hbb-bs` ==0))
chd5.28d.tac=subset(chd5.28d.tac,cells = WhichCells(chd5.28d.tac,expression = Ptprc ==0))
#---normalize and findvariable features-----------------------
chd5.28d.tac <- NormalizeData(chd5.28d.tac, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
chd5.28d.tac <- FindVariableFeatures(chd5.28d.tac, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#------save individual filtered objects as RDS
saveRDS(chd5.28d.tac, file = paste0(sample, "_seurat_object_after_bacis_filters.rds"))

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------NG2 samples--------------------------------------------------
sample="NG2_Sham&TAC_integrated"
ctrl_sample = 'FP23_NG2_Sham'
stim_sample = 'FP22_NG2_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP23_NG2_Sham/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP22_NG2_TAC/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')
VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
ctrl.ng2 = ctrl

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')
VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
stim <- subset(stim, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 6)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
stim.ng2 = stim

#-------------------------------------------------------------------------------------------------------------------------------
#load for -------------------------------------------------------Myh11 samples--------------------------------------------------
sample="Myh11_Sham&TAC_integrated"
ctrl_sample = 'FP56_Myh11_Sham'
stim_sample = 'FP57_Myh11_TAC'
ctrl.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP56_Myh11_Sham/filtered_feature_bc_matrix")
stim.data <- Read10X(data.dir = "/media/kramannworkstation/HD3/HD3/PVM_cellranger_outs/FP57_Myh11_TAC/filtered_feature_bc_matrix")

#######-CTRL-##############Set up ctrl object - standard filter + normalization
ctrl <- CreateSeuratObject(counts = ctrl.data, project = ctrl_sample, min.cells = 3, min.features = 200)
ctrl$stim <- ctrl_sample
ctrl[['percent.mt']] = PercentageFeatureSet(ctrl, pattern = '^mt-')
VlnPlot(ctrl, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(ctrl_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
ctrl.myh11 = ctrl

#######-STIM-##############Set up stimulated object - standard filter + normalization
stim <- CreateSeuratObject(counts = stim.data, project = stim_sample, min.cells = 3, min.features = 200)
stim$stim <- stim_sample
stim[['percent.mt']] = PercentageFeatureSet(stim, pattern = '^mt-')
VlnPlot(stim, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.5)
ggsave(filename = paste0(stim_sample,"_Gene&RNA_count_unfiltered.jpeg"), width=10 , height = 10)
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
#-----transfer to individual variable for later integration with Gli1 and Col1a1 (Myh11)
stim.myh11 = stim

############################################################################################################################################
##################################################-----integrations--------------------#####################################################
############################################################################################################################################
#all datasets containing tdTomato tranced Endothelial cells
#perform integration of Cdh5 Sham, TAC 14 days and TAC 28 days
sample="Endothelial_datasets"
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl.cdh5, stim.cdh5 ,chd5.28d.tac), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_integrated_after_bacis_filters+TAC28d.rds"))

#all datasets containing tdTomato tranced Fibroblasts cells
#perform integration of Pdgfrb, Col1a1, Gli1 Sham+TAC each (+28d gli1 tac)
sample="Fibroblast_datasets"
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl.pdgfrb, stim.pdgfrb, ctrl.col1a1, stim.col1a1, ctrl.gli1, stim.gli1, gli1.28d.tac), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_integrated_after_bacis_filters+TAC28d.rds"))

#all datasets containing tdTomato tranced mural cells (VSMC+Pericytes) + neuronal cells
#futher analysis showed that for the mural cells a seperate prefiltering before integration is necessary
ctrl.pdgfrb <- readRDS("~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering/FP3_PDGFRb_Sham_seurat_object_after_bacis_filters.rds")
stim.pdgfrb <- readRDS("~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering/FP14_PDGFRb_TAC_seurat_object_after_bacis_filters.rds")
ctrl.ng2 <- readRDS("~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering/FP23_NG2_Sham_seurat_object_after_bacis_filters.rds")
stim.ng2 <- readRDS("~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering/FP22_NG2_TAC_seurat_object_after_bacis_filters.rds")
ctrl.myh11 <- readRDS("~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering/FP56_Myh11_Sham_seurat_object_after_bacis_filters.rds")
stim.myh11 <- readRDS("~/Documents/FP_scRNA/R stuff/PVM RDS/RDS_basic_filtering/FP57_Myh11_TAC_seurat_object_after_bacis_filters.rds")
source("~/Documents/FP_scRNA/R stuff/scripts/new PVM scripts/usefull.fun.sc.PVM.R")
ctrl.pdgfrb <- recluster_RNA(ctrl.pdgfrb,res = 0.1)
FeaturePlot(ctrl.pdgfrb,features = c("nFeature_RNA","tdTomato-Green","Pecam1","Pdgfra"))
ctrl.pdgfrb <- subset(ctrl.pdgfrb,idents=2) #only keep mural cells
stim.pdgfrb <- recluster_RNA(stim.pdgfrb,res = 0.1)
FeaturePlot(stim.pdgfrb,features = c("nFeature_RNA","tdTomato-Green","Pecam1","Pdgfra"))
stim.pdgfrb <- subset(stim.pdgfrb,idents=2) #only keep mural cells
ctrl.ng2 <- recluster_RNA(ctrl.ng2,res = 0.2)
FeaturePlot(ctrl.ng2,features = c("nFeature_RNA","tdTomato-Green","Pecam1","Pdgfra"))
ctrl.ng2 <- subset(ctrl.ng2,idents=c(1,2)) #only keep mural cells
stim.ng2 <- recluster_RNA(stim.ng2,res = 0.2)
FeaturePlot(stim.ng2,features = c("nFeature_RNA","tdTomato-Green","Pecam1","Pdgfra"))
stim.ng2 <- subset(stim.ng2,idents=c(1,2)) #only keep mural cells
ctrl.myh11 <- recluster_RNA(ctrl.myh11,res = 0.1)
FeaturePlot(ctrl.myh11,features = c("nFeature_RNA","tdTomato-Green","Pecam1","Pdgfra"))
ctrl.myh11 <- subset(ctrl.myh11,idents=1) #only keep mural cells
stim.myh11 <- recluster_RNA(stim.myh11,res = 0.1)
FeaturePlot(stim.myh11,features = c("nFeature_RNA","tdTomato-Green","Pecam1","Pdgfra"))
stim.myh11 <- subset(stim.myh11,idents=1) #only keep mural cells
#perform integration of Pdgfrb, Ng2, Myh11 Sham+TAC each
sample="Mural_datasets"
Samples.anchors <- FindIntegrationAnchors(object.list = list(ctrl.myh11, stim.myh11,ctrl.ng2, stim.ng2,ctrl.pdgfrb, stim.pdgfrb), dims = 1:20)
#exclude tdtomato from integration features
anchor.list = Samples.anchors@anchor.features
anchor.list = anchor.list[!anchor.list %in% "tdTomato-Green"]
Samples.anchors@anchor.features=anchor.list
Samples.combined <- IntegrateData(anchorset = Samples.anchors, dims = 1:20)
saveRDS(Samples.combined, file = paste0(sample, "_integrated_after_bacis_filters.rds"))



