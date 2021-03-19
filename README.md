# Mapping the cardiac vascular niche in heart failure
## Summary 

The cardiac vasculature and perivasculature are of major importance in homeostasis and during disease but we lack a complete understanding of its cellular heterogeneity and alteration in response to injury as a major driver of heart failure. Using combined genetic fate tracing with confocal imaging and single cell RNA sequencing of this niche in homeostasis and during heart failure, we unraveled a previously unappreciated heterogeneity of vascular and perivascular cell types including various distinct endothelial, fibroblast, pericyte and vascular smooth muscle cell populations with specific transcriptional changes in response to injury. Interstitial fibrosis was driven by various fibroblast subsets, with Thbs4+ fibroblasts showing the highest extracellular matrix gene expression. We identified that injury causes a loss of cellular cross-talk towards endothelium in the vascular niche as a potential driver of endothelial dysfunction. Our analysis of mural cells identified a previously unappreciated role of pericytes contributing to cardiac fibrosis and remodeling in heart failure.

## Data availabitly
Processed and raw single cell RNA sequencing data will be availablea after publication via the Gene Expression Omnibus under the series number GSE166403.

## Script discription

script for filtering of PVM datasets
This script was used for the initial first loading of the cell ranger out put (filtered feature barcode matrix) and filtering of each individual data set for high quality cells, excluding extreme high and low feature count cells, cells with high mitochondrial percentage, cell with hemoglobin contamination (erythrocytes) and cells with reads for Ptprc (immune cell contamination).

pairwise_integration (script Cdh5 / Col1a1 / Gli1 / Myh11 / NG2 / PDGFRb Sham and TAC pairwise integrated)
These scripts were used for further filtering and analysis of pairwise integrated (seurat) dataset, proving the input for the script for the integration of fibroblast and mural cells. To be noticed, the actually pairwise integration function is already run in the script for filtering of PVM datasets and is prerequisite for these scripts. The Cdh5 script here is more complex as all analysis specific for endothelial cells were run here (required reference here: EC-atlas Heart heatmap table.csv)

integration fibroblasts Col1a1 Gli1 Pdgfrb
This script was used to integrate and process dataset from Col1a1 Gli1 Pdgfrb samples (sham and TAC each) based on seurat. Pairwise integration scripts were run beforehand to provide the input objects.

integration VSMCs and pericytes Pdgfrb Myh11 Ng2
This script was used to integrate and process dataset from Pdgfrb Myh11 Ng2 samples (sham and TAC each) based on seurat. Pairwise integration scripts were run beforehand to provide the input objects.

integration of all samples_harmony
This script was used for the integration of all 12 data sets based on harmony and subsequent analysis. The script for first filtering of PVM datasets was run beforehand to provide the input objects.

filterLowTdtomCl_function_harmony
This script provides the function that was used to filter cells based on tdTomato reads and is called in the harmony integration script.

cell cycle analysis
This script was used to perform cell cycle phase scoring per cell.

dorothea testing
This script was used to estimate transcription factor activity per cell.

ecm_score
This script was used to perform scoring based on gene sets from the matrisome database (http://matrisomeproject.mit.edu/). The gene sets are also provided under references.

gProfileR2 for GO on DE_sig_filter_add
This script was used to calculate differentially expressed genes per cluster and marker genes per cluster. Based on the results, association with GO, KEGG and Reactome was tested.

Progeny testing
This script was used to estimate signal pathway activity per cell.
