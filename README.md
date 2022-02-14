# Mapping the cardiac vascular niche in heart failure
## Summary 
Will be provided after publication.

## Data availabitly
Processed and raw single cell RNA sequencing data will be availablea after publication via the Gene Expression Omnibus under the series number GSE166403.

## Software
Raw sequencing files were processed with cell ranger (v3.0.2). All scripts provided were run with R (version 4.0.2). Most processing and analysis was performed with Seurat (v3.2.2). For futher information on packages used, see references/R_sessionInfo.txt . 

## Script discription
#### script for first filtering of PVM datasets
This script was used to load the datasets after cell ranger alignment in order to perform first quality filtering before any subsequent analysis.
Samples were filtered for nfeature, mt-reads, immune cell and erythocyte contaminations.
Integration of fibroblast, endothelial cell and mural cell data sets was performed afterwards by seurat cca.

#### 14 sample integration based on harmony
This script was used for the integration of all 14 data sets based on harmony and subsequent analysis. The script for first filtering of PVM datasets was run beforehand to provide the input objects.

#### filterLowTdtomCl_function_harmony
This script provides the function that was used to filter cells based on tdTomato reads and is called in the harmony integration script.

#### script Endothelial cells integrated analysis
This script was used to analyse the integrated dataset for endothelial cells from Cdh5 samples (sham and TAC each) based on seurat. The script for first filtering of PVM datasets was run beforehand to provide the input objects.

#### script Fibroblast cells integrated analysis
This script was used to analyse the integrated dataset for fibroblasts cells from Col1a1 Gli1 Pdgfrb samples (sham and TAC each) based on seurat. The script for first filtering of PVM datasets was run beforehand to provide the input objects.

#### script for Mural cell integrated analysis
This script was used to analyse the integrated dataset for mural cells from Myh11, Ng2 and Pdgfrb samples (sham and TAC each) based on seurat. The script for first filtering of PVM datasets was run beforehand to provide the input objects.

#### helper_functions_DEG_GO_ploting
This script provides several function to perfrom differential gene expression and gene set enrichment analysis, as well as some small function for plotting specific data.

#### pathway activity prediction_progeny
This script was used to estimate signal pathway activity per cell.

#### pathway activity prediction_dorothea
This script was used to estimate transcription factor activity per cell.

#### ECMscore correlation
This script was used to perform a correlation analysis of the ECMscore to gene expression in the fibroblast dataset.

#### python-correlation
This script provides a helper function for the correlation analysis in python.

#### script for Forte MI data








