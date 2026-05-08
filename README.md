# Mapping the cardiac vascular niche in heart failure
## Summary 
Abstract
The cardiac vascular and perivascular niche are of major importance in homeostasis and during disease, but we lack a complete understanding of its cellular heterogeneity and alteration in response to injury as a major driver of heart failure. Using combined genetic fate tracing with confocal imaging and single-cell RNA sequencing of this niche in homeostasis and during heart failure, we unravel cell type specific transcriptomic changes in fibroblast, endothelial, pericyte and vascular smooth muscle cell subtypes. We characterize a specific fibroblast subpopulation that exists during homeostasis, acquires Thbs4 expression and expands after injury driving cardiac fibrosis, and identify the transcription factor TEAD1 as a regulator of fibroblast activation. Endothelial cells display a proliferative response after injury, which is not sustained in later remodeling, together with transcriptional changes related to hypoxia, angiogenesis, and migration. Collectively, our data provides an extensive resource of transcriptomic changes in the vascular niche in hypertrophic cardiac remodeling.

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

#### transcription factor activity prediction_dorothea
This script was used to estimate transcription factor activity per cell.

#### ECMscore correlation
This script was used to perform a correlation analysis of the ECMscore to gene expression in the fibroblast dataset.

#### python-correlation
This script provides a helper function for the correlation analysis in python.

#### trajectory pseudotime analysis monocle3
This script was used to perform trajectory analysis with monocle3 and identify pseudotime asocitade genes.

#### scVelo
This script is python based and was used to perform RNA velocity analysis with scVelo.
Spliced/unspliced alignment was analysed with command line function from velocyto.

#### script for Forte MI data
This script was used to process the dataset from Forte et al. as discribed in their study.

#### script for McLellan data AngII
This script was used to process the dataset from McLellan et al. as discribed in their study.

#### script for integration PVM forteMI
This script was used to integrated fibroblasts from this study and fibroblast from the Forte et al. data.

#### script for integration PVM McLellanAngII
This script was used to integrated fibroblasts from this study and fibroblast from the McLellan et al. data.


