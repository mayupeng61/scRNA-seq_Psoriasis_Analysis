# scRNA-seq_Psoriasis_Analysis
# Analysis Code for [Celastrol disrupts the inflammatory feedback loop driven by Tenascin-C⁺ fibroblast subtype in psoriasis: Insights from single-cell RNA sequencing and experimental validation]

This repository contains the source code used for the single-cell RNA-seq analysis presented in our manuscript.

## Repository Structure

The analysis pipeline is organized into numbered scripts to ensure reproducibility. Please execute them in the following order:

### R Analysis Pipeline (Seurat & Downstream)
* `01_load_qc_clustering_markers.R`: Initial data loading, quality control, normalization, and clustering.
* `02_celltype_annotation_and_visualization.R`: Cell type annotation and generation of UMAP plots.
* `03_CellChat_Multisample_Analysis.R`: Cell-cell communication analysis using CellChat.
* `04_fibroblasts_subclustering.R`: Sub-clustering analysis specifically for the Fibroblast population.
* `05_fibroblasts_expression_and_mechanism_analysis.R`: Differential expression and mechanistic studies within fibroblast subtypes.
* `06_pathway_enrichment_analysis.R`: GO/KEGG and other pathway enrichment analyses (using clusterProfiler, irGSEA, etc.).
* `07_Regulon_AUCell_RSS_Fibroblasts.R`: Post-processing of SCENIC results, including AUCell scoring and RSS visualization.

### Python Analysis Pipeline (GRN Inference)
* `run_scenic.py`: The main script for running pySCENIC (GRNBoost2 -> cisTarget -> AUCell).
    * *Note: This script requires high-performance computing resources due to memory usage.*

## Environment & Dependencies

### R Environment
The full list of R packages and versions used in this analysis is detailed in `sessionInfo.txt`. 
Key packages include:
* Seurat v5.3.0
* CellChat v2.1.2
* clusterProfiler
* irGSEA

### Python Environment
The specific versions for the Python environment (pySCENIC) are listed in `requirements.txt`.
To install dependencies:
```bash
pip install -r requirements.txt
