**Mast cell stabilization inhibits benign prostatic hyperplasia — analysis code**

This repository contains the code used to process single-cell RNA-seq data and generate the main and supplementary figures for our study on mast cell–driven stromal remodeling in benign prostatic hyperplasia (BPH), including trajectory inference (Monocle3) and cell–cell communication analysis (CellChat).

Note: This repo provides analysis code only. Raw data are not included.

Data processing pipeline
	•	raw data process/Quality Control.R
Quality control, filtering, and basic QC visualizations.
	•	raw data process/Integration.R
Integration and dimensionality reduction (Harmony-based workflow) and UMAP generation.
	•	raw data process/Major Annotation.R
Major lineage annotation (broad cell types).
	•	raw data process/Minor Annotation/
Subclustering and minor annotations for specific lineages.


**Requirements**
	•	R (>= 4.2 recommended)
	•	Key R packages:
	•	Seurat
	•	SeuratWrappers
	•	harmony
	•	monocle3
	•	CellChat
	•	ggplot2, dplyr, tidyr, patchwork, Matrix

**Getting started**

1) Install packages (example)

Install base packages:
	•	ggplot2, dplyr, tidyr, patchwork, Matrix

Then install:
	•	Seurat
	•	harmony
	•	monocle3 (follow official installation for your OS)
	•	CellChat (follow official installation)
	•	SeuratWrappers (for conversions if needed)

2) Configure paths

These scripts assume you have local access to:
	•	a gene expression matrix / Seurat objects produced after basic preprocessing
	•	an output directory for intermediate .rds objects and figures

Recommended: create a config.R in the repo root:
	•	DATA_DIR: path to input data
	•	OUT_DIR: path to outputs

Then add near the top of each script:
	•	source(“config.R”)

Suggested run order
	1.	Quality control
Run: raw data process/Quality Control.R
Output: QC-filtered Seurat object(s) saved as .rds
	2.	Integration + UMAP
Run: raw data process/Integration.R
Output: integrated Seurat object(s) with Harmony reduction and UMAP embeddings
	3.	Annotation
Run: raw data process/Major Annotation.R
Run: scripts under raw data process/Minor Annotation/
Output: annotated Seurat object(s), marker tables, and annotation metadata
	4.	Downstream analyses (used in figures)
Monocle3 pseudotime/trajectory inference (root selection supported)
CellChat ligand–receptor communication inference
Gene/module score visualization and feature plots
	5.	Reproduce figures
Run: Figure 1.R, Figure 2.R, Figure 3.R, Figure 4.R
Run: Figure S*.R

**Notes on reproducibility**
	•	Set a seed where applicable (e.g., UMAP, clustering)
	•	Record session info (sessionInfo())

**Contact**
For questions or reproducibility issues, please open an issue in this repository.
