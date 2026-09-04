# IS_LRG
# Code for: Immunological Signatures of Lactylation-related Genes in Ischemic Stroke: a Bulk and Single-cell Transcriptomic Study

This repository contains the custom R scripts and analysis pipelines used to generate the results presented in the manuscript: **"Immunological Signatures of Lactylation-related Genes in Ischemic Stroke: a Bulk and Single-cell Transcriptomic Study"**.

## 1. Overview

This study integrates microarray expression data and single-cell RNA-sequencing (scRNA-seq) data to investigate the role of lactylation-related genes in ischemic stroke. The analysis pipeline includes data preprocessing, batch effect correction, differential expression analysis, immune infiltration estimation, single-cell clustering, and trajectory inference, etc.

## 2. System Requirements and Dependencies

All analyses were performed using **R (version 4.2.2)**.

The following major R packages are required to run the scripts:

- **Bulk RNA-seq Analysis:** limma, sva, preprocessCore, pheatmap, dplyr, ggplot2, ggrepel, clusterProfiler, glmnet, randomForest, kernlab, pROC, GSVA, ConsensusClusterPlus, etc.
- **scRNA-seq Analysis:** Seurat, DoubletFinder, Harmony, Monocle 2, etc.

## 3. Data Availability

The raw datasets analyzed in this study are publicly available from the Gene Expression Omnibus (GEO) database. To run the code, please download the following datasets and place them in the `data/` directory:

- **Bulk RNA-seq cohorts:** GSE16561, GSE22255, GSE37587
- **scRNA-seq cohort:** GSE225948
