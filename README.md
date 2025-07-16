# COVID-19 Bulk RNA-seq Analysis with Comprehensive GSEA

A complete bioinformatics pipeline for analysing bulk RNA-seq data from COVID-19 patients across different disease severities, featuring comprehensive Gene Set Enrichment Analysis (GSEA) and pathway interpretation.

## 🔬 Project Overview

This pipeline analyzes RNA-seq count data from COVID-19 patients to:
- Compare gene expression across disease severities (Mild, Severe, Critical vs Normal)
- Identify differentially expressed genes and pathways
- Perform comprehensive GSEA using multiple gene set databases
- Generate publication-ready visualizations and reports

## 📊 Dataset Information

- **Study**: GSE293708 (COVID-19 severity analysis)
- **Sample Groups**: 
  - Normal (H1-H4)
  - COVID-19 Mild (CM1-CM4)
  - COVID-19 Severe (CS1-CS4)
  - COVID-19 Critical (CI1-CI4)
- **Total Samples**: 16
- **Analysis Type**: Bulk RNA-seq

## 🚀 Pipeline Features

### Core Analysis
- ✅ **Data Loading & Preprocessing**: Automatic sample classification and quality filtering
- ✅ **Normalization**: DESeq2-based size factor normalization
- ✅ **Gene Annotation**: Automatic gene ID to symbol conversion
- ✅ **Quality Control**: PCA, sample correlation, clustering analysis
- ✅ **Differential Expression**: Multi-group comparison using DESeq2

### Advanced GSEA Analysis
- ✅ **Multiple Gene Set Databases**:
  - Hallmark gene sets (MSigDB)
  - GO Biological Process
  - Reactome pathways
- ✅ **Comprehensive Visualizations**:
  - Dotplots, enrichment maps, ridge plots
  - Individual pathway enrichment plots
  - Heatmaps and comparison matrices
- ✅ **Pathway Interpretation**:
  - Immune/inflammatory pathway focus
  - Disease progression analysis
  - COVID-specific pathway categories

### Specialized Features
- 🔍 **Immune Response Analysis**: Targeted analysis of immune and inflammatory pathways
- 📈 **Progression Analysis**: Pathway changes across disease severities
- 🎯 **COVID-Specific Pathways**: Focus on interferon, cytokine storm, complement, coagulation
- 📊 **Comprehensive Reporting**: Automated summary reports and statistics

## 📋 Requirements

### R Packages
```r
# Core analysis
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Visualization
library(pheatmap)
library(ggplot2)
library(plotly)
library(RColorBrewer)

# Data manipulation
library(tidyverse)
library(dplyr)

# GSEA analysis
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

# Additional utilities
library(DT)
library(knitr)
```

### Input Data
- `GSE293708_raw_counts.tsv`: Raw count matrix (genes × samples)
- Sample IDs should follow the naming convention: H1-H4, CM1-CM4, CS1-CS4, CI1-CI4

## 🔧 Installation & Setup

1. **Install R packages**:
```r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", 
                       "clusterProfiler", "enrichplot"))

# Install CRAN packages
install.packages(c("pheatmap", "ggplot2", "tidyverse", "plotly", 
                   "RColorBrewer", "msigdbr", "DT", "knitr"))
```

2. **Set up directory structure**:
```
Project/
├── GSE293708_raw_counts.tsv
├── covid_rnaseq_pipeline.R
└── Results/
    └── GSE293708_Simple_Analysis/
```

## 📖 Usage

### Quick Start
```r
# 1. Set working directory
setwd("path/to/your/project")

# 2. Run the complete pipeline
source("covid_rnaseq_pipeline.R")
```

### Step-by-Step Execution
The pipeline is organized into logical sections:

1. **Data Loading & Setup** (Lines 1-50)
2. **Preprocessing & Filtering** (Lines 51-100)
3. **Normalization & Annotation** (Lines 101-150)
4. **Quality Control Analysis** (Lines 151-250)
5. **Exploratory Data Analysis** (Lines 251-350)
6. **Differential Expression** (Lines 351-450)
7. **GSEA Preparation** (Lines 451-550)
8. **Comprehensive GSEA** (Lines 551-800)
9. **Results Interpretation** (Lines 801-1000)

## 📁 Output Files

### Primary Results
- `normalized_counts.csv`: Normalized gene expression matrix
- `sample_statistics.csv`: Per-sample quality metrics
- `top_expressed_genes.csv`: Highest expressed genes
- `key_pathways_across_severities.csv`: Pathways consistent across severities

### GSEA Results
- `gsea_results_[comparison]_[geneset].csv`: Detailed GSEA results
- `immune_inflammatory_pathways.csv`: Immune-focused pathway analysis
- `pathway_progression_analysis.csv`: Disease progression patterns
- `specific_[category]_pathways.csv`: COVID-specific pathway categories

### Visualizations
- `PCA_plot.png`: Principal component analysis
- `sample_correlation_heatmap.png`: Sample correlation matrix
- `gsea_dotplot_*.png`: GSEA enrichment dotplots
- `gsea_emapplot_*.png`: Pathway enrichment maps
- `immune_pathways_heatmap.png`: Immune pathway heatmap
- `pathway_progression_plot.png`: Disease progression visualization
- `covid_pathway_categories.png`: COVID-specific pathway analysis

### Reports
- `GSEA_Analysis_Report.txt`: Comprehensive analysis summary
- `sample_dendrogram.png`: Hierarchical clustering
- `top_variable_genes_heatmap.png`: Most variable genes

## 🔬 Analysis Highlights

### Differential Expression Comparisons
- COVID-19 Mild vs Normal
- COVID-19 Severe vs Normal  
- COVID-19 Critical vs Normal

### Key Pathway Categories Analyzed
- **Immune Response**: Interferon signaling, cytokine responses
- **Inflammation**: Cytokine storm, inflammatory responses
- **Host Defense**: Complement activation, innate immunity
- **Disease Complications**: Coagulation, hypoxia response
- **Cellular Stress**: Apoptosis, cell death pathways

### Specialized Analyses
- **Progression Analysis**: How pathways change with disease severity
- **Immune Focus**: Detailed analysis of immune/inflammatory pathways
- **COVID-Specific**: Targeted analysis of COVID-relevant biological processes

## 📊 Interpretation Guide

### Understanding Results
- **NES > 0**: Pathways upregulated in COVID vs Normal
- **NES < 0**: Pathways downregulated in COVID vs Normal
- **|NES| > 1.5**: Strong enrichment (recommended threshold)
- **p.adjust < 0.05**: Statistically significant enrichment

### Key Files to Review
1. `key_pathways_barplot.png`: Most consistently altered pathways
2. `immune_pathways_heatmap.png`: Immune response patterns
3. `pathway_progression_plot.png`: Disease severity progression
4. `GSEA_Analysis_Report.txt`: Complete analysis summary

## 🎯 Next Steps & Validation

### Recommended Follow-up
1. **Literature Validation**: Compare with published COVID-19 studies
2. **Functional Validation**: Experimental validation of key pathways
3. **Clinical Correlation**: Link pathway enrichment to clinical outcomes
4. **Drug Discovery**: Identify targetable pathways for therapeutics
5. **Network Analysis**: Pathway crosstalk and regulatory networks


## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- Gene set databases from MSigDB
- DESeq2 and clusterProfiler development teams
- COVID-19 research community
- GEO database (GSE293708)

