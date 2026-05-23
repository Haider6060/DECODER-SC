# DECODER-SC: Sign-Aware TF–Target Regulatory Discordance Analysis in Single-Cell Transcriptomes

## Overview

**DECODER-SC** is an R-based computational framework for quantifying **TF–target regulatory discordance** in single-cell RNA-sequencing (scRNA-seq) data. Unlike conventional transcription factor (TF) activity inference methods that primarily estimate regulatory potential, DECODER-SC evaluates whether inferred TF activity is concordantly reflected in downstream target gene expression.

DECODER-SC identifies two complementary forms of regulatory inconsistency:

- **Direct regulatory discordance** (TF active → target gene suppressed)
- **Inverse regulatory discordance** (TF inactive → target gene expressed)

These sign-aware TF–target mismatches are aggregated across cellular states to generate statistically grounded TF-level, interaction-level, and cell-level measures of regulatory inconsistency.

---

# Key Features

- Quantifies sign-aware TF–target regulatory discordance
- Detects both direct and inverse regulatory inconsistencies
- Generates TF-level discordance scores and TF–gene discordance maps
- Computes cell-level regulatory violation landscapes
- Performs permutation-based statistical significance testing with FDR correction
- Produces publication-ready visualizations and output tables
- Supports large-scale scRNA-seq datasets and cross-dataset analysis

---

# Method Summary

DECODER-SC operates through the following analytical steps:

1. Quality control and preprocessing of scRNA-seq data
2. TF activity inference using AUCell and curated DoRothEA regulons
3. Quantile-based discretization of TF activity and gene expression states
4. Sign-aware detection of TF–target regulatory discordance
5. Aggregation of discordance scores across TFs, genes, and cells
6. Permutation-based statistical inference and FDR correction

This framework introduces **regulatory execution fidelity** as an orthogonal analytical dimension beyond conventional TF activity estimation.

---

# Input Requirements

- Processed Seurat object (`.rds`)
- Normalized single-cell gene expression matrix
- TF regulons (recommended: DoRothEA confidence levels A–C)

### Example

```r
library(Seurat)

seurat_obj <- CreateSeuratObject(counts = expression_matrix)
seurat_obj <- NormalizeData(seurat_obj)

saveRDS(seurat_obj, file = "input_seurat_object.rds")
```

---

# Output

DECODER-SC generates:

- TF-level regulatory discordance scores
- TF–gene direct discordance maps
- TF–gene inverse discordance maps
- Cell-level regulatory violation scores
- Permutation-based significance estimates and FDR-adjusted results
- Publication-ready figures and summary tables

---

# Datasets Evaluated

| Dataset | Source | Description |
|---|---|---|
| Lung adenocarcinoma | GEO (GSE185204) | Primary discovery dataset |
| Independent lung cancer cohort | GEO | Cross-cohort validation |
| Breast cancer | GEO | Generalizability analysis |
| Pancreatic cancer | GEO | Tumor microenvironment analysis |
| Brain tumor | GEO | Glioma regulatory analysis |
| Developmental hypothalamus | GEO (GSE132355) | Non-cancer developmental validation |

---

# Example Usage

```r
results <- run_DECODER_SC(
  seurat_object = seurat_obj,
  regulons = dorothea_regulons,
  activity_quantile = 0.7,
  expression_quantile = 0.7
)
```

---

# Installation

```r
# Required packages
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("AUCell", "dorothea"))
```

---

# License

MIT License

---

# Contact

**Ali Haider**  
Department of Biomedical Engineering  
College of Chemistry and Life Sciences  
Beijing University of Technology, China  

Email: haider@emails.bjut.edu.cn

---

# Citation

If you use **DECODER-SC**, please cite the corresponding publication describing the DECODER-SC framework for TF–target regulatory discordance analysis in single-cell transcriptomics.