# DECODER-SC: Regulatory Execution Failure Analysis for Single-Cell Transcriptomes

## Overview

**DECODER-SC** is a computational **R-based algorithm** for quantifying **regulatory execution failure** in single-cell RNA-sequencing (scRNA-seq) data.
Unlike conventional transcription factor (TF) activity inference methods that focus on regulatory potential, DECODER-SC explicitly evaluates whether inferred TF activity is faithfully executed at the level of downstream target gene expression.

DECODER-SC identifies two complementary regulatory inconsistencies:
- **Regulatory execution failure** (TF ON → target gene OFF)
- **Regulatory escape** (TF OFF → target gene ON)

These events are aggregated across cluster-resolved transcriptional states to produce statistically grounded TF-level, interaction-level, and cell-level measures of transcriptional dysregulation.

---

## Key Capabilities

- Quantifies discordance between TF activity and target gene expression
- Identifies regulatory execution failure and regulatory escape at TF–gene resolution
- Computes TF-level failure scores and TF–gene failure/escape maps
- Generates per-cell regulatory violation scores for visualization
- Applies permutation-based null models for statistical inference
- Scales efficiently to large scRNA-seq datasets

---

## Method Summary

DECODER-SC operates through the following steps:

1. Cell-state definition via clustering of single-cell expression profiles
2. TF activity inference using AUCell and curated DoRothEA regulons
3. Binarization of TF activity and gene expression using data-driven thresholds
4. Detection of TF–gene regulatory inconsistencies
5. Aggregation and permutation-based statistical testing

This design introduces regulatory execution fidelity as an analytical dimension orthogonal to TF activity.

---

## Input Requirements

- Processed Seurat object (.rds)
- Normalized gene expression matrix
- TF–target regulons (e.g., DoRothEA A–C)

### Example

```r
library(Seurat)
seurat_obj <- CreateSeuratObject(counts = expression_matrix)
seurat_obj <- NormalizeData(seurat_obj)
saveRDS(seurat_obj, file = "input_seurat_object.rds")
```

---

## Output

- TF failure scores with significance estimates
- TF–gene regulatory failure and escape matrices
- Per-cell regulatory violation scores
- Publication-ready plots and tables

---

## Datasets Tested

| Dataset | Source | Description |
|--------|--------|------------|
| Lung adenocarcinoma | GEO (GSE185204) | Primary use case |
| Breast cancer | GEO | Validation |
| Pancreatic cancer | GEO | Tumor microenvironment |
| Brain tumor | GEO | Glioma |
| Hypothalamus development | GEO | Non-cancer system |

---

## Example Usage

```r
source("R/decoder_sc.R")
results <- run_DECODER_SC(
  seurat_object = seurat_obj,
  regulons = dorothea_regulons,
  activity_quantile = 0.7,
  expression_quantile = 0.7
)
```

---

## License

MIT License

---

## Contact

**Ali Haider**  
Department of Biomedical Engineering  
College of Chemistry and Life Science  
Beijing University of Technology, China  

Email: haider@emails.bjut.edu.cn

---

## Citation

If you use **DECODER-SC**, please cite the corresponding publication describing the regulatory execution failure algorithm.
