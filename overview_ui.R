overview_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    
    h1("DECODER-SC"),
    h4("Decoupling and Evaluation of Cellular Regulatory Logic in Single Cells"),
    
    hr(),
    
    h3("What does this method do?"),
    p(
      "DECODER-SC is an automated single-cell RNA-seq analysis method that quantifies",
      "regulatory rule violations at single-cell resolution. Instead of measuring only",
      "gene expression, it evaluates whether genes obey or disobey their upstream",
      "transcriptional regulators within individual cells."
    ),
    
    h3("Why is this important in single-cell studies?"),
    p(
      "Most existing single-cell analysis approaches implicitly assume intact transcriptional",
      "regulation. However, in cancer and other disease states, regulatory control is frequently",
      "disrupted. DECODER-SC directly measures this disruption by identifying cellular states",
      "where regulatory logic breaks down."
    ),
    
    h3("What kind of outputs does DECODER-SC generate?"),
    tags$ul(
      tags$li("Per-cell regulatory violation score"),
      tags$li("Identification and ranking of transcription factors with regulatory failure"),
      tags$li("Identification of genes escaping normal regulatory control"),
      tags$li("Visualization of regulatory instability across UMAP and trajectories")
    ),
    
    h3("Why is this useful for cancer research?"),
    p(
      "Cancer cells frequently activate oncogenic programs independent of normal regulatory",
      "constraints. DECODER-SC identifies tumor cells with extreme regulatory breakdown,",
      "captures transitional and unstable cellular states, and highlights genes potentially",
      "associated with aggressiveness, plasticity, or therapy resistance."
    ),
    
    h3("Dataset compatibility"),
    p(
      "DECODER-SC is dataset-agnostic and can be applied to any single-cell RNA-seq dataset,",
      "including datasets downloaded from GEO. The method operates on a processed Seurat",
      "object stored as an .rds file."
    ),
    
    h3("Prepare your input data"),
    p(
      "Please upload a Seurat object saved as an .rds file. The object should contain normalized",
      "expression data and basic quality control. Optional annotations such as cell-type labels",
      "may also be included."
    ),
    
    h3("How does the analysis run?"),
    p(
      "DECODER-SC is executed as a fully automated pipeline. After a dataset is uploaded, the user",
      "initiates the analysis once, and all internal modules execute sequentially without further",
      "interaction."
    ),
    
    tags$ul(
      tags$li("Data upload and validation"),
      tags$li("Quality control and data preparation"),
      tags$li("Regulatory reference loading"),
      tags$li("Regulatory rule violation scoring"),
      tags$li("Result aggregation and visualization")
    )
  )
}
