qc_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h1("QC & Data Preparation"),
    
    p(
      "This module performs quality control and data preparation automatically as part of the",
      "DECODER-SC pipeline. No manual interaction is required."
    ),
    
    tags$ul(
      tags$li("Cell and gene filtering"),
      tags$li("Mitochondrial content calculation"),
      tags$li("Normalization"),
      tags$li("Variable feature selection"),
      tags$li("Scaling and PCA")
    ),
    
    h3("QC status"),
    verbatimTextOutput(ns("qc_status"))
  )
}
