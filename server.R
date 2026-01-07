# ---- LOAD SERVER MODULES ONLY ----
library(shiny)
library(Seurat)

source("qc_server.R")
source("decoder_analysis_server.R")
source("decoder_results_server.R")

server <- function(input, output, session) {
  
  # ---------- Read uploaded RDS ----------
  seurat_obj <- reactive({
    req(input$rds_file)
    readRDS(input$rds_file$datapath)
  })
  
  # ---------- Validate object ----------
  observeEvent(seurat_obj(), {
    obj <- seurat_obj()
    if (!inherits(obj, "Seurat")) {
      showNotification("Uploaded file is NOT a Seurat object!",
                       type = "error", duration = 10)
    } else {
      showNotification("Seurat object loaded successfully",
                       type = "message", duration = 5)
    }
  })
  
  # ---------- Dataset summary ----------
  output$dataset_summary <- renderPrint({
    obj <- seurat_obj()
    req(obj)
    
    cat("===== DATASET SUMMARY =====\n\n")
    cat("Cells:", ncol(obj), "\n")
    cat("Genes:", nrow(obj), "\n")
    cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
    cat("Default assay:", DefaultAssay(obj), "\n")
    cat("Metadata columns:", length(colnames(obj@meta.data)), "\n")
  })
  
  # =====================================================
  # STEP 1: Manual QC module
  # =====================================================
  qc_result <- qc_server(
    id = "qc1",
    seurat_obj = seurat_obj
  )
  
  # =====================================================
  # FULL AUTOMATIC DECODER-SC ANALYSIS
  # =====================================================
  decoder_analysis_results <- decoder_analysis_server(
    id = "decoder",
    seurat_obj = seurat_obj
  )
  
  # =====================================================
  # RESULTS & DOWNLOADS (NO COMPUTATION)
  # =====================================================
  decoder_results_server(
    id = "results",
    analysis_results = decoder_analysis_results
  )
}
