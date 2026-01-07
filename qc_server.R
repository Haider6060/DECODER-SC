qc_server <- function(id, seurat_obj) {
  moduleServer(id, function(input, output, session) {
    
    qc_result <- reactive({
      
      req(seurat_obj())
      
      obj <- seurat_obj()
      
      # Add mitochondrial percentage (generic)
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
      
      # Generic QC filtering
      obj <- subset(
        obj,
        subset =
          nFeature_RNA > 200 &
          nFeature_RNA < 6000 &
          percent.mt < 20
      )
      
      # Normalization & scaling
      obj <- NormalizeData(obj)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
      obj <- ScaleData(obj)
      obj <- RunPCA(obj, features = VariableFeatures(obj))
      
      obj
    })
    
    output$qc_status <- renderPrint({
      req(qc_result())
      
      obj <- qc_result()
      
      cat("QC & Data Preparation completed successfully\n\n")
      cat("Number of cells:", ncol(obj), "\n")
      cat("Number of genes:", nrow(obj), "\n")
      cat("Default assay:", DefaultAssay(obj), "\n")
    })
    
    return(qc_result)
  })
}
