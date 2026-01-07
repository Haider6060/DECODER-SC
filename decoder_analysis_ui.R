decoder_analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    
    shiny::h1("Run DECODER-SC (Automatic)"),
    
    shiny::p(
      "This module runs the full pipeline automatically: QC → clustering → TF activity → violation/failure/escape scoring → UMAP → saves PNG/CSV."
    ),
    
    # ============================
    # NEW: EXPLANATION SECTION
    # ============================
    shiny::h3("What does DECODER-SC do?"),
    
    shiny::p(
      "DECODER-SC is designed to identify dysfunctional transcriptional regulation at single-cell resolution. 
       It quantifies how often transcription factors (TFs) fail to control their target genes or how genes escape TF regulation across cell states."
    ),
    
    shiny::tags$ul(
      shiny::tags$li(
        shiny::strong("TF failure (TF ON → Gene OFF): "),
        "The transcription factor is active, but its target gene is not expressed. 
         This indicates broken or ineffective regulation."
      ),
      shiny::tags$li(
        shiny::strong("TF escape (TF OFF → Gene ON): "),
        "The transcription factor is inactive, but the target gene is expressed. 
         This suggests alternative or compensatory regulatory mechanisms."
      ),
      shiny::tags$li(
        shiny::strong("Violation score: "),
        "A per-cell quantitative measure summarizing regulatory inconsistency across all TF–gene rules."
      )
    ),
    
    shiny::p(
      "By operating on cluster-averaged single-cell states, DECODER-SC scales to very large datasets 
       and highlights regulatory instability that is not detectable using expression alone."
    ),
    
    shiny::p(
      shiny::strong("Why this matters in single-cell analysis: "),
      "Regulatory failures and escapes reveal transcriptional dysregulation, cellular plasticity, 
       and hidden regulatory rewiring associated with disease progression, cell-state transitions, 
       and tumor heterogeneity."
    ),
    
    shiny::hr(),
    
    # ============================
    # EXISTING CONTENT (UNCHANGED)
    # ============================
    shiny::h3("Outputs"),
    
    shiny::tags$ul(
      shiny::tags$li("Per-cell violation score (UMAP)"),
      shiny::tags$li("Top TF failure summary (CSV)"),
      shiny::tags$li("TF–gene failure heatmap + dot plot (PNG)"),
      shiny::tags$li("TF–gene escape dot plot (PNG)")
    ),
    
    shiny::hr(),
    
    shiny::textInput(
      ns("out_dir"),
      "Output directory (PNG and CSV files will be saved here)",
      value = "DECODER_SC_outputs"
    ),
    
    shiny::fluidRow(
      shiny::column(
        4,
        shiny::numericInput(ns("top_n_tf"), "Number of top TFs", value = 20, min = 5, max = 100)
      ),
      shiny::column(
        4,
        shiny::numericInput(ns("top_k_genes"), "Top genes per TF", value = 5, min = 3, max = 50)
      ),
      shiny::column(
        4,
        shiny::numericInput(ns("tf_quantile"), "Activity threshold (quantile)",
                            value = 0.75, min = 0.5, max = 0.95, step = 0.05)
      )
    ),
    
    shiny::checkboxGroupInput(
      ns("dorothea_conf"),
      "DoRothEA confidence levels",
      choices = c("A","B","C","D","E"),
      selected = c("A","B","C"),
      inline = TRUE
    ),
    
    shiny::br(),
    
    shiny::actionButton(
      ns("run_analysis"),
      "Run DECODER-SC Analysis",
      class = "btn-primary"
    ),
    
    shiny::br(), shiny::br(),
    shiny::verbatimTextOutput(ns("run_status"))
  )
}
