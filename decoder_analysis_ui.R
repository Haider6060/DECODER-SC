decoder_analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    
    shiny::h1("Run DECODER-SC (Automatic)"),
    
    shiny::p(
      "This module runs the full pipeline automatically: QC → clustering → TF activity → sign-aware regulatory discordance scoring → UMAP → saves PNG/CSV."
    ),
    
    # ============================
    # NEW: EXPLANATION SECTION
    # ============================
    shiny::h3("Sign-aware regulatory discordance analysis"),
    
    shiny::p(
      "DECODER-SC is designed to identify regulatory inconsistency and transcription factor–target discordance at single-cell resolution. 
       It quantifies transcription factor–target regulatory discordance across single-cell states by identifying activity-expression mismatches between TFs and their target genes."
    ),
    
    shiny::tags$ul(
      shiny::tags$li(
        shiny::strong("Direct regulatory discordance: "),
        "For activator TFs, discordance occurs when the TF is active but the target gene is not expressed. 
   For repressor TFs, discordance occurs when both the TF and target gene are simultaneously active."
      ),
      shiny::tags$li(
        shiny::strong("Inverse regulatory discordance: "),
        "For activator TFs, inverse discordance occurs when the TF is inactive but the target gene is expressed. 
   For repressor TFs, inverse discordance occurs when both the TF and target gene are simultaneously inactive."
      ),
      shiny::tags$li(
        shiny::strong("Discordance score: "),
        "A per-cell quantitative measure summarizing regulatory inconsistency across all TF–gene rules."
      )
    ),
    
    shiny::p(
      "By operating on cluster-averaged single-cell states, DECODER-SC scales to very large datasets 
       and highlights regulatory inconsistency that is not detectable using expression alone."
    ),
    
    shiny::p(
      shiny::strong("Why this matters in single-cell analysis: "),
      "Regulatory discordance patterns reveal transcriptional dysregulation, cellular plasticity, 
       and hidden regulatory rewiring associated with disease progression, cell-state transitions, 
       and tumor heterogeneity."
    ),
    
    shiny::hr(),
    
    # ============================
    # EXISTING CONTENT (UNCHANGED)
    # ============================
    shiny::h3("Outputs"),
    
    shiny::tags$ul(
      shiny::tags$li("TF regulatory discordance summary (CSV)"),
      shiny::tags$li("Direct regulatory discordance heatmap (PNG)"),
      shiny::tags$li("Inverse regulatory discordance dot plot (PNG)"),
      shiny::tags$li("TF regulatory discordance significance plot (PNG)"),
      shiny::tags$li("TF–gene discordance score tables (CSV)")
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
        shiny::numericInput(ns("tf_quantile"), "TF activity binarization threshold (quantile)",
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
