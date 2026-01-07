# decoder_results_ui.R
decoder_results_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    
    shiny::h1("Results & Downloads"),
    
    # =========================
    # FIGURES
    # =========================
    shiny::h3("Figures"),
    
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::actionButton(ns("preview_heatmap"), "Preview Failure Heatmap"),
        shiny::br(), shiny::br(),
        shiny::downloadButton(
          ns("download_heatmap"),
          "Download Failure Heatmap (PNG)"
        )
      ),
      shiny::column(
        6,
        shiny::actionButton(ns("preview_escape_dot"), "Preview Escape Dot Plot"),
        shiny::br(), shiny::br(),
        shiny::downloadButton(
          ns("download_escape_dot"),
          "Download Escape Dot Plot (PNG)"
        )
      )
    ),
    
    shiny::br(),
    
    # =========================
    # NEW: LOLLIPOP PLOT
    # =========================
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::actionButton(
          ns("preview_lollipop"),
          "Preview TF Failure Lollipop Plot"
        ),
        shiny::br(), shiny::br(),
        shiny::downloadButton(
          ns("download_lollipop"),
          "Download TF Failure Lollipop (PNG)"
        )
      )
    ),
    
    shiny::hr(),
    
    # =========================
    # TABLES
    # =========================
    shiny::h3("Tables (CSV)"),
    
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_failure_summary"),
          "TF Failure Summary (TF-level CSV)"
        )
      ),
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_failure_all"),
          "TF–Gene Failure (ALL pairs CSV)"
        )
      )
    ),
    
    shiny::br(),
    
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_escape_all"),
          "TF–Gene Escape (ALL pairs CSV)"
        )
      ),
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_failure_plot"),
          "TF–Gene Failure (Top TFs, Plot CSV)"
        )
      )
    ),
    
    shiny::br(),
    
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_escape_plot"),
          "TF–Gene Escape (Top TFs, Plot CSV)"
        )
      )
    ),
    
    shiny::br(),
    
    # =========================
    # NEW: TF-LEVEL SIGNIFICANCE CSV
    # =========================
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_failure_significance"),
          "TF Failure Significance (Permutation + FDR CSV)"
        )
      )
    )
    
  )
}
