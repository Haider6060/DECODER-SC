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
        shiny::actionButton(ns("preview_heatmap"), "Preview Direct Regulatory Discordance Map"),
        shiny::br(), shiny::br(),
        shiny::downloadButton(
          ns("download_heatmap"),
          "Download Direct Regulatory Discordance Map (PNG)"
        )
      ),
      shiny::column(
        6,
        shiny::actionButton(ns("preview_escape_dot"), "Preview Inverse Regulatory Discordance Plot"),
        shiny::br(), shiny::br(),
        shiny::downloadButton(
          ns("download_escape_dot"),
          "Download Inverse Regulatory Discordance Plot (PNG)"
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
          "Preview TF Discordance Significance Plot"
        ),
        shiny::br(), shiny::br(),
        shiny::downloadButton(
          ns("download_lollipop"),
          "Download TF Discordance Significance Plot (PNG)"
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
          ns("download_tf_discordance_summary"),
          "TF Regulatory Discordance Scores (TF-level CSV)"
        )
      ),
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_discordance_all"),
          "TF–Gene Direct Discordance (ALL pairs CSV)"
        )
      )
    ),
    
    shiny::br(),
    
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_inverse_discordance_all"),
          "TF–Gene Inverse Discordance (ALL pairs CSV)"
        )
      ),
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_discordance_plot"),
          "TF–Gene Direct Discordance (Top TFs CSV)"
        )
      )
    ),
    
    shiny::br(),
    
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::downloadButton(
          ns("download_tf_gene_inverse_discordance_plot"),
          "TF–Gene Inverse Discordance (Top TFs CSV)"
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
          ns("download_tf_discordance_significance"),
          "TF Regulatory Discordance Significance (Permutation + FDR CSV)"
        )
      )
    )
    
  )
}
