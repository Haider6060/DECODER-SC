# ---- LOAD UI MODULES ONLY ----
library(shiny)

source("overview_ui.R")
source("qc_ui.R")
source("decoder_analysis_ui.R")
source("decoder_results_ui.R")

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      body {
        background-color: #1f2933;
        font-family: Arial, Helvetica, sans-serif;
      }
      .sidebar {
        background-color: #111827;
        height: 100vh;
        padding-top: 25px;
      }
      .sidebar h2 {
        color: #e5e7eb;
        text-align: center;
        margin-bottom: 30px;
        font-weight: 700;
      }
      .sidebar .radio label {
        color: #d1d5db;
        font-size: 15px;
        padding: 8px;
      }
      .content-box {
        background-color: #f0f4f8;
        padding: 35px;
        border-radius: 16px;
        margin: 30px;
      }
      h1 { color: #0f172a; font-weight: 800; }
      h4 { color: #334155; }
      h3 { color: #1e3a8a; margin-top: 25px; }
      p, li { color: #1f2937; font-size: 16px; line-height: 1.7; }
      pre { background-color: #e5eaf0; padding: 15px; border-radius: 10px; font-size: 14px; }
    "))
  ),
  
  fluidRow(
    column(
      width = 3,
      class = "sidebar",
      h2("DECODER-SC"),
      
      radioButtons(
        inputId = "module_select",
        label = NULL,
        choices = c(
          "Method Overview" = "overview",
          "Data Upload" = "upload",
          "QC & Data Preparation" = "step1",
          "Run Full Analysis" = "run",
          "Results & Downloads" = "results"
        )
      )
    ),
    
    column(
      width = 9,
      
      conditionalPanel(
        condition = "input.module_select == 'overview'",
        div(class = "content-box", overview_ui("overview"))
      ),
      
      conditionalPanel(
        condition = "input.module_select == 'upload'",
        div(
          class = "content-box",
          h1("Data Upload"),
          p("Upload a processed Seurat object (.rds). This dataset will be used for all downstream analyses."),
          fileInput("rds_file", "Upload Seurat Object (.rds)", accept = ".rds"),
          h3("Seurat Object Information"),
          verbatimTextOutput("dataset_summary")
        )
      ),
      
      conditionalPanel(
        condition = "input.module_select == 'step1'",
        div(class = "content-box", qc_ui("qc1"))
      ),
      
      # Run Full Analysis = Decoder Analysis module UI (automatic QC + analysis)
      conditionalPanel(
        condition = "input.module_select == 'run'",
        div(class = "content-box", decoder_analysis_ui("decoder"))
      ),
      
      conditionalPanel(
        condition = "input.module_select == 'results'",
        div(class = "content-box", decoder_results_ui("results"))
      )
    )
  )
)
