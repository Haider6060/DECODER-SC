# decoder_results_server.R
decoder_results_server <- function(id, analysis_results) {
  shiny::moduleServer(id, function(input, output, session) {
    
    analysis_done <- shiny::reactive({
      ar <- analysis_results()
      isTRUE(ar$done)
    })
    
    # =====================================================
    # PREVIEW: FAILURE HEATMAP
    # =====================================================
    shiny::observeEvent(input$preview_heatmap, {
      if (!analysis_done()) {
        shiny::showNotification("Run analysis first.", type = "warning", duration = 4)
        return()
      }
      
      ar <- analysis_results()
      shiny::req(ar$p_heat)
      
      shiny::showModal(
        shiny::modalDialog(
          title = "Failure Heatmap (Top TFs)",
          size = "l",
          easyClose = TRUE,
          footer = shiny::modalButton("Close"),
          shiny::plotOutput(session$ns("plot_heat"), height = "650px")
        )
      )
    }, ignoreInit = TRUE)
    
    output$plot_heat <- shiny::renderPlot({
      analysis_results()$p_heat
    })
    
    # =====================================================
    # PREVIEW: ESCAPE DOT PLOT
    # =====================================================
    shiny::observeEvent(input$preview_escape_dot, {
      if (!analysis_done()) {
        shiny::showNotification("Run analysis first.", type = "warning", duration = 4)
        return()
      }
      
      ar <- analysis_results()
      shiny::req(ar$p_escape)
      
      shiny::showModal(
        shiny::modalDialog(
          title = "Escape Dot Plot (Top TFs)",
          size = "l",
          easyClose = TRUE,
          footer = shiny::modalButton("Close"),
          shiny::plotOutput(session$ns("plot_escape"), height = "650px")
        )
      )
    }, ignoreInit = TRUE)
    
    output$plot_escape <- shiny::renderPlot({
      analysis_results()$p_escape
    })
    
    # =====================================================
    # NEW: PREVIEW LOLLIPOP PLOT
    # =====================================================
    shiny::observeEvent(input$preview_lollipop, {
      if (!analysis_done()) {
        shiny::showNotification("Run analysis first.", type = "warning", duration = 4)
        return()
      }
      
      ar <- analysis_results()
      shiny::req(ar$p_lollipop)
      
      shiny::showModal(
        shiny::modalDialog(
          title = "TF Failure Lollipop Plot",
          size = "l",
          easyClose = TRUE,
          footer = shiny::modalButton("Close"),
          shiny::plotOutput(session$ns("plot_lollipop"), height = "650px")
        )
      )
    }, ignoreInit = TRUE)
    
    output$plot_lollipop <- shiny::renderPlot({
      analysis_results()$p_lollipop
    })
    
    # =====================================================
    # DOWNLOAD PNGs
    # =====================================================
    output$download_heatmap <- shiny::downloadHandler(
      filename = function() "Figure_Failure_Heatmap_TopTFs.png",
      content = function(file) {
        file.copy(analysis_results()$heatmap_path, file, overwrite = TRUE)
      },
      contentType = "image/png"
    )
    
    output$download_escape_dot <- shiny::downloadHandler(
      filename = function() "Figure_Escape_DotPlot_TopTFs.png",
      content = function(file) {
        file.copy(analysis_results()$escape_path, file, overwrite = TRUE)
      },
      contentType = "image/png"
    )
    
    output$download_lollipop <- shiny::downloadHandler(
      filename = function() "Figure_TF_Failure_Lollipop.png",
      content = function(file) {
        file.copy(analysis_results()$lollipop_path, file, overwrite = TRUE)
      },
      contentType = "image/png"
    )
    
    # =====================================================
    # DOWNLOAD CSVs
    # =====================================================
    
    # TF-level failure summary
    output$download_tf_failure_summary <- shiny::downloadHandler(
      filename = function() "TF_failure_frequency_all_TFs.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_failure_summary,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Failure (ALL)
    output$download_tf_gene_failure_all <- shiny::downloadHandler(
      filename = function() "TF_Gene_Failure_ALL.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_gene_failure_all,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Escape (ALL)
    output$download_tf_gene_escape_all <- shiny::downloadHandler(
      filename = function() "TF_Gene_Escape_ALL.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_gene_escape_all,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Failure (Top TFs)
    output$download_tf_gene_failure_plot <- shiny::downloadHandler(
      filename = function() "TF_Gene_Failure_TopTFs_NONZERO.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_gene_failure_all[
            analysis_results()$tf_gene_failure_all$tf_gene_failure_freq > 0, ],
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Escape (Top TFs)
    output$download_tf_gene_escape_plot <- shiny::downloadHandler(
      filename = function() "TF_Gene_Escape_TopTFs_NONZERO.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_gene_escape_all[
            analysis_results()$tf_gene_escape_all$tf_gene_escape_freq > 0, ],
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # =====================================================
    # NEW: TF-LEVEL PERMUTATION + FDR CSV
    # =====================================================
    output$download_tf_failure_significance <- shiny::downloadHandler(
      filename = function() "DECODER_SC_TF_failure_significance.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_failure_significance,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
  })
}
