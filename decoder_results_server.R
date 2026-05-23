# decoder_results_server.R
decoder_results_server <- function(id, analysis_results) {
  shiny::moduleServer(id, function(input, output, session) {
    
    analysis_done <- shiny::reactive({
      ar <- analysis_results()
      isTRUE(ar$done)
    })
    
    # =====================================================
    # PREVIEW: DIRECT DISCORDANCE HEATMAP
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
          title = "Direct Regulatory Discordance Map",
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
    # PREVIEW: INVERSE DISCORDANCE DOT PLOT
    # =====================================================
    shiny::observeEvent(input$preview_escape_dot, {
      if (!analysis_done()) {
        shiny::showNotification("Run analysis first.", type = "warning", duration = 4)
        return()
      }
      
      ar <- analysis_results()
      shiny::req(ar$p_inverse)
      
      shiny::showModal(
        shiny::modalDialog(
          title = "Inverse Regulatory Discordance Plot",
          size = "l",
          easyClose = TRUE,
          footer = shiny::modalButton("Close"),
          shiny::plotOutput(session$ns("plot_escape"), height = "650px")
        )
      )
    }, ignoreInit = TRUE)
    
    output$plot_escape <- shiny::renderPlot({
      analysis_results()$p_inverse
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
          title = "TF Regulatory Discordance Landscape",
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
      filename = function() "Figure_Direct_Regulatory_Discordance_Map.png",
      content = function(file) {
        file.copy(analysis_results()$heatmap_path, file, overwrite = TRUE)
      },
      contentType = "image/png"
    )
    
    output$download_escape_dot <- shiny::downloadHandler(
      filename = function() "Figure_Inverse_Regulatory_Discordance.png",
      content = function(file) {
        file.copy(analysis_results()$inverse_path, file, overwrite = TRUE)
      },
      contentType = "image/png"
    )
    
    output$download_lollipop <- shiny::downloadHandler(
      filename = function() "Figure_TF_Regulatory_Discordance.png",
      content = function(file) {
        file.copy(analysis_results()$lollipop_path, file, overwrite = TRUE)
      },
      contentType = "image/png"
    )
    
    # =====================================================
    # DOWNLOAD CSVs
    # =====================================================
    
    # TF-level discordance summary
    output$download_tf_discordance_summary <- shiny::downloadHandler(
      filename = function() "TF_regulatory_discordance_all_TFs.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_discordance_summary,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Direct Discordance (ALL)
    output$download_tf_gene_discordance_all <- shiny::downloadHandler(
      filename = function() "TF_Gene_Direct_Discordance_ALL.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_gene_discordance_all,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Inverse Discordance (ALL)
    output$download_tf_gene_inverse_discordance_all <- shiny::downloadHandler(
      filename = function() "TF_Gene_Inverse_Discordance_ALL.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_gene_inverse_discordance_all,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Direct Discordance (Top TFs)
    output$download_tf_gene_discordance_plot <- shiny::downloadHandler(
      filename = function() "TF_Gene_Direct_Discordance_TopTFs.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$direct_top_plot,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # TF–Gene Inverse Discordance (Top TFs)
    output$download_tf_gene_inverse_discordance_plot <- shiny::downloadHandler(
      filename = function() "TF_Gene_Inverse_Discordance_TopTFs.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$inverse_top_plot,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
    # =====================================================
    # NEW: TF-LEVEL PERMUTATION + FDR CSV
    # =====================================================
    output$download_tf_discordance_significance <- shiny::downloadHandler(
      filename = function() "DECODER_SC_TF_discordance_significance.csv",
      content = function(file) {
        utils::write.csv(
          analysis_results()$tf_discordance_significance,
          file,
          row.names = FALSE
        )
      },
      contentType = "text/csv"
    )
    
  })
}
