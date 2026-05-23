decoder_analysis_server <- function(id, seurat_obj) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # ---------------------------
    # STATUS LOGGER (UI)
    # ---------------------------
    status_lines <- shiny::reactiveVal(character())
    add_status <- function(msg) {
      status_lines(c(
        status_lines(),
        paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
      ))
    }
    
    output$run_status <- shiny::renderText({
      paste(status_lines(), collapse = "\n")
    })
    
    # ---------------------------
    # RESULTS CONTAINER
    # ---------------------------
    res <- shiny::reactiveValues(
      done = FALSE,
      p_heat = NULL,
      p_inverse = NULL,
      tf_discordance_summary = NULL,
      tf_gene_discordance_all = NULL,
      tf_gene_inverse_discordance_all = NULL,
      heatmap_path = NULL,
      inverse_path  = NULL,
      p_lollipop = NULL,
      tf_discordance_significance = NULL,
      lollipop_path = NULL,
      tf_sig_csv_path = NULL
    )
    
    # =====================================================
    # RUN ANALYSIS
    # =====================================================
    shiny::observeEvent(input$run_analysis, {
      
      status_lines(character())
      res$done <- FALSE
      add_status("⏳ Analysis started")
      message("DECODER-SC: analysis started")
      
      tryCatch({
        shiny::withProgress(
          message = "Running DECODER-SC analysis",
          value = 0, {
            
            # ---------------------------
            # LOAD OBJECT
            # ---------------------------
            shiny::incProgress(0.05, detail = "Loading Seurat object")
            message("Loading Seurat object")
            obj <- seurat_obj()
            shiny::req(obj)
            
            out_dir <- input$out_dir
            if (out_dir == "") out_dir <- "DECODER_SC_outputs"
            if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
            
            # ---------------------------
            # QC
            # ---------------------------
            shiny::incProgress(0.15, detail = "QC and normalization")
            message("Running QC and normalization")
            if (!"percent.mt" %in% colnames(obj@meta.data)) {
              obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, "^MT-")
            }
            keep.cells <- WhichCells(
              obj,
              expression = nFeature_RNA > 200 & percent.mt < 20
            )
            
            obj <- obj[, keep.cells]
            rm(keep.cells)
            gc()
            obj <- NormalizeData(obj, verbose = FALSE)
            obj <- FindVariableFeatures(obj, verbose = FALSE)
            obj <- ScaleData(
              obj,
              features = VariableFeatures(obj),
              verbose = FALSE
            )
            obj <- RunPCA(obj, verbose = FALSE)
            
            # ---------------------------
            # CLUSTERING
            # ---------------------------
            shiny::incProgress(0.30, detail = "Clustering cells")
            message("Clustering cells")
            obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
            obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
            
            # ---------------------------
            # TF ACTIVITY
            # ---------------------------
            shiny::incProgress(0.45, detail = "Inferring TF activity")
            message("Inferring TF activity")
            regulons_df <- dorothea::dorothea_hs[
              dorothea::dorothea_hs$confidence %in% input$dorothea_conf, ]
            regulons <- split(regulons_df$target, regulons_df$tf)
            
            Idents(obj) <- obj$seurat_clusters
            avg_expr <- Seurat::AggregateExpression(
              obj,
              group.by = "seurat_clusters",
              assays = "RNA",
              slot = "data"
            )$RNA
            rm(obj)
            gc()
            rankings <- AUCell::AUCell_buildRankings(avg_expr, plotStats = FALSE)
            tf_activity <- AUCell::AUCell_calcAUC(regulons, rankings)
            tf_activity_mat <- AUCell::getAUC(tf_activity)
            
            # ---------------------------
            # DIRECT / INVERSE REGULATORY DISCORDANCE
            # ---------------------------
            shiny::incProgress(0.65, detail = "Computing regulatory discordance")
            message("Computing regulatory discordance")
            common_genes <- intersect(rownames(avg_expr), regulons_df$target)
            state_expr <- avg_expr[common_genes, , drop = FALSE]
            
            rules <- regulons_df[regulons_df$target %in% common_genes, ]
            rules$gene <- rules$target
            
            q <- input$tf_quantile
            tf_state   <- tf_activity_mat > apply(tf_activity_mat, 1, quantile, probs = q)
            gene_state <- state_expr > apply(state_expr, 1, quantile, probs = q)
            gene_z <- t(scale(t(state_expr)))
            tf_discordance_mat <- matrix(0, nrow = nrow(rules), ncol = ncol(state_expr))
            gene_inverse_discordance_mat <- tf_discordance_mat
            
            for (i in seq_len(nrow(rules))) {
              
              tf  <- rules$tf[i]
              tg  <- rules$target[i]
              mor <- rules$mor[i]
              
              if (tf %in% rownames(tf_state)) {
                
                tf_on   <- tf_state[tf, ]
                gene_on <- gene_state[tg, ]
                
                # =====================================================
                # ACTIVATOR TF
                # mor = +1
                # Expected:
                # TF ON  -> Gene ON
                # TF OFF -> Gene OFF
                # =====================================================
                
                if (mor == 1) {
                  
                  # Discordance:
                  # TF ON + Gene OFF
                  
                  observed <- ifelse(
                    tf_on & !gene_on,
                    abs(gene_z[tg, ]),
                    0
                  )
                  
                  expected <- mean(tf_on) * mean(!gene_on)
                  
                  tf_discordance_mat[i, ] <- pmax(observed - expected, 0)
                  
                  # Inverse discordance:
                  # TF OFF + Gene ON
                  
                  observed <- ifelse(
                    !tf_on & gene_on,
                    abs(gene_z[tg, ]),
                    0
                  )
                  
                  expected <- mean(!tf_on) * mean(gene_on)
                  
                  gene_inverse_discordance_mat[i, ] <- pmax(observed - expected, 0)
                }
                
                # =====================================================
                # REPRESSOR TF
                # mor = -1
                # Expected:
                # TF ON  -> Gene OFF
                # TF OFF -> Gene ON
                # =====================================================
                
                if (mor == -1) {
                  
                  # Discordance:
                  # TF ON + Gene ON
                  
                  observed <- ifelse(
                    tf_on & gene_on,
                    abs(gene_z[tg, ]),
                    0
                  )
                  
                  expected <- mean(tf_on) * mean(gene_on)
                  
                  tf_discordance_mat[i, ] <- pmax(observed - expected, 0)
                  
                  # Inverse discordance:
                  # TF OFF + Gene OFF
                  
                  observed <- ifelse(
                    !tf_on & !gene_on,
                    abs(gene_z[tg, ]),
                    0
                  )
                  
                  expected <- mean(!tf_on) * mean(!gene_on)
                  
                  gene_inverse_discordance_mat[i, ] <- pmax(observed - expected, 0)
                }
              }
            }
            
            rules$tf_gene_direct_discordance_score <- rowMeans(tf_discordance_mat)
            rules$tf_gene_inverse_discordance_score  <- rowMeans(gene_inverse_discordance_mat)
            res$tf_gene_discordance_all <- rules[, c(
              "tf",
              "confidence",
              "gene",
              "tf_gene_direct_discordance_score"
            )]
            
            res$tf_gene_inverse_discordance_all <- rules[, c(
              "tf",
              "confidence",
              "gene",
              "tf_gene_inverse_discordance_score"
            )]
            # ---------------------------
            # SAVE CSVs
            # ---------------------------
            shiny::incProgress(0.80, detail = "Saving result tables")
            message("Saving CSV outputs")
            write.csv(
              rules[, c("tf","confidence","gene","tf_gene_direct_discordance_score")],
              file.path(out_dir, "TF_Gene_Direct_Discordance_ALL.csv"),
              row.names = FALSE
            )
            
            write.csv(
              rules[, c("tf","confidence","gene","tf_gene_inverse_discordance_score")],
              file.path(out_dir, "TF_Gene_Inverse_Discordance_ALL.csv"),
              row.names = FALSE
            )
            
            # ---------------------------
            # TF RANKINGS
            # ---------------------------
            shiny::incProgress(0.85, detail = "Ranking transcription factors")
            message("Ranking transcription factors")
            tf_discordance_rank <- aggregate(tf_gene_direct_discordance_score ~ tf, rules, mean)
            tf_discordance_rank <- tf_discordance_rank[order(-tf_discordance_rank$tf_gene_direct_discordance_score), ]
            
            tf_inverse_discordance_rank <- aggregate(tf_gene_inverse_discordance_score ~ tf, rules, mean)
            tf_inverse_discordance_rank <- tf_inverse_discordance_rank[order(-tf_inverse_discordance_rank$tf_gene_inverse_discordance_score), ]
            
            res$tf_discordance_summary <- tf_discordance_rank
            
            write.csv(
              tf_discordance_rank,
              file.path(out_dir, "TF_regulatory_discordance_summary.csv"),
              row.names = FALSE
            )
            
            # =====================================================
            # LIMIT TFs FOR PLOTS ONLY
            # =====================================================
            plot_top_n_tf <- min(input$top_n_tf, 12)
            top_tfs_direct   <- head(tf_discordance_rank$tf,   plot_top_n_tf)
            top_tfs_inverse <- head(tf_inverse_discordance_rank$tf, plot_top_n_tf)
            
            # ---------------------------
            # TOP 5 GENES PER TF
            # ---------------------------
            shiny::incProgress(0.90, detail = "Preparing plot data")
            message("Preparing plot data")
            direct_data <- rules[
              rules$tf %in% top_tfs_direct &
                rules$tf_gene_direct_discordance_score > 0,
            ]
            
            inverse_data <- rules[
              rules$tf %in% top_tfs_inverse &
                rules$tf_gene_inverse_discordance_score > 0,
            ]
            
            direct_top <- direct_data |>
              dplyr::group_by(tf) |>
              dplyr::arrange(
                dplyr::desc(tf_gene_direct_discordance_score),
                .by_group = TRUE
              ) |>
              dplyr::slice_head(n = 10) |>
              dplyr::ungroup()
            
            inverse_top <- inverse_data |>
              dplyr::group_by(tf) |>
              dplyr::arrange(
                dplyr::desc(tf_gene_inverse_discordance_score),
                .by_group = TRUE
              ) |>
              dplyr::slice_head(n = 10) |>
              dplyr::ungroup()
            direct_top <- direct_top[!is.na(direct_top$tf) & !is.na(direct_top$gene), ]
            inverse_top  <- inverse_top[!is.na(inverse_top$tf)  & !is.na(inverse_top$gene), ]
            
            direct_top$tf <- droplevels(factor(direct_top$tf))
            inverse_top$tf  <- droplevels(factor(inverse_top$tf))
            
            write.csv(
              direct_top,
              file.path(out_dir, "TF_Gene_Direct_Discordance_TopTFs_TOP5.csv"),
              row.names = FALSE
            )
            
            write.csv(
              inverse_top,
              file.path(out_dir, "TF_Gene_Inverse_Discordance_TopTFs_TOP5.csv"),
              row.names = FALSE
            )
            
            # ---------------------------
            # PERMUTATION + FDR
            # ---------------------------
            shiny::incProgress(0.93, detail = "Permutation testing")
            message("Running permutation testing")
            tf_discordance_score <- tapply(rowMeans(tf_discordance_mat), rules$tf, mean)
            tf_discordance_score <- sort(tf_discordance_score, decreasing = TRUE)
            
            set.seed(1)
            n_perm <- 2000
            
            tf_names <- names(tf_discordance_score)
            perm_p <- setNames(rep(NA_real_, length(tf_names)), tf_names)
            
            for (tf in tf_names) {
              idx <- which(rules$tf == tf)
              if (length(idx) < 10) next
              observed <- mean(rowMeans(tf_discordance_mat[idx, , drop = FALSE]))
              perm_vals <- replicate(n_perm, {
                samp <- sample(seq_len(nrow(tf_discordance_mat)), length(idx), replace = FALSE)
                mean(rowMeans(tf_discordance_mat[samp, , drop = FALSE]))
              })
              perm_p[tf] <- mean(perm_vals >= observed)
            }
            
            perm_p[is.na(perm_p)] <- 1
            perm_fdr <- p.adjust(perm_p, method = "BH")
            
            tf_sig <- data.frame(
              TF = tf_names,
              TF_discordance_score = as.numeric(tf_discordance_score),
              permutation_p = as.numeric(perm_p[tf_names]),
              FDR = as.numeric(perm_fdr[tf_names])
            )
            
            tf_sig <- tf_sig[order(tf_sig$FDR, -tf_sig$TF_discordance_score), ]
            res$tf_discordance_significance <- tf_sig
            
            write.csv(
              tf_sig,
              file.path(out_dir, "DECODER_SC_tf_discordance_significance.csv"),
              row.names = FALSE
            )
            
            # ---------------------------
            # GENERATE PLOTS
            # ---------------------------
            shiny::incProgress(0.97, detail = "Generating plots")
            message("Generating plots")
            
            # ==========================================================
            # DIRECT REGULATORY DISCORDANCE HEATMAP (TF on X, Gene on Y)
            # ==========================================================
            # Order factors for clean display
            
            direct_top_plot <- direct_top
            direct_top_plot$tf <- factor(
              direct_top_plot$tf,
              levels = unique(direct_top_plot$tf)
            )
            
            direct_top_plot$gene <- factor(
              direct_top_plot$gene,
              levels = rev(unique(direct_top_plot$gene))
            )
            
            # ---------------------------
            # HEATMAP (FACETED, SAFE SIZE)
            # ---------------------------
            n_tf <- length(unique(direct_top_plot$tf))
            heat_h <- max(6, min(14, 3.5 * n_tf))
            
            res$p_heat <- ggplot2::ggplot(
              direct_top_plot,
              ggplot2::aes(
                x = "Discordance",
                y = gene,
                fill = tf_gene_direct_discordance_score
              )
            ) +
              ggplot2::geom_tile(
                color = "white",
                linewidth = 0.5,
                width = 0.95,
                height = 0.95
              ) +
              ggplot2::facet_wrap(~ tf, ncol = 3, scales = "free_y") +
              ggplot2::scale_fill_gradient(low = "white", high = "red") +
              ggplot2::theme_minimal(base_size = 12) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                panel.grid   = ggplot2::element_blank(),
                panel.spacing.x = unit(0.15, "lines"),
                strip.text   = ggplot2::element_text(face = "bold", size = 11),
                strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70")
              ) +
              ggplot2::labs(
                title = "Direct Regulatory Discordance",
                y = "Target genes"
              )
            
            res$heatmap_path <- file.path(out_dir, "Figure_Direct_Regulatory_Discordance.png")
            ggplot2::ggsave(
              res$heatmap_path,
              res$p_heat,
              width  = 12,
              height = heat_h,
              dpi    = 600
            )
            
            # ==========================================================
            # INVERSE REGULATORY DISCORDANCE DOTPLOT (TF on X, Gene on Y)
            # ==========================================================
            inverse_top_plot <- inverse_top
            
            inverse_top_plot$tf <- factor(
              inverse_top_plot$tf,
              levels = unique(inverse_top_plot$tf)
            )
            
            inverse_top_plot$gene <- factor(
              inverse_top_plot$gene,
              levels = rev(unique(inverse_top_plot$gene))
            )
            
            # ---------------------------
            # INVERSE REGULATORY DISCORDANCE DOTPLOT
            # ---------------------------
            n_tf <- length(unique(inverse_top_plot$tf))
            esc_h <- max(6, min(14, 3.5 * n_tf))
            
            res$p_inverse <- ggplot2::ggplot(
              inverse_top_plot,
              ggplot2::aes(
                x = 1,
                y = gene,
                size  = tf_gene_inverse_discordance_score,
                color = tf_gene_inverse_discordance_score
              )
            ) +
              ggplot2::geom_point(alpha = 0.9) +
              ggplot2::facet_wrap(~ tf, ncol = 3, scales = "free_y") +
              ggplot2::scale_color_gradient(
                low = "#c7e9c0",
                high = "darkgreen"
              ) +
              ggplot2::scale_size(range = c(3, 8)) +
              ggplot2::theme_classic(base_size = 12) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                panel.grid   = ggplot2::element_blank(),
                panel.spacing = unit(0.3, "lines"),
                strip.text   = ggplot2::element_text(face = "bold", size = 11),
                strip.background = ggplot2::element_rect(
                  fill = "grey95",
                  color = "grey70"
                )
              ) +
              ggplot2::labs(
                title = "Inverse Regulatory Discordance",
                y = "Target genes",
                color = "Inverse discordance score",
                size  = "Inverse discordance score"
              )
            
            res$inverse_path <- file.path(out_dir, "Figure_Inverse_Regulatory_Discordance.png")
            ggplot2::ggsave(
              res$inverse_path,
              res$p_inverse,
              width  = 12,
              height = esc_h,
              dpi    = 600
            )
            tf_sig_plot <- res$tf_discordance_significance
            top_n <- min(20, nrow(tf_sig_plot))
            plot_df <- tf_sig_plot[order(-tf_sig_plot$TF_discordance_score), ][1:top_n, , drop = FALSE]
            plot_df$significant <- plot_df$FDR < 0.05
            plot_df$TF <- factor(plot_df$TF, levels = rev(plot_df$TF))
            
            res$p_lollipop <- ggplot2::ggplot(
              plot_df,
              ggplot2::aes(x = TF, y = TF_discordance_score, color = significant)
            ) +
              ggplot2::geom_segment(
                ggplot2::aes(xend = TF, y = 0, yend = TF_discordance_score),
                linewidth = 1
              ) +
              ggplot2::geom_point(size = 4) +
              ggplot2::coord_flip() +
              ggplot2::scale_color_manual(
                values = c("FALSE" = "grey60", "TRUE" = "#D55E00"),
                name = "FDR < 0.05"
              ) +
              ggplot2::theme_classic(base_size = 14) +
              ggplot2::labs(
                title = "TF Regulatory Discordance Significance (Permutation + FDR)",
                x = "Transcription Factor",
                y = "TF Regulatory Discordance Score"
              )
            
            res$lollipop_path <- file.path(out_dir, "Figure_TF_Regulatory_Discordance_Lollipop.png")
            ggplot2::ggsave(res$lollipop_path, res$p_lollipop,
                            width = 7, height = 5, dpi = 600)
            
            shiny::incProgress(1, detail = "Done")
            message("DECODER-SC: analysis completed")
            add_status("✅ Analysis completed")
            res$done <- TRUE
          })
        
      }, error = function(e) {
        message("ERROR:", conditionMessage(e))
        add_status(paste("❌ ERROR:", conditionMessage(e)))
      })
    })
    
    shiny::reactive({ res })
  })
}
