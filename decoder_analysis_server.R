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
      p_escape = NULL,
      tf_failure_summary = NULL,
      tf_gene_failure_all = NULL,
      tf_gene_escape_all = NULL,
      heatmap_path = NULL,
      escape_path  = NULL,
      p_lollipop = NULL,
      tf_failure_significance = NULL,
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
            obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)
            obj <- NormalizeData(obj, verbose = FALSE)
            obj <- FindVariableFeatures(obj, verbose = FALSE)
            obj <- ScaleData(obj, verbose = FALSE)
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
            avg_expr <- Seurat::AverageExpression(obj, layer = "scale.data")[[1]]
            
            rankings <- AUCell::AUCell_buildRankings(avg_expr, plotStats = FALSE)
            tf_activity <- AUCell::AUCell_calcAUC(regulons, rankings)
            tf_activity_mat <- as.matrix(AUCell::getAUC(tf_activity))
            
            # ---------------------------
            # FAILURE / ESCAPE
            # ---------------------------
            shiny::incProgress(0.65, detail = "Computing regulatory failure and escape")
            message("Computing regulatory failure and escape")
            common_genes <- intersect(rownames(avg_expr), regulons_df$target)
            state_expr <- avg_expr[common_genes, , drop = FALSE]
            
            rules <- regulons_df[regulons_df$target %in% common_genes, ]
            rules$gene <- rules$target
            
            q <- input$tf_quantile
            tf_state   <- tf_activity_mat > apply(tf_activity_mat, 1, quantile, probs = q)
            gene_state <- state_expr > apply(state_expr, 1, quantile, probs = q)
            
            tf_failure_mat <- matrix(0, nrow = nrow(rules), ncol = ncol(state_expr))
            gene_escape_mat <- tf_failure_mat
            
            for (i in seq_len(nrow(rules))) {
              tf <- rules$tf[i]
              tg <- rules$target[i]
              if (tf %in% rownames(tf_state)) {
                tf_failure_mat[i, ] <- tf_state[tf, ] & !gene_state[tg, ]
                gene_escape_mat[i, ] <- !tf_state[tf, ] & gene_state[tg, ]
              }
            }
            
            rules$tf_gene_failure_freq <- rowMeans(tf_failure_mat)
            rules$tf_gene_escape_freq  <- rowMeans(gene_escape_mat)
            
            # ---------------------------
            # SAVE CSVs
            # ---------------------------
            shiny::incProgress(0.80, detail = "Saving result tables")
            message("Saving CSV outputs")
            write.csv(
              rules[, c("tf","confidence","gene","tf_gene_failure_freq")],
              file.path(out_dir, "TF_Gene_Failure_ALL.csv"),
              row.names = FALSE
            )
            
            write.csv(
              rules[, c("tf","confidence","gene","tf_gene_escape_freq")],
              file.path(out_dir, "TF_Gene_Escape_ALL.csv"),
              row.names = FALSE
            )
            
            # ---------------------------
            # TF RANKINGS
            # ---------------------------
            shiny::incProgress(0.85, detail = "Ranking transcription factors")
            message("Ranking transcription factors")
            tf_fail_rank <- aggregate(tf_gene_failure_freq ~ tf, rules, mean)
            tf_fail_rank <- tf_fail_rank[order(-tf_fail_rank$tf_gene_failure_freq), ]
            
            tf_escape_rank <- aggregate(tf_gene_escape_freq ~ tf, rules, mean)
            tf_escape_rank <- tf_escape_rank[order(-tf_escape_rank$tf_gene_escape_freq), ]
            
            res$tf_failure_summary <- tf_fail_rank
            
            write.csv(
              tf_fail_rank,
              file.path(out_dir, "TF_failure_frequency_all_TFs.csv"),
              row.names = FALSE
            )
            
            # =====================================================
            # LIMIT TFs FOR PLOTS ONLY
            # =====================================================
            plot_top_n_tf <- min(input$top_n_tf, 10)
            top_tfs_fail   <- head(tf_fail_rank$tf,   plot_top_n_tf)
            top_tfs_escape <- head(tf_escape_rank$tf, plot_top_n_tf)
            
            # ---------------------------
            # TOP 5 GENES PER TF
            # ---------------------------
            shiny::incProgress(0.90, detail = "Preparing plot data")
            message("Preparing plot data")
            fail_top <- do.call(
              rbind,
              lapply(split(
                rules[rules$tf %in% top_tfs_fail &
                        rules$tf_gene_failure_freq > 0, ],
                rules$tf
              ), function(df) head(df[order(-df$tf_gene_failure_freq), ], 5))
            )
            
            esc_top <- do.call(
              rbind,
              lapply(split(
                rules[rules$tf %in% top_tfs_escape &
                        rules$tf_gene_escape_freq > 0, ],
                rules$tf
              ), function(df) head(df[order(-df$tf_gene_escape_freq), ], 5))
            )
            
            fail_top <- fail_top[!is.na(fail_top$tf) & !is.na(fail_top$gene), ]
            esc_top  <- esc_top[!is.na(esc_top$tf)  & !is.na(esc_top$gene), ]
            
            fail_top$tf <- droplevels(factor(fail_top$tf))
            esc_top$tf  <- droplevels(factor(esc_top$tf))
            
            write.csv(
              fail_top,
              file.path(out_dir, "TF_Gene_Failure_TopTFs_TOP5.csv"),
              row.names = FALSE
            )
            
            write.csv(
              esc_top,
              file.path(out_dir, "TF_Gene_Escape_TopTFs_TOP5.csv"),
              row.names = FALSE
            )
            
            # ---------------------------
            # PERMUTATION + FDR
            # ---------------------------
            shiny::incProgress(0.93, detail = "Permutation testing")
            message("Running permutation testing")
            tf_fail_score <- tapply(rowMeans(tf_failure_mat), rules$tf, mean)
            tf_fail_score <- sort(tf_fail_score, decreasing = TRUE)
            
            set.seed(1)
            n_perm <- 2000
            
            tf_names <- names(tf_fail_score)
            perm_p <- setNames(rep(NA_real_, length(tf_names)), tf_names)
            
            for (tf in tf_names) {
              idx <- which(rules$tf == tf)
              if (length(idx) < 10) next
              observed <- mean(rowMeans(tf_failure_mat[idx, , drop = FALSE]))
              perm_vals <- replicate(n_perm, {
                samp <- sample(seq_len(nrow(tf_failure_mat)), length(idx), replace = FALSE)
                mean(rowMeans(tf_failure_mat[samp, , drop = FALSE]))
              })
              perm_p[tf] <- mean(perm_vals >= observed)
            }
            
            perm_p[is.na(perm_p)] <- 1
            perm_fdr <- p.adjust(perm_p, method = "BH")
            
            tf_sig <- data.frame(
              TF = tf_names,
              TF_failure_score = as.numeric(tf_fail_score),
              permutation_p = as.numeric(perm_p[tf_names]),
              FDR = as.numeric(perm_fdr[tf_names])
            )
            
            tf_sig <- tf_sig[order(tf_sig$FDR, -tf_sig$TF_failure_score), ]
            res$tf_failure_significance <- tf_sig
            
            write.csv(
              tf_sig,
              file.path(out_dir, "DECODER_SC_TF_failure_significance.csv"),
              row.names = FALSE
            )
            
            # ---------------------------
            # GENERATE PLOTS
            # ---------------------------
            shiny::incProgress(0.97, detail = "Generating plots")
            message("Generating plots")
            
            # ==========================================================
            # REFINED FAILURE HEATMAP DESIGN (TF on X, Gene on Y)
            # ==========================================================
            # Order factors for clean display
            fail_top_plot <- fail_top |>
              dplyr::group_by(tf) |>
              dplyr::slice_max(tf_gene_failure_freq, n = 5, with_ties = FALSE) |>
              dplyr::ungroup()
            
            fail_top_plot$tf <- factor(
              fail_top_plot$tf,
              levels = unique(fail_top_plot$tf)
            )
            
            fail_top_plot$gene <- factor(
              fail_top_plot$gene,
              levels = rev(unique(fail_top_plot$gene))
            )
            
            # ---------------------------
            # HEATMAP (FACETED, SAFE SIZE)
            # ---------------------------
            n_tf <- length(unique(fail_top_plot$tf))
            heat_h <- max(6, min(14, 3.5 * n_tf))
            
            res$p_heat <- ggplot2::ggplot(
              fail_top_plot,
              ggplot2::aes(x = 1, y = gene, fill = tf_gene_failure_freq)
            ) +
              ggplot2::geom_tile(color = "white") +
              ggplot2::facet_wrap(~ tf, ncol = 3, scales = "free_y") +
              ggplot2::scale_fill_gradient(low = "white", high = "red") +
              ggplot2::theme_minimal(base_size = 12) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_blank(),
                panel.grid   = ggplot2::element_blank(),
                strip.text   = ggplot2::element_text(face = "bold", size = 11)
              ) +
              ggplot2::labs(
                title = "Regulatory Execution Failure (TF ON → Gene OFF)",
                y = "Target genes"
              )
            
            res$heatmap_path <- file.path(out_dir, "Figure_Failure_Faceted_TF_Gene.png")
            ggplot2::ggsave(
              res$heatmap_path,
              res$p_heat,
              width  = 12,
              height = heat_h,
              dpi    = 600
            )
            
            # ==========================================================
            # REFINED ESCAPE DOTPLOT DESIGN (TF on X, Gene on Y)
            # ==========================================================
            esc_top_plot <- esc_top |>
              dplyr::group_by(tf) |>
              dplyr::slice_max(tf_gene_escape_freq, n = 5, with_ties = FALSE) |>
              dplyr::ungroup()
            
            esc_top_plot$tf <- factor(
              esc_top_plot$tf,
              levels = unique(esc_top_plot$tf)
            )
            
            esc_top_plot$gene <- factor(
              esc_top_plot$gene,
              levels = rev(unique(esc_top_plot$gene))
            )
            
            # ---------------------------
            # ESCAPE DOT PLOT (FACETED, SAFE SIZE)
            # ---------------------------
            n_tf <- length(unique(esc_top_plot$tf))
            esc_h <- max(6, min(14, 3.5 * n_tf))
            
            res$p_escape <- ggplot2::ggplot(
              esc_top_plot,
              ggplot2::aes(
                x = 1,
                y = gene,
                size  = tf_gene_escape_freq,
                color = tf_gene_escape_freq
              )
            ) +
              ggplot2::geom_point(alpha = 0.85) +
              ggplot2::facet_wrap(~ tf, ncol = 3, scales = "free_y") +
              ggplot2::scale_color_gradient(low = "lightyellow", high = "darkgreen") +
              ggplot2::theme_classic(base_size = 12) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_blank(),
                panel.grid   = ggplot2::element_blank(),
                strip.text   = ggplot2::element_text(face = "bold", size = 11)
              ) +
              ggplot2::labs(
                title = "Regulatory Escape (TF OFF → Gene ON)",
                y = "Target genes"
              )
            
            res$escape_path <- file.path(out_dir, "Figure_Escape_Faceted_TF_Gene.png")
            ggplot2::ggsave(
              res$escape_path,
              res$p_escape,
              width  = 12,
              height = esc_h,
              dpi    = 600
            )
            tf_sig_plot <- res$tf_failure_significance
            top_n <- min(20, nrow(tf_sig_plot))
            plot_df <- tf_sig_plot[order(-tf_sig_plot$TF_failure_score), ][1:top_n, , drop = FALSE]
            plot_df$significant <- plot_df$FDR < 0.05
            plot_df$TF <- factor(plot_df$TF, levels = rev(plot_df$TF))
            
            res$p_lollipop <- ggplot2::ggplot(
              plot_df,
              ggplot2::aes(x = TF, y = TF_failure_score, color = significant)
            ) +
              ggplot2::geom_segment(
                ggplot2::aes(xend = TF, y = 0, yend = TF_failure_score),
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
                title = "TFs with Regulatory Execution Failure (Permutation + FDR)",
                x = "Transcription Factor",
                y = "TF Failure Score"
              )
            
            res$lollipop_path <- file.path(out_dir, "Figure_TF_Failure_Lollipop.png")
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
