#Ablation Analysis



#Ablation Analysis 1: No sign handling


run_decoder_no_sign <- function(obj) {
  
  library(Seurat)
  library(AUCell)
  library(dplyr)
  library(dorothea)
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^MT-")
  }
  
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  regulons_df <- dorothea_hs[
    dorothea_hs$confidence %in% c("A","B","C"),
  ]
  
  regulons <- split(regulons_df$target, regulons_df$tf)
  
  Idents(obj) <- obj$seurat_clusters
  
  avg_expr <- AverageExpression(
    obj,
    layer = "scale.data"
  )[[1]]
  
  rankings <- AUCell_buildRankings(
    avg_expr,
    plotStats = FALSE
  )
  
  tf_activity <- AUCell_calcAUC(
    regulons,
    rankings
  )
  
  tf_activity_mat <- as.matrix(
    getAUC(tf_activity)
  )
  
  common_genes <- intersect(
    rownames(avg_expr),
    regulons_df$target
  )
  
  state_expr <- avg_expr[common_genes, , drop = FALSE]
  
  rules <- regulons_df[
    regulons_df$target %in% common_genes,
  ]
  
  q <- 0.75
  
  tf_state <- tf_activity_mat >
    apply(tf_activity_mat, 1, quantile, probs = q)
  
  gene_state <- state_expr >
    apply(state_expr, 1, quantile, probs = q)
  
  gene_z <- t(scale(t(state_expr)))
  
  discordance_mat <- matrix(
    0,
    nrow = nrow(rules),
    ncol = ncol(state_expr)
  )
  
  for (i in seq_len(nrow(rules))) {
    
    tf <- rules$tf[i]
    tg <- rules$target[i]
    
    if (tf %in% rownames(tf_state)) {
      
      tf_on <- tf_state[tf, ]
      gene_on <- gene_state[tg, ]
      
      # NO SIGN HANDLING:
      # same rule for all TFs
      observed <- ifelse(
        tf_on & !gene_on,
        abs(gene_z[tg, ]),
        0
      )
      
      expected <- mean(tf_on) * mean(!gene_on)
      
      discordance_mat[i, ] <-
        pmax(observed - expected, 0)
    }
  }
  
  rules$discordance_score <- rowMeans(discordance_mat)
  
  tf_summary <- aggregate(
    discordance_score ~ tf,
    rules,
    mean
  )
  
  tf_summary <- tf_summary[
    order(-tf_summary$discordance_score),
  ]
  
  return(tf_summary)
}

no_sign_tf <- run_decoder_no_sign(sim_seu_perturbed)
no_sign_tf[no_sign_tf$tf %in% simulated_tfs, ]

perturbed_tf[perturbed_tf$tf %in% simulated_tfs, ]
ablation_df <- data.frame(
  
  tf = simulated_tfs,
  
  Full_DECODER_SC =
    perturbed_tf[
      match(simulated_tfs, perturbed_tf$tf),
      "discordance_score"
    ],
  
  No_Sign_Handling =
    no_sign_tf[
      match(simulated_tfs, no_sign_tf$tf),
      "discordance_score"
    ]
)
library(reshape2)

ablation_long <- melt(
  ablation_df,
  id.vars = "tf"
)
p_ablation <- ggplot(
  ablation_long,
  aes(
    x = tf,
    y = value,
    fill = variable
  )
) +
  
  geom_bar(
    stat = "identity",
    position = "dodge",
    width = 0.7
  ) +
  
  theme_classic(base_size = 16) +
  
  labs(
    title = "Ablation analysis of DECODER-SC",
    subtitle = "Sign-aware modeling improves TF-target rewiring detection",
    x = "Simulated transcription factors",
    y = "Regulatory discordance score",
    fill = ""
  ) +
  
  scale_fill_manual(
    values = c(
      "tomato",
      "skyblue3"
    )
  ) +
  
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

p_ablation
ggsave(
  "Figure_Ablation_Analysis.png",
  p_ablation,
  width = 8,
  height = 5,
  dpi = 1000,
  bg = "white"
)














setwd("D:/Interdisciplinary2222/New Method")

library(Seurat)

seu <- readRDS(
  "GSE185204_lung cancer_dataset.rds"
)

seu



set.seed(123)

cells_use <- sample(
  colnames(seu),
  20000
)

sim_seu <- subset(
  seu,
  cells = cells_use
)

sim_seu



DefaultAssay(sim_seu) <- "RNA"

sim_seu <- JoinLayers(sim_seu)

sim_seu






sim_seu <- NormalizeData(sim_seu)

sim_seu <- FindVariableFeatures(sim_seu)

sim_seu <- ScaleData(sim_seu)

sim_seu <- RunPCA(sim_seu)

sim_seu <- FindNeighbors(sim_seu, dims = 1:20)

sim_seu <- FindClusters(sim_seu, resolution = 0.5)



Idents(sim_seu) <- sim_seu$seurat_clusters

avg_expr <- AverageExpression(
  sim_seu,
  layer = "scale.data"
)[[1]]


library(dorothea)

regulons_df <- dorothea_hs[
  dorothea_hs$confidence %in% c("A","B","C"),
]


regulons <- split(
  regulons_df$target,
  regulons_df$tf
)
library(AUCell)

rankings <- AUCell_buildRankings(
  avg_expr,
  plotStats = FALSE
)

tf_activity <- AUCell_calcAUC(
  regulons,
  rankings
)

tf_activity_mat <- as.matrix(
  getAUC(tf_activity)
)



common_genes <- intersect(
  rownames(avg_expr),
  regulons_df$target
)

state_expr <- avg_expr[
  common_genes,
  ,
  drop = FALSE
]

rules <- regulons_df[
  regulons_df$target %in% common_genes,
]

rules$gene <- rules$target











run_decoder_q <- function(q_value) {
  
  tf_state <- tf_activity_mat >
    apply(tf_activity_mat, 1, quantile, probs = q_value)
  
  gene_state <- state_expr >
    apply(state_expr, 1, quantile, probs = q_value)
  
  gene_z <- t(scale(t(state_expr)))
  
  tf_discordance_mat <- matrix(
    0,
    nrow = nrow(rules),
    ncol = ncol(state_expr)
  )
  
  for(i in seq_len(nrow(rules))) {
    
    tf  <- rules$tf[i]
    tg  <- rules$target[i]
    mor <- rules$mor[i]
    
    if(tf %in% rownames(tf_state)) {
      
      tf_on   <- tf_state[tf, ]
      gene_on <- gene_state[tg, ]
      
      if(mor == 1) {
        
        observed <- ifelse(
          tf_on & !gene_on,
          abs(gene_z[tg, ]),
          0
        )
        
        expected <- mean(tf_on) * mean(!gene_on)
        
        tf_discordance_mat[i, ] <- pmax(
          observed - expected,
          0
        )
      }
      
      if(mor == -1) {
        
        observed <- ifelse(
          tf_on & gene_on,
          abs(gene_z[tg, ]),
          0
        )
        
        expected <- mean(tf_on) * mean(gene_on)
        
        tf_discordance_mat[i, ] <- pmax(
          observed - expected,
          0
        )
      }
    }
  }
  
  rules$discordance <- rowMeans(tf_discordance_mat)
  
  tf_rank <- aggregate(
    discordance ~ tf,
    rules,
    mean
  )
  
  tf_rank <- tf_rank[
    order(-tf_rank$discordance),
  ]
  
  return(tf_rank)
}

res_q05 <- run_decoder_q(0.5)

res_q06 <- run_decoder_q(0.6)

res_q07 <- run_decoder_q(0.7)

res_q08 <- run_decoder_q(0.8)

res_q09 <- run_decoder_q(0.9)


cor(
  res_q07$discordance,
  res_q05$discordance,
  method = "spearman"
)

cor(
  res_q07$discordance,
  res_q06$discordance,
  method = "spearman"
)

cor(
  res_q07$discordance,
  res_q08$discordance,
  method = "spearman"
)

cor(
  res_q07$discordance,
  res_q09$discordance,
  method = "spearman"
)




sens_df <- data.frame(
  Comparison = c(
    "0.7 vs 0.5",
    "0.7 vs 0.6",
    "0.7 vs 0.8",
    "0.7 vs 0.9"
  ),
  Spearman = c(
    0.9977493,
    1,
    1,
    0.9977493
  )
)

sens_df




library(ggplot2)

p <- ggplot(
  sens_df,
  aes(
    x = Comparison,
    y = Spearman,
    group = 1
  )
) +
  geom_point(
    size = 5,
    color = "#2C7FB8"
  ) +
  geom_line(
    linewidth = 1,
    color = "#2C7FB8"
  ) +
  geom_text(
    aes(label = round(Spearman, 3)),
    vjust = -1,
    size = 5
  ) +
  ylim(0.995, 1.001) +
  theme_classic(base_size = 15) +
  labs(
    title = "Threshold Sensitivity Analysis",
    x = "Quantile threshold comparison",
    y = "Spearman correlation"
  )

p

ggsave(
  filename = "Supplementary_Threshold_Sensitivity.png",
  plot = p,
  width = 6,
  height = 4.5,
  dpi = 1000,
  bg = "white"
)