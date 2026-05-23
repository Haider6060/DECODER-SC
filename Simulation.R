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


simulated_tfs <- c(
  "TBX21",
  "BACH1",
  "MAFB",
  "RUNX3",
  "PPARA"
)

simulated_tfs



expr_mat <- GetAssayData(
  sim_seu,
  layer = "counts"
)

dim(expr_mat)
sim_expr <- expr_mat
set.seed(123)

for (tf in simulated_tfs) {
  
  tf_targets <- dorothea_hs %>%
    dplyr::filter(tf == !!tf) %>%
    dplyr::pull(target)
  
  tf_targets <- intersect(
    tf_targets,
    rownames(sim_expr)
  )
  
  if (length(tf_targets) < 10) next
  
  selected_targets <- sample(
    tf_targets,
    min(20, length(tf_targets))
  )
  
  selected_cells <- sample(
    colnames(sim_expr),
    floor(0.3 * ncol(sim_expr))
  )
  
  sim_expr[
    selected_targets,
    selected_cells
  ] <- 0
  
  cat(
    "Simulated discordance for:",
    tf,
    "\n"
  )
}



sim_seu_perturbed <- CreateSeuratObject(
  counts = sim_expr
)

sim_seu_perturbed

original_tf <- run_decoder_validation(sim_seu)
head(original_tf)
perturbed_tf <- run_decoder_validation(sim_seu_perturbed)
head(perturbed_tf)

original_tf[original_tf$tf %in% simulated_tfs, ]

perturbed_tf[perturbed_tf$tf %in% simulated_tfs, ]




plot_df$improvement <-
  plot_df$perturbed - plot_df$original

plot_df <- plot_df[
  order(plot_df$improvement),
]

p_rank <- ggplot(
  plot_df,
  aes(
    x = improvement,
    y = reorder(tf, improvement)
  )
) +
  
  geom_segment(
    aes(
      x = 0,
      xend = improvement,
      y = tf,
      yend = tf
    ),
    linewidth = 2,
    color = "skyblue3"
  ) +
  
  geom_point(
    size = 7,
    color = "tomato"
  ) +
  
  theme_classic(base_size = 16) +
  
  labs(
    title = "Rank improvement of simulated transcription factors",
    subtitle = "Selective recovery of simulated TF-target regulatory rewiring by DECODER-SC",
    x = "Regulatory discordance increase",
    y = "Simulated TFs"
  ) +
  
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

p_rank


ggsave(
  "Figure_Simulation_PositiveControl_RankImprovement.png",
  p_rank,
  width = 8,
  height = 5,
  dpi = 1000,
  bg = "white"
)





#Simulation 2 Random Perturbation Negative Control


random_expr <- expr_mat

set.seed(456)

random_genes <- sample(
  rownames(random_expr),
  100
)

random_cells <- sample(
  colnames(random_expr),
  floor(0.3 * ncol(random_expr))
)

random_expr[
  random_genes,
  random_cells
] <- 0

cat(
  "Random perturbation completed\n"
)


random_seu <- CreateSeuratObject(
  counts = random_expr
)

random_seu


random_tf <- run_decoder_validation(
  random_seu
)
head(random_tf)

mean(random_tf$discordance_score)

mean(original_tf$discordance_score)







#Simulation 3 nOisy



set.seed(789)

noise_tf <- original_tf

noise_tf$discordance_score <-
  noise_tf$discordance_score +
  rnorm(
    nrow(noise_tf),
    mean = 0,
    sd = 0.01
  )

noise_tf$discordance_score[
  noise_tf$discordance_score < 0
] <- 0

cor(
  original_tf$discordance_score,
  noise_tf$discordance_score,
  method = "spearman"
)

mean(original_tf$discordance_score)
mean(noise_tf$discordance_score)
library(ggplot2)
library(patchwork)
plot_df$improvement <-
  plot_df$perturbed - plot_df$original

plot_df <- plot_df[
  order(plot_df$improvement),
]

p1 <- ggplot(
  plot_df,
  aes(
    x = improvement,
    y = reorder(tf, improvement)
  )
) +
  
  geom_segment(
    aes(
      x = 0,
      xend = improvement,
      y = tf,
      yend = tf
    ),
    linewidth = 2,
    color = "skyblue3"
  ) +
  
  geom_point(
    size = 6,
    color = "tomato"
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "A. Positive-control rewiring",
    x = "Regulatory discordance increase",
    y = "Simulated TFs"
  ) +
  
  theme(
    plot.title = element_text(face = "bold")
  )



spec_df <- data.frame(
  group = c("Original", "Random perturbation"),
  score = c(
    mean(original_tf$discordance_score),
    mean(random_tf$discordance_score)
  )
)

p2 <- ggplot(
  spec_df,
  aes(
    x = group,
    y = score,
    fill = group
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.6
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "B. Random perturbation control",
    x = "",
    y = "Mean discordance"
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

rho <- cor(
  original_tf$discordance_score,
  noise_tf$discordance_score,
  method = "spearman"
)

robust_df <- data.frame(
  original = original_tf$discordance_score,
  noisy = noise_tf$discordance_score
)

p3 <- ggplot(
  robust_df,
  aes(
    x = original,
    y = noisy
  )
) +
  
  geom_point(
    alpha = 0.6,
    color = "darkgreen",
    size = 2
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black"
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = paste0(
      "C. Noise robustness (ρ = ",
      round(rho, 3),
      ")"
    ),
    x = "Original discordance",
    y = "Noisy discordance"
  ) +
  
  theme(
    plot.title = element_text(face = "bold")
  )


combined_simulation <- 
  (p1 | p2 | p3)

combined_simulation


ggsave(
  "Figure_Combined_Simulation_Validation.png",
  combined_simulation,
  width = 16,
  height = 5,
  dpi = 1000,
  bg = "white"
)
