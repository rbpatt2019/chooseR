source("R/pipeline.R")
library(Seurat)
library(ggplot2)
`%>%` <- magrittr::`%>%`

# Load Seurat object
# Make sure this has already been normalised and that PCA has been performed
# The data included here is the Ding, et al., Nature Biotechnology, 2020
# human PBMC Smart-Seq data from our manuscript
obj <- readRDS("data/ding_smartseq_pbmc_preprocessed.rds")

# Define the number of PCs to use, and which assay and reduction to use.
# We recommend testing a broad range of resolutions
# For more on picking the correct number of PCs, see:
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
npcs <- 100
resolutions <- c(0.8, 1, 1.2, 1.6, 2, 4, 6, 8, 12, 16)
assay <- "SCT"
reduction <- "pca"
results_path <- "results/"

# Run pipeline
for (res in resolutions) {
  message(paste0("Clustering ", res, "..."))
  message("\tFinding ground truth...")

  # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
  obj <- find_clusters(
                       obj,
                       reduction = reduction,
                       assay = assay,
                       resolution = res
  )
  clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]

  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(
                              obj,
                              n = 100,
                              size = 0.8,
                              npcs = npcs,
                              res = res,
                              reduction = reduction,
                              assay = assay
  )

  # Now calculate the co-clustering frequencies
  message(paste0("Tallying ", res, "..."))
  # This is the more time efficient vectorisation
  # However, it exhausts vector memory for (nearly) all datasets
  # matches <- purrr::map(columns, find_matches, df = results)
  # matches <- purrr::reduce(matches, `+`)
  columns <- colnames(dplyr::select(results, -cell))
  mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
  i <- 1 # Counter
  for (col in columns) {
    message(paste0("\tRound ", i, "..."))
    mtchs <- Reduce("+", list(
      mtchs,
      find_matches(col, df = results)
    ))
    i <- i + 1
  }

  message(paste0("Scoring ", res, "..."))
  mtchs <- dplyr::mutate_all(
    dplyr::as_tibble(mtchs),
    function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
  )

  # Now calculate silhouette scores
  message(paste0("Silhouette ", res, "..."))
  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(clusters))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))

  # Finally, calculate grouped metrics
  message(paste0("Grouping ", res, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
  sil <- group_sil(sil, res)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
}

# Save original data, with ground truth labels
saveRDS(obj, paste0(results_path, "clustered_data.rds"))

# Create silhouette plot
# Read in scores and calculate CIs
scores <- purrr::map(
  paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
  readRDS
)
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()
meds <- scores %>%
  dplyr::group_by(res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)

writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))

# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
  meds %>%
  dplyr::filter(med >= threshold) %>%
  dplyr::arrange(n_clusters) %>%
  tail(n = 1) %>%
  dplyr::pull(res)
)

# And plot!
ggplot(meds, aes(factor(res), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = "grey",
    size = 0.25
  ) +
  geom_hline(aes(yintercept = threshold), colour = "blue") +
  geom_vline(aes(xintercept = choice), colour = "red") +
  geom_jitter(
    data = scores,
    aes(factor(res), avg_sil),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = paste0(results_path, "silhouette_distribution_plot.png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

# Finally, a dot plot of silhouette scores to help identify less robust clusters
# The initial pipe is to order the clusters by silhouette score
scores %>%
  dplyr::filter(res == choice) %>%
  dplyr::arrange(dplyr::desc(avg_sil)) %>%
  dplyr::mutate_at("cluster", ordered, levels = .$cluster) %>%
  ggplot(aes(factor(cluster), avg_sil)) +
    geom_point() +
    scale_x_discrete("Cluster") +
    scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
    ) +
    cowplot::theme_minimal_grid() +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
    )

ggsave(
  filename = paste0(results_path, "silhouette_point_plot_", choice, ".png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)
