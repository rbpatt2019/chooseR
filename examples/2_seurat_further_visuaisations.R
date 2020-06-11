library(Seurat)
library(tidyverse)

# Currently, running ggsave with Rscript produces spurious pdfs.
# To prevent this, all plots are caught and saved directly
# Call print(plot) at any point to see your plot!

# Define sommon common variables
# choice is the res elected by the pipeline in examples/1_seurat_pipeline.R
# Be sure to change your path as necessary!
reduction <- "pca"
assay <- "SCT"
choice <- 1.6
results_path <- "results/"

# Load in the object containing the clustered results
obj <- readRDS(paste0(results_path, "clustered_data.rds"))

# First is a cluster average co-clustering heatmap
# Read the data
grp <- readRDS(paste0(results_path, "frequency_grouped_", choice, ".rds"))

# As the data is symmetrical, we do not need the upper triangle
grp <- grp %>%
  pivot_wider(names_from = "cell_2", values_from = "avg_percent") %>%
  select(str_sort(colnames(.), numeric = T)) %>%
  column_to_rownames("cell_1")
grp[lower.tri(grp)] <- NA
grp <- grp %>%
  as_tibble(rownames = "cell_1") %>%
  pivot_longer(-cell_1, names_to = "cell_2", values_to = "avg_percent") %>%
  mutate_at("cell_2", ordered, levels = unique(.$cell_1)) %>%
  mutate_at("cell_1", ordered, levels = unique(.$cell_1))

# And plot!
plot <- ggplot(grp, aes(factor(cell_1), cell_2, fill = avg_percent)) +
  geom_tile() +
  scale_x_discrete("Cluster", expand = c(0, 0)) +
  scale_y_discrete(
    "Cluster",
    limits = rev(levels(grp$cell_2)),
    expand = c(0, 0)
  ) +
  scale_fill_distiller(
    " ",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    palette = "RdYlBu",
    na.value = "white"
  ) +
  coord_fixed() +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = c(0.9, 0.9)
  ) +
  guides(fill = guide_colorbar(barheight = 3, barwidth = 1))

ggsave(
  plot = plot,
  filename = paste0(results_path, "coclustering_heatmap_", choice, ".png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

# Let's add the silhouette scores to the Seurat object!
sil_scores <- readRDS(paste0(results_path, "silhouette_", choice, ".rds"))
sil_scores <- as.data.frame(sil_scores[, 3], row.names = Seurat::Cells(obj))
colnames(sil_scores) <- c("sil_score")
obj <- AddMetaData(obj, metadata = sil_scores)

# Let's visualise the selected cluster
# If your data has known  clusters, you could also visualise those!
# Remember, truths are in "glue({reduction}.{assay}_res.{choice})"
# Seurat Changes color scheme if you order your data, so we provide
# the following helper function to restore defaults
gg_color <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  colours <- hcl(h = hues, c = 100, l = 65)[1:n]
  return(colours)
}

plot <- DimPlot(
  obj,
  reduction = "umap",
  group.by = glue::glue("{reduction}.{assay}_res.{choice}"),
  pt.size = 0.5,
  # cols = gg_color(6) # Only necessary if you have ordered your clusters
)

ggsave(
  plot = plot,
  filename = paste0(results_path, choice, "_cluster_umap.png"),
  dpi = 300,
  height = 5,
  width = 5,
  units = "in"
)

# We also find it useful to visualise the silhouette scores on the UMAP!
plot <- FeaturePlot(
  obj,
  "sil_score",
  reduction = "umap",
  pt.size = 0.5,
  min.cutoff = -1,
  max.cutoff = 1
) +
  scale_colour_distiller(
    palette = "RdYlBu",
    labels = c(-1, 0, 1),
    breaks = c(-1, 0, 1),
    limits = c(-1, 1)
  )

ggsave(
  plot = plot,
  filename = paste0(results_path, choice, "_silhouette_umap.png"),
  dpi = 300,
  height = 5,
  width = 5,
  units = "in"
)

# This dataset does have annotated cell types, so we can create a heatmap
# to see how well they line up. Remember, "truths" will always be in
# glue::glue("{reduction}.{assay}_res.{choice}"), but annotated cell types
# might be in different locations. Be sure to change the names as needed.
ids <- as_tibble(
    obj[[c(glue::glue("{reduction}.{assay}_res.{choice}"), "CellType")]]
  ) %>%
  mutate_at(
    c(glue::glue("{reduction}.{assay}_res.{choice}"), "CellType"),
    as.factor
  ) %>%
  group_by(pca.SCT_res.1.6, CellType) %>%
  summarise("count" = n()) %>%
  ungroup() %>%
  complete(pca.SCT_res.1.6, CellType, fill = list("count" = 0)) %>%
  # Count suggested clusters
  group_by(pca.SCT_res.1.6) %>%
  mutate("n_suggested" = sum(count)) %>%
  ungroup() %>%
  # Count known clusters
  group_by(CellType) %>%
  mutate("n_known" = sum(count)) %>%
  ungroup() %>%
  # Calculate statistics
  mutate("dice" = (2 * count) / (n_suggested + n_known)) %>%
  mutate(
    "pca.SCT_res.1.6" = fct_reorder(pca.SCT_res.1.6, dice, max, .desc = TRUE),
    "CellType" = fct_reorder(CellType, dice, max, .desc = TRUE)
  )

plot <- ggplot(ids, aes(pca.SCT_res.1.6, CellType, fill = dice)) +
  geom_tile() +
  scale_x_discrete("Suggested Clusters", expand = c(0, 0)) +
  scale_y_discrete(
    "Known Clusters",
    limits = rev(levels(ids$CellType)),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  scale_fill_distiller(
    NULL,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    palette = "RdYlBu"
  ) +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 7),
  ) +
  guides(fill = guide_colorbar(barheight = 3, barwidth = 0.5))

ggsave(
  plot = plot,
  filename = paste0(results_path, "cluster_relation_heatmap.png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)
