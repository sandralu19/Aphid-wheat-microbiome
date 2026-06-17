###############################################################################
# Procrustes analyses linking chemical and microbial datasets
# Manuscript: Aphid herbivory transiently restructures rhizosphere chemistry
#             and bacterial communities in wheat
# Author: Sandra Cortes-Patiño
# Tidied for GitHub/reproducibility
###############################################################################

# This script documents the Procrustes analyses used for:
#   - Main Fig. 4: leaf VOC profiles vs rhizosphere microbiome composition
#   - Supplementary Fig. 4: bootstrap distribution for Fig. 4
#   - Supplementary Fig. 7: leaf VOCs vs rhizosphere VOCs/non-volatile metabolites

# TODO: Confirm that input ordination score files are the final versions used in the manuscript.
# TODO: If ordination scores are generated elsewhere, point users to those scripts in README.

source("scripts/00_project_setup.R")

packages <- c("vegan", "ggplot2", "dplyr", "patchwork")
invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. Input/output paths
# -----------------------------------------------------------------------------

procrustes_dir <- file.path(DATA_DIR, "procrustes")
procrustes_out <- file.path(OUTPUT_DIR, "procrustes")
procrustes_fig <- file.path(FIGURE_DIR, "procrustes")

dir.create(procrustes_out, recursive = TRUE, showWarnings = FALSE)
dir.create(procrustes_fig, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

read_scores <- function(filename) {
  read.csv(file.path(procrustes_dir, filename), row.names = 1, check.names = FALSE)
}

match_rows <- function(target, rotated) {
  # Match rows of rotated dataset to target dataset.
  common_samples <- intersect(rownames(target), rownames(rotated))
  if (length(common_samples) < 3) {
    stop("Fewer than three shared samples between ordination matrices.")
  }
  list(
    target = target[common_samples, , drop = FALSE],
    rotated = rotated[common_samples, , drop = FALSE]
  )
}

run_procrustes_pair <- function(target_scores, rotated_scores, permutations = 9999) {
  matched <- match_rows(target_scores, rotated_scores)
  proc <- vegan::procrustes(
    X = as.matrix(matched$target),
    Y = as.matrix(matched$rotated)
  )
  protest_res <- vegan::protest(
    X = as.matrix(matched$target),
    Y = as.matrix(matched$rotated),
    permutations = permutations
  )

  list(
    procrustes = proc,
    protest = protest_res,
    target = matched$target,
    rotated = matched$rotated
  )
}

extract_procrustes_plot_df <- function(proc_result, metadata, target_label, rotated_label, group_col = "Treatment") {
  proc <- proc_result$procrustes

  df_segments <- data.frame(
    SampleID = rownames(proc$X),
    X1 = proc$X[, 1],
    X2 = proc$X[, 2],
    Y1 = proc$Yrot[, 1],
    Y2 = proc$Yrot[, 2]
  )

  metadata$SampleID <- rownames(metadata)
  df_segments <- dplyr::left_join(df_segments, metadata, by = "SampleID")

  df_points <- dplyr::bind_rows(
    data.frame(
      SampleID = rownames(proc$X),
      Dim1 = proc$X[, 1],
      Dim2 = proc$X[, 2],
      Dataset = target_label
    ),
    data.frame(
      SampleID = rownames(proc$Yrot),
      Dim1 = proc$Yrot[, 1],
      Dim2 = proc$Yrot[, 2],
      Dataset = rotated_label
    )
  ) |>
    dplyr::left_join(metadata, by = "SampleID")

  list(segments = df_segments, points = df_points)
}

plot_procrustes_pair <- function(proc_result,
                                 metadata,
                                 target_label,
                                 rotated_label,
                                 group_col = "Treatment",
                                 label_text = NULL,
                                 show_legend = TRUE) {
  plot_data <- extract_procrustes_plot_df(proc_result, metadata, target_label, rotated_label, group_col)

  p <- ggplot() +
    geom_segment(
      data = plot_data$segments,
      aes(x = X1, y = X2, xend = Y1, yend = Y2, colour = .data[[group_col]]),
      arrow = arrow(length = unit(0.16, "cm")),
      alpha = 0.55,
      linewidth = 0.6
    ) +
    geom_point(
      data = plot_data$points,
      aes(x = Dim1, y = Dim2, colour = .data[[group_col]], shape = Dataset),
      size = 3.2,
      stroke = 1
    ) +
    scale_colour_manual(values = condition_colours, na.value = "grey60") +
    coord_equal() +
    labs(x = "Dimension 1", y = "Dimension 2", colour = "Group", shape = "Dataset") +
    manuscript_theme(base_size = 12) +
    guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2))

  p <- p + scale_shape_manual(
    values = setNames(c(16, 17), c(target_label, rotated_label))
  )

  if (!is.null(label_text)) {
    p <- p + annotate("text", x = Inf, y = -Inf, label = label_text, hjust = 1.1, vjust = -0.8, size = 4)
  }

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}

bootstrap_procrustes_r <- function(target_scores, rotated_scores, n_boot = 10000, seed = 123) {
  matched <- match_rows(target_scores, rotated_scores)
  n_samples <- nrow(matched$target)

  set.seed(seed)
  boot_r <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n_samples), replace = TRUE)
    proc_sub <- vegan::procrustes(
      X = as.matrix(matched$target[idx, , drop = FALSE]),
      Y = as.matrix(matched$rotated[idx, , drop = FALSE])
    )
    boot_r[i] <- cor(as.vector(proc_sub$X), as.vector(proc_sub$Yrot), method = "pearson")
  }

  data.frame(r = boot_r)
}

plot_bootstrap_histogram <- function(boot_df) {
  mean_r <- mean(boot_df$r, na.rm = TRUE)
  ci <- quantile(boot_df$r, probs = c(0.025, 0.975), na.rm = TRUE)

  ggplot(boot_df, aes(x = r)) +
    geom_histogram(bins = 50, fill = "skyblue", colour = "white") +
    geom_vline(xintercept = mean_r, colour = "red", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = ci, colour = "darkgreen", linetype = "dotted", linewidth = 0.8) +
    labs(
      title = "Bootstrap Procrustes correlation",
      x = "Correlation r",
      y = "Frequency"
    ) +
    manuscript_theme(base_size = 12) +
    theme(legend.position = "none")
}

save_pdf <- function(plot, filename, width = 16, height = 10) {
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, units = "cm", device = cairo_pdf)
}

# -----------------------------------------------------------------------------
# 3. Load ordination scores and metadata
# -----------------------------------------------------------------------------

ord_micro <- read_scores("pcoa_distances_procrusts.csv")
ord_leaf_vocs_all <- read_scores("pca_distances_procrusts_above.csv")
ord_ecoplates <- read_scores("pca_distances_procrusts_ecoplates.csv")

ord_leaf_vocs_2w <- read_scores("pca_distances_procrusts_above2weeks.csv")
ord_rhizo_vocs_2w <- read_scores("pca_distances_procrusts_below_2weeks.csv")
ord_nonvolatiles_2w <- read_scores("pca_distances_procrusts_metabolites.csv")

metadata_all <- read.csv(file.path(procrustes_dir, "meta_proc_all.csv"), row.names = 1, check.names = FALSE)
metadata_2w <- read.csv(file.path(procrustes_dir, "meta_proc_2weeks.csv"), row.names = 1, check.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Main Fig. 4: leaf VOCs vs rhizosphere microbiome
# -----------------------------------------------------------------------------

leaf_micro <- run_procrustes_pair(ord_leaf_vocs_all, ord_micro, permutations = 9999)
leaf_micro_summary <- data.frame(
  comparison = "Leaf VOCs vs rhizosphere microbiome",
  r = leaf_micro$protest$t0,
  p_value = leaf_micro$protest$signif,
  permutations = leaf_micro$protest$permutations
)
write.csv(leaf_micro_summary, file.path(procrustes_out, "procrustes_leaf_vocs_microbiome.csv"), row.names = FALSE)

fig4 <- plot_procrustes_pair(
  proc_result = leaf_micro,
  metadata = metadata_all,
  target_label = "Leaf VOCs",
  rotated_label = "Rhizosphere microbiome",
  group_col = "Treatment",
  label_text = paste0("r = ", round(leaf_micro$protest$t0, 2), ", p = ", signif(leaf_micro$protest$signif, 3))
)

save_pdf(fig4, file.path(procrustes_fig, "figure4_leaf_vocs_microbiome.pdf"), width = 18, height = 12)

# Bootstrap distribution for Supplementary Fig. 4.
boot_leaf_micro <- bootstrap_procrustes_r(ord_leaf_vocs_all, ord_micro, n_boot = 10000, seed = 123)
write.csv(boot_leaf_micro, file.path(procrustes_out, "bootstrap_leaf_vocs_microbiome.csv"), row.names = FALSE)

supp_boot <- plot_bootstrap_histogram(boot_leaf_micro)
save_pdf(supp_boot, file.path(procrustes_fig, "supplementary_bootstrap_leaf_vocs_microbiome.pdf"), width = 14, height = 10)

# -----------------------------------------------------------------------------
# 5. Supplementary Fig. 7: leaf VOCs vs rhizosphere VOCs/non-volatiles
# -----------------------------------------------------------------------------

leaf_rhizo_vocs <- run_procrustes_pair(ord_leaf_vocs_2w, ord_rhizo_vocs_2w, permutations = 9999)
leaf_nonvolatiles <- run_procrustes_pair(ord_leaf_vocs_2w, ord_nonvolatiles_2w, permutations = 9999)

supp_summary <- data.frame(
  comparison = c("Leaf VOCs vs rhizosphere VOCs", "Leaf VOCs vs rhizosphere non-volatiles"),
  r = c(leaf_rhizo_vocs$protest$t0, leaf_nonvolatiles$protest$t0),
  p_value = c(leaf_rhizo_vocs$protest$signif, leaf_nonvolatiles$protest$signif),
  permutations = c(leaf_rhizo_vocs$protest$permutations, leaf_nonvolatiles$protest$permutations)
)
write.csv(supp_summary, file.path(procrustes_out, "procrustes_leaf_vocs_rhizosphere_metabolites.csv"), row.names = FALSE)

p_leaf_rhizo_vocs <- plot_procrustes_pair(
  proc_result = leaf_rhizo_vocs,
  metadata = metadata_2w,
  target_label = "Leaf VOCs",
  rotated_label = "Rhizosphere VOCs",
  group_col = "Treatment",
  label_text = paste0("r = ", round(leaf_rhizo_vocs$protest$t0, 2), ", p = ", signif(leaf_rhizo_vocs$protest$signif, 3)),
  show_legend = TRUE
)

p_leaf_nonvolatiles <- plot_procrustes_pair(
  proc_result = leaf_nonvolatiles,
  metadata = metadata_2w,
  target_label = "Leaf VOCs",
  rotated_label = "Rhizosphere non-volatiles",
  group_col = "Treatment",
  label_text = paste0("r = ", round(leaf_nonvolatiles$protest$t0, 2), ", p = ", signif(leaf_nonvolatiles$protest$signif, 3)),
  show_legend = TRUE
)

supp7 <- p_leaf_rhizo_vocs + p_leaf_nonvolatiles + patchwork::plot_annotation(tag_levels = "a")
save_pdf(supp7, file.path(procrustes_fig, "supplementary_leaf_vocs_rhizosphere_metabolites.pdf"), width = 22, height = 10)

# -----------------------------------------------------------------------------
# 6. Optional exploratory comparisons
# -----------------------------------------------------------------------------

# These comparisons were part of exploratory analysis and can be uncommented if needed.
# eco_micro <- run_procrustes_pair(ord_micro, ord_ecoplates, permutations = 9999)
# leaf_eco <- run_procrustes_pair(ord_leaf_vocs_all, ord_ecoplates, permutations = 9999)

# -----------------------------------------------------------------------------
# 7. Session information
# -----------------------------------------------------------------------------

write_session_info(file.path(procrustes_out, "session_info_procrustes.txt"))
