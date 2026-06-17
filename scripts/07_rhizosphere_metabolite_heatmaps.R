###############################################################################
# Rhizosphere VOC and non-volatile metabolite heatmaps
# Manuscript: Aphid herbivory transiently restructures rhizosphere chemistry
#             and bacterial communities in wheat
# Author: Sandra Cortes-Patiño
# Purpose: Generate supplementary heatmaps for rhizosphere VOCs and LC-MS
#          non-volatile metabolite features.
###############################################################################

source("scripts/00_project_setup.R")

library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -----------------------------------------------------------------------------
# 1. Input files
# -----------------------------------------------------------------------------

# Expected files:
# data/rhizosphere_metabolomics/data_norm_visual.csv
#   Z-score or normalised LC-MS feature table for Supplementary Figure 6.
#   First column should contain feature IDs; remaining columns are samples.
#
# data/rhizosphere_vocs/data_norm_VOCs_heat.csv
#   Z-score or normalised rhizosphere VOC feature table for Supplementary Figure 5.
#   First column should contain feature IDs; remaining columns are samples.
#
# TODO: Confirm whether these matrices are already Z-score normalised. If not,
# apply scaling before plotting and update the figure legend accordingly.

lcms_file <- file.path(DATA_DIR, "rhizosphere_metabolomics", "data_norm_visual.csv")
voc_file <- file.path(DATA_DIR, "rhizosphere_vocs", "data_norm_VOCs_heat.csv")

# -----------------------------------------------------------------------------
# 2. Shared heatmap function
# -----------------------------------------------------------------------------

heatmap_colours <- colorRamp2(
  c(-2, 0, 2),
  c("white", "#FFE1D4", "#F53D2A")
)

condition_colours_heatmap <- c(
  "Herbivory" = "#fc8d59",
  "NoHerb" = "#91bfdb",
  "NoH" = "#91bfdb"
)

read_feature_matrix <- function(path) {
  mat <- read.csv(path, check.names = FALSE) %>%
    as.data.frame()

  rownames(mat) <- mat[[1]]
  mat <- mat[, -1, drop = FALSE]
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat
}

make_condition_annotation <- function(conditions) {
  HeatmapAnnotation(
    Condition = conditions,
    col = list(Condition = condition_colours_heatmap)
  )
}

save_heatmap <- function(mat, conditions, output_file, width = 7, height = 9) {
  if (ncol(mat) != length(conditions)) {
    stop("Number of condition labels does not match number of matrix columns.")
  }

  ha <- make_condition_annotation(conditions)

  ht <- Heatmap(
    mat,
    name = "Z-score",
    col = heatmap_colours,
    top_annotation = ha,
    show_row_dend = TRUE,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE
  )

  pdf(output_file, width = width, height = height)
  draw(ht)
  dev.off()

  invisible(ht)
}

# -----------------------------------------------------------------------------
# 3. Supplementary Figure 5: rhizosphere VOC heatmap
# -----------------------------------------------------------------------------

voc_matrix <- read_feature_matrix(voc_file)

# Original exploratory script used five Herbivory and four NoH labels.
# TODO: Confirm final sample order in `data_norm_VOCs_heat.csv`.
voc_conditions <- c(
  rep("Herbivory", 5),
  rep("NoHerb", 4)
)

save_heatmap(
  mat = voc_matrix,
  conditions = voc_conditions,
  output_file = file.path(FIGURE_DIR, "Supplementary_Figure5_Rhizosphere_VOCs_heatmap.pdf"),
  width = 3,
  height = 5
)

# -----------------------------------------------------------------------------
# 4. Supplementary Figure 6: rhizosphere non-volatile LC-MS heatmap
# -----------------------------------------------------------------------------

lcms_matrix <- read_feature_matrix(lcms_file)

# Original exploratory script used three NoH and three Herbivory labels.
# TODO: Confirm final sample order in `data_norm_visual.csv`.
lcms_conditions <- c(
  rep("NoHerb", 3),
  rep("Herbivory", 3)
)

save_heatmap(
  mat = lcms_matrix,
  conditions = lcms_conditions,
  output_file = file.path(FIGURE_DIR, "Supplementary_Figure6_Rhizosphere_nonvolatile_heatmap.pdf"),
  width = 7,
  height = 9
)

write_session_info(file.path(OUTPUT_DIR, "session_info_rhizosphere_metabolite_heatmaps.txt"))
