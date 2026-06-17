###############################################################################
# Project setup and shared plotting style
# Manuscript: Aphid herbivory transiently restructures rhizosphere chemistry
#             and bacterial communities in wheat
# Author: Sandra Cortes-Patiño
# Purpose: Shared paths, colours, and plotting theme used by analysis scripts.
###############################################################################

# This repository is intended to document the analyses used in the manuscript.
# It is not an R package. Scripts assume that input data are placed in the
# folder structure described in data/README_data.md.

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

DATA_DIR <- "data"
OUTPUT_DIR <- "outputs"
FIGURE_DIR <- "figures"

for (dir in c(OUTPUT_DIR, FIGURE_DIR)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------------------------------------------------------
# 2. Shared colours
# -----------------------------------------------------------------------------

condition_colours <- c(
  "Herbivory" = "#fc8d59",
  "NoHerb" = "#91bfdb",
  "Herbivory-2weeks" = "#fc8d59",
  "Herbivory-4weeks" = "tomato",
  "NoHerb-2weeks" = "#91bfdb",
  "NoHerb-4weeks" = "steelblue",
  "Before" = "#61576A",
  "Bulk" = "grey50"
)

# -----------------------------------------------------------------------------
# 3. Shared ggplot theme
# -----------------------------------------------------------------------------

manuscript_theme <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      legend.position = "right",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# -----------------------------------------------------------------------------
# 4. Utility functions
# -----------------------------------------------------------------------------

write_session_info <- function(path = file.path(OUTPUT_DIR, "session_info.txt")) {
  sink(path)
  print(sessionInfo())
  sink()
}

check_sample_order <- function(data, metadata) {
  if (!all(rownames(data) == rownames(metadata))) {
    stop("Sample names/order do not match between data matrix and metadata.")
  }
  invisible(TRUE)
}
