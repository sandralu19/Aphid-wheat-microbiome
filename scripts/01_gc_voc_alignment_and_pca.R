###############################################################################
# GC VOC alignment, normalisation, PCA and PERMANOVA
# Manuscript: Aphid herbivory transiently restructures rhizosphere chemistry...
# Author: Sandra Cortes
# Tidied version for reproducibility/GitHub
###############################################################################

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

packages <- c(
  "GCalignR",
  "readxl",
  "dplyr",
  "vegan",
  "FactoMineR",
  "factoextra",
  "ggplot2"
)

invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. User-defined paths
# -----------------------------------------------------------------------------
# Change these paths to match your local folder structure.
# Avoid setwd() so the script can be run reproducibly from any location.

aboveground_dir <- "data/gc_aboveground"
rhizosphere_dir <- "data/gc_rhizosphere_vocs"
output_dir <- "outputs/gc_vocs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

load_chromatograms <- function(data_dir, sample_files, blank_file = NULL) {
  chrom_list <- lapply(file.path(data_dir, sample_files), read.delim)
  chrom_list <- lapply(chrom_list, as.data.frame)
  names(chrom_list) <- names(sample_files)
  
  if (!is.null(blank_file)) {
    blank <- readxl::read_xlsx(file.path(data_dir, blank_file)) |>
      as.data.frame()
    chrom_list <- c(chrom_list, list(Blank = blank))
  }
  
  chrom_list
}

align_gc_peaks <- function(chrom_list,
                           reference = NULL,
                           blank_name = NULL,
                           rt_low = 5,
                           rt_high = 40,
                           max_linear_shift = 0.05,
                           max_diff_peak2mean = 0.03,
                           min_diff_peak2peak = 0.03) {
  
  check_input(chrom_list, plot = TRUE)
  
  peak_interspace(
    data = chrom_list,
    rt_col_name = "RT",
    quantile_range = c(0, 0.8),
    quantiles = 0.05
  )
  
  align_args <- list(
    data = chrom_list,
    rt_col_name = "RT",
    rt_cutoff_low = rt_low,
    rt_cutoff_high = rt_high,
    reference = reference,
    max_linear_shift = max_linear_shift,
    max_diff_peak2mean = max_diff_peak2mean,
    min_diff_peak2peak = min_diff_peak2peak,
    delete_single_peak = TRUE,
    write_output = NULL
  )
  
  if (!is.null(blank_name)) {
    align_args$blanks <- blank_name
  }
  
  do.call(GCalignR::align_chromatograms, align_args)
}

normalise_gc_peaks <- function(aligned_object) {
  norm_peaks(
    aligned_object,
    conc_col_name = "Area",
    rt_col_name = "RT",
    out = "data.frame"
  ) |>
    log1p()
}

match_metadata_order <- function(data_matrix, metadata) {
  matched <- data_matrix[match(rownames(metadata), rownames(data_matrix)), , drop = FALSE]
  stopifnot(all(rownames(matched) == rownames(metadata)))
  matched
}

run_pca_plot <- function(data_matrix, metadata, colour_col, shape_col = NULL, title = NULL) {
  pca_result <- FactoMineR::PCA(data_matrix, graph = FALSE)
  
  eig <- factoextra::get_eig(pca_result)
  print(eig)
  
  p <- factoextra::fviz_pca_ind(
    pca_result,
    geom.ind = "point",
    col.ind = metadata[[colour_col]],
    addEllipses = TRUE,
    ellipse.type = "confidence",
    legend.title = "Groups"
  ) +
    scale_colour_manual(
      values = c(
        "Herbivory" = "#fc8d59",
        "NoHerb" = "#91bfdb",
        "Herbivory-2weeks" = "#fc8d59",
        "Herbivory-4weeks" = "tomato",
        "NoHerb-2weeks" = "#91bfdb",
        "NoHerb-4weeks" = "steelblue"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Herbivory" = "#fc8d59",
        "NoHerb" = "#91bfdb",
        "Herbivory-2weeks" = "#fc8d59",
        "Herbivory-4weeks" = "tomato",
        "NoHerb-2weeks" = "#91bfdb",
        "NoHerb-4weeks" = "steelblue"
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "right"
    ) +
    labs(title = title)
  
  if (!is.null(shape_col)) {
    p <- p + aes(shape = metadata[[shape_col]])
  }
  
  list(pca = pca_result, plot = p, eigenvalues = eig)
}

run_permanova <- function(data_matrix, metadata, formula_text, method = "bray") {
  formula_obj <- as.formula(formula_text)
  vegan::adonis2(
    formula_obj,
    data = metadata,
    permutations = 999,
    method = method
  )
}

# -----------------------------------------------------------------------------
# 3. Aboveground VOCs
# -----------------------------------------------------------------------------

aboveground_files <- c(
  AG9A  = "AGA9.txt",
  AG5A  = "AG5.txt",
  AG14A = "A14.txt",
  AG2A  = "A2.txt",
  AG4C  = "C4.txt",
  AG10C = "C10.txt",
  AG15C = "C15.txt",
  AG12C = "A12.txt",
  LA7   = "LA7.txt",
  LA8   = "LA8.txt",
  LA12  = "LA12.txt",
  LA13  = "LA13.txt",
  LC6   = "LC6.txt",
  LC7   = "LC7.txt",
  LC8   = "LC8.txt",
  LC13  = "LC13.txt",
  LC14  = "LC14.txt"
)

aboveground_chrom <- load_chromatograms(
  data_dir = aboveground_dir,
  sample_files = aboveground_files,
  blank_file = "Blank.xlsx"
)

aboveground_aligned <- align_gc_peaks(
  chrom_list = aboveground_chrom,
  reference = "AG10C",
  blank_name = "Blank",
  rt_low = 5,
  rt_high = 40,
  max_linear_shift = 0.05,
  max_diff_peak2mean = 0.03,
  min_diff_peak2peak = 0.03
)

# Quality-check plots
GCalignR::gc_heatmap(aboveground_aligned)
plot(aboveground_aligned, which_plot = "all")

aboveground_area <- aboveground_aligned$aligned$Area
write.csv(
  aboveground_area,
  file.path(output_dir, "aboveground_vocs_aligned_area.csv"),
  row.names = FALSE
)

aboveground_norm <- normalise_gc_peaks(aboveground_aligned)
write.csv(
  aboveground_norm,
  file.path(output_dir, "aboveground_vocs_log_normalised.csv")
)

aboveground_metadata <- read.csv(
  file.path(aboveground_dir, "metadata.csv"),
  row.names = 1
)

aboveground_norm <- match_metadata_order(aboveground_norm, aboveground_metadata)

aboveground_permanova <- run_permanova(
  data_matrix = aboveground_norm,
  metadata = aboveground_metadata,
  formula_text = "aboveground_norm ~ Treatment * Sampling",
  method = "bray"
)

write.csv(
  as.data.frame(aboveground_permanova),
  file.path(output_dir, "aboveground_vocs_permanova.csv")
)

aboveground_pca <- run_pca_plot(
  data_matrix = aboveground_norm,
  metadata = aboveground_metadata,
  colour_col = "Treatment",
  shape_col = "Sampling",
  title = "Leaf VOC profiles"
)

aboveground_pca$plot

ggsave(
  filename = file.path(output_dir, "aboveground_vocs_pca.pdf"),
  plot = aboveground_pca$plot,
  width = 16,
  height = 10,
  units = "cm",
  device = cairo_pdf
)

# -----------------------------------------------------------------------------
# 4. Rhizosphere VOCs
# -----------------------------------------------------------------------------

rhizosphere_files <- c(
  b14a = "Below_14.txt",
  b2a  = "NBelow_2A.txt",
  b5a  = "Below_A5.txt",
  b6a  = "Below_A6.txt",
  b9a  = "Below_A9.txt",
  b3c  = "Below_C3.txt",
  b10c = "Below_C10.txt",
  b12c = "Below_C12.txt",
  b14c = "Below_C14.txt",
  b15c = "Below_C15.txt"
)

rhizosphere_chrom <- load_chromatograms(
  data_dir = rhizosphere_dir,
  sample_files = rhizosphere_files,
  blank_file = NULL
)

rhizosphere_aligned <- align_gc_peaks(
  chrom_list = rhizosphere_chrom,
  reference = NULL,
  blank_name = NULL,
  rt_low = 5,
  rt_high = 40,
  max_linear_shift = 0.05,
  max_diff_peak2mean = 0.03,
  min_diff_peak2peak = 0.03
)

GCalignR::gc_heatmap(rhizosphere_aligned)
plot(rhizosphere_aligned, which_plot = "all")

rhizosphere_area <- rhizosphere_aligned$aligned$Area
write.csv(
  rhizosphere_area,
  file.path(output_dir, "rhizosphere_vocs_aligned_area.csv"),
  row.names = FALSE
)

rhizosphere_norm <- normalise_gc_peaks(rhizosphere_aligned)
write.csv(
  rhizosphere_norm,
  file.path(output_dir, "rhizosphere_vocs_log_normalised.csv")
)

rhizosphere_metadata <- read.csv(
  file.path(rhizosphere_dir, "metadata.csv"),
  row.names = 1
)

rhizosphere_norm <- match_metadata_order(rhizosphere_norm, rhizosphere_metadata)

rhizosphere_permanova <- run_permanova(
  data_matrix = rhizosphere_norm,
  metadata = rhizosphere_metadata,
  formula_text = "rhizosphere_norm ~ Treatment",
  method = "bray"
)

write.csv(
  as.data.frame(rhizosphere_permanova),
  file.path(output_dir, "rhizosphere_vocs_permanova.csv")
)

rhizosphere_pca <- run_pca_plot(
  data_matrix = rhizosphere_norm,
  metadata = rhizosphere_metadata,
  colour_col = "Treatment",
  shape_col = NULL,
  title = "Rhizosphere VOC profiles"
)

rhizosphere_pca$plot

ggsave(
  filename = file.path(output_dir, "rhizosphere_vocs_pca.pdf"),
  plot = rhizosphere_pca$plot,
  width = 14,
  height = 10,
  units = "cm",
  device = cairo_pdf
)

# -----------------------------------------------------------------------------
# 5. Session information
# -----------------------------------------------------------------------------

sink(file.path(output_dir, "session_info.txt"))
sessionInfo()
sink()
