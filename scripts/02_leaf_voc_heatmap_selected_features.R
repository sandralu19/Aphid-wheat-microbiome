###############################################################################
# Heatmap of selected leaf VOC features
# Author: Sandra Cortes
# Tidied version for reproducibility/GitHub
###############################################################################

library(dplyr)
library(ComplexHeatmap)
library(circlize)

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

input_file <- "data/gc_aboveground/100424_selectedpeaks.csv"
output_dir <- "outputs/gc_vocs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load selected peak matrix
# -----------------------------------------------------------------------------
# Expected format:
# - first column = metabolite/feature name
# - remaining columns = samples

selected_vocs <- read.csv(input_file, check.names = FALSE)
rownames(selected_vocs) <- selected_vocs[[1]]
selected_vocs <- selected_vocs[, -1, drop = FALSE]

# Z-score by feature, then transpose so rows = metabolites and columns = samples
scaled_matrix <- t(scale(t(as.matrix(selected_vocs))))

# -----------------------------------------------------------------------------
# 3. Feature annotations
# -----------------------------------------------------------------------------
# Check that the order below matches rownames(scaled_matrix).

feature_annotations <- data.frame(
  Metabolite = c(
    "Unknown_739", "E-2Hexen-1-ol", "1,3 Dimethylbenzene",
    "5-Methyl-3-methylene-5-hexen-2-one", "Hexyl acetate", "Benzyl alcohol",
    "Limonene", "Z-Ocimene", "Acetophenone", "Isophorone",
    "alpha-Methylenephenylacetaldehyde", "Dodecene", "Myrtenol",
    "Z-3-Hexen-1-yl-3-methylbutyrate", "Indole", "alpha-cubebene",
    "3-Methylindole", "(Z)-Jasmone", "Citronellyl propionate",
    "Methyl jasmonate", "Unknown_1631", "Caryophyllene oxide", "Unknown_1867"
  ),
  Pathway = c(
    "Unknown", "LOX", "Shikimate", "LOX", "LOX", "LOX",
    "MEV/Non-MEV", "MEV/Non-MEV", "Shikimate", "LOX", "Shikimate",
    "LOX", "MEV/Non-MEV", "LOX", "Shikimate", "MEV/Non-MEV",
    "Shikimate", "LOX", "MEV/Non-MEV", "LOX", "Unknown",
    "MEV/Non-MEV", "Unknown"
  ),
  Function = c(
    "Unknown", "Signalling", "No", "Limited", "Signalling", "Antimicrobial",
    "Direct defence", "Indirect defence", "Limited", "Limited", "Limited",
    "Limited", "Direct defence", "Signalling", "Multifunctional",
    "Indirect defence", "Multifunctional", "Multifunctional", "Limited",
    "Signalling", "Unknown", "Direct defence", "Unknown"
  )
)

# Optional safety check if your metabolite names match row names exactly
# stopifnot(all(feature_annotations$Metabolite == rownames(scaled_matrix)))

# -----------------------------------------------------------------------------
# 4. Sample annotations
# -----------------------------------------------------------------------------

sample_metadata <- data.frame(
  Sample = c(
    "Herb1", "Herb2", "Herb3", "Herb4",
    "NoHerb1", "NoHerb2", "NoHerb3", "NoHerb4"
  ),
  Condition = c(
    rep("Herbivory", 4),
    rep("NoHerb", 4)
  )
)

# If column names match sample names, use this check.
# stopifnot(all(colnames(scaled_matrix) == sample_metadata$Sample))

# -----------------------------------------------------------------------------
# 5. Colours and annotations
# -----------------------------------------------------------------------------

heatmap_colours <- colorRamp2(
  c(-2, 0, 2),
  c("white", "#FFE1D4", "#F53D2A")
)

pathway_colours <- c(
  "Unknown" = "#99bc84",
  "LOX" = "#008ea8",
  "Shikimate" = "#ee6a66",
  "MEV/Non-MEV" = "#a54fc5"
)

function_colours <- c(
  "Signalling" = "#5671d4",
  "No" = "#8a92a2",
  "Indirect defence" = "#93003a",
  "Direct defence" = "#78dac9",
  "Multifunctional" = "#ffa77b",
  "Antimicrobial" = "#ee6a66",
  "Limited" = "#ca2f50",
  "Unknown" = "#8db6a7"
)

condition_colours <- c(
  "Herbivory" = "#fc8d59",
  "NoHerb" = "#91bfdb"
)

row_anno <- rowAnnotation(
  Pathway = feature_annotations$Pathway,
  Function = feature_annotations$Function,
  col = list(
    Pathway = pathway_colours,
    Function = function_colours
  )
)

col_anno <- HeatmapAnnotation(
  Condition = sample_metadata$Condition,
  col = list(Condition = condition_colours)
)

# -----------------------------------------------------------------------------
# 6. Plot and save heatmap
# -----------------------------------------------------------------------------

heatmap_plot <- Heatmap(
  scaled_matrix,
  name = "Z-score",
  top_annotation = col_anno,
  left_annotation = row_anno,
  col = heatmap_colours,
  show_row_dend = FALSE,
  cluster_rows = FALSE,
  row_split = feature_annotations$Pathway,
  column_split = sample_metadata$Condition,
  show_row_names = TRUE,
  show_column_names = FALSE
)

pdf(
  file.path(output_dir, "leaf_voc_selected_features_heatmap.pdf"),
  width = 8,
  height = 7
)
draw(heatmap_plot)
dev.off()

# -----------------------------------------------------------------------------
# 7. Session information
# -----------------------------------------------------------------------------

sink(file.path(output_dir, "leaf_voc_heatmap_session_info.txt"))
sessionInfo()
sink()
