###############################################################################
# Root dry weight supplementary figure
# Manuscript: Aphid herbivory transiently restructures rhizosphere chemistry
#             and bacterial communities in wheat
# Author: Sandra Cortes-Patiño
# Tidied for GitHub/reproducibility
###############################################################################

source("scripts/00_project_setup.R")

packages <- c("dplyr", "ggplot2", "ggpubr", "rstatix")
invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. Input/output paths
# -----------------------------------------------------------------------------

root_file <- file.path(DATA_DIR, "microbiome", "roots_dryweight.csv")
root_out <- file.path(OUTPUT_DIR, "root_dry_weight")
root_fig <- file.path(FIGURE_DIR, "root_dry_weight")

dir.create(root_out, recursive = TRUE, showWarnings = FALSE)
dir.create(root_fig, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------

root_df <- read.csv(root_file, check.names = FALSE)

# Expected columns:
#   Time:      Two weeks / Four weeks
#   Treatment: NoHerb / Herbivory
#   Root:      root dry weight in grams
# TODO: Check that these column names match the final archived CSV file.

root_df <- root_df |>
  dplyr::mutate(
    Time = factor(Time, levels = c("Two weeks", "Four weeks")),
    Treatment = factor(Treatment, levels = c("NoHerb", "Herbivory"))
  )

# -----------------------------------------------------------------------------
# 3. Statistics
# -----------------------------------------------------------------------------

root_stats <- root_df |>
  dplyr::group_by(Time) |>
  rstatix::wilcox_test(Root ~ Treatment) |>
  dplyr::mutate(
    label = paste0("p = ", signif(p, 2)),
    x = 1.5,
    y = max(root_df$Root, na.rm = TRUE) * 1.05
  )

write.csv(root_stats, file.path(root_out, "root_dry_weight_wilcox_tests.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Plot
# -----------------------------------------------------------------------------

p_root <- ggplot(root_df, aes(x = Treatment, y = Root, fill = Treatment)) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.85,
    colour = "black",
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = Treatment),
    width = 0.08,
    size = 3,
    alpha = 0.9
  ) +
  geom_text(
    data = root_stats,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~Time, nrow = 1) +
  scale_fill_manual(values = condition_colours) +
  scale_colour_manual(values = condition_colours) +
  labs(x = NULL, y = "Root dry weight (g)") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    strip.text = element_text(face = "bold")
  )

p_root

ggsave(
  filename = file.path(root_fig, "supplementary_root_dry_weight.pdf"),
  plot = p_root,
  width = 16,
  height = 10,
  units = "cm",
  device = cairo_pdf
)

write_session_info(file.path(root_out, "session_info_root_dry_weight.txt"))
