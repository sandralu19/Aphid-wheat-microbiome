###############################################################################
# Biolog EcoPlate analysis
# Manuscript: Aphid herbivory transiently restructures rhizosphere chemistry
#             and bacterial communities in wheat
# Author: Sandra Cortes-Patiño
# Purpose: Analyse AWCD dynamics and substrate-level utilisation for Figure 6.
###############################################################################

# This script was tidied from exploratory EcoPlate analyses. It keeps the final
# analysis steps used in the manuscript: AWCD summaries, treatment x time model,
# post-hoc grouping, substrate-level tests at 72 h, and Figure 6 export.

source("scripts/00_project_setup.R")

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(agricolae)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Input files
# -----------------------------------------------------------------------------

# Expected files:
# data/ecoplates/231023_summary.xlsx
#   Long-format AWCD data with columns similar to:
#   id, Treatment, Hour, AWCD
#
# data/ecoplates/Ecoplates_complete.xlsx, sheet = "Mean_nosoil"
#   Long-format substrate-level data with columns similar to:
#   Treatment, Type, Name, AWCD
#
# TODO: Confirm final file names and sheet names before making the repository public.

ecoplate_dir <- file.path(DATA_DIR, "ecoplates")

awcd_file <- file.path(ecoplate_dir, "231023_summary.xlsx")
substrate_file <- file.path(ecoplate_dir, "Ecoplates_complete.xlsx")

# -----------------------------------------------------------------------------
# 2. Load and format AWCD data
# -----------------------------------------------------------------------------

awcd_raw <- read_excel(awcd_file) %>%
  as.data.frame()

# Standardise common treatment labels used during exploratory analysis.
awcd <- awcd_raw %>%
  mutate(
    Treatment = recode(
      as.character(Treatment),
      "Insect" = "Herbivory",
      "NoI" = "NoHerb",
      "Bsoil" = "Bulk soil",
      .default = as.character(Treatment)
    ),
    Treatment = factor(Treatment, levels = c("Bulk soil", "NoHerb", "Herbivory")),
    Hour = as.numeric(Hour),
    id = factor(id)
  )

# -----------------------------------------------------------------------------
# 3. Summarise AWCD by treatment and time
# -----------------------------------------------------------------------------

awcd_summary <- awcd %>%
  group_by(Treatment, Hour) %>%
  summarise(
    n = sum(!is.na(AWCD)),
    mean_awcd = mean(AWCD, na.rm = TRUE),
    sd = sd(AWCD, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

write.csv(
  awcd_summary,
  file.path(OUTPUT_DIR, "ecoplate_awcd_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Treatment x time model and post-hoc grouping
# -----------------------------------------------------------------------------

# The manuscript reports a two-way ANOVA with treatment, hour and interaction.
# The original exploratory script also checked repeated-measures/mixed models.
# Here we keep the final ANOVA/SNK workflow used for the manuscript figure.

awcd_aov <- aov(AWCD ~ Treatment * Hour, data = awcd)
awcd_anova <- as.data.frame(summary(awcd_aov)[[1]])

write.csv(
  awcd_anova,
  file.path(OUTPUT_DIR, "ecoplate_awcd_two_way_anova.csv")
)

# Student-Newman-Keuls post hoc grouping for Treatment x Hour combinations.
snk <- SNK.test(
  awcd_aov,
  trt = c("Treatment", "Hour"),
  console = FALSE
)

snk_groups <- snk$groups %>%
  tibble::rownames_to_column("Treatment_Hour") %>%
  rename(letters = groups)

write.csv(
  snk_groups,
  file.path(OUTPUT_DIR, "ecoplate_awcd_snk_groups.csv"),
  row.names = FALSE
)

# Optional: join letters to summary if labels are required in the AWCD line plot.
# The combined label format from agricolae may depend on factor coding. If labels
# do not join automatically, manually check `ecoplate_awcd_snk_groups.csv`.
awcd_summary <- awcd_summary %>%
  mutate(Treatment_Hour = paste(Treatment, Hour, sep = ":")) %>%
  left_join(snk_groups, by = "Treatment_Hour")

# -----------------------------------------------------------------------------
# 5. Substrate-level utilisation at 72 h
# -----------------------------------------------------------------------------

substrates <- read_excel(substrate_file, sheet = "Mean_nosoil") %>%
  as.data.frame() %>%
  mutate(
    Treatment = recode(
      as.character(Treatment),
      "Insect" = "Herbivory",
      "NoI" = "NoHerb",
      .default = as.character(Treatment)
    ),
    Treatment = factor(Treatment, levels = c("NoHerb", "Herbivory"))
  )

# Substrate-level tests. This assumes the input sheet contains replicate-level or
# suitable summary-level AWCD values at 72 h. If the sheet contains only group
# means, replace this section with the replicate-level substrate table.
substrate_stats <- substrates %>%
  filter(Treatment %in% c("NoHerb", "Herbivory")) %>%
  group_by(Type, Name) %>%
  wilcox_test(AWCD ~ Treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  ungroup()

write.csv(
  substrate_stats,
  file.path(OUTPUT_DIR, "ecoplate_substrate_wilcox_tests.csv"),
  row.names = FALSE
)

substrates_plot <- substrates %>%
  left_join(
    substrate_stats %>% select(Type, Name, p.adj, p.adj.signif),
    by = c("Type", "Name")
  )

# -----------------------------------------------------------------------------
# 6. Generate Figure 6 panels
# -----------------------------------------------------------------------------

awcd_colours <- c(
  "Bulk soil" = "#61576A",
  "NoHerb" = "#91bfdb",
  "Herbivory" = "#fc8d59"
)

p_awcd <- ggplot(
  awcd_summary,
  aes(x = Hour, y = mean_awcd, colour = Treatment, group = Treatment)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean_awcd - se, ymax = mean_awcd + se),
    width = 2,
    linewidth = 0.6
  ) +
  scale_colour_manual(values = awcd_colours, drop = FALSE) +
  labs(x = "Time (h)", y = "AWCD", colour = "Treatment") +
  manuscript_theme(base_size = 13)

p_substrate <- ggplot(
  substrates_plot,
  aes(x = Name, y = Treatment, fill = AWCD)
) +
  geom_tile(colour = "black", linewidth = 0.2) +
  geom_text(
    data = substrates_plot %>% filter(!is.na(p.adj.signif), p.adj < 0.05),
    aes(label = p.adj.signif),
    colour = "black",
    size = 4
  ) +
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "white", high = "#F53D2A", name = "AWCD") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

figure6 <- p_awcd / p_substrate +
  plot_annotation(tag_levels = "a")

figure6

ggsave(
  filename = file.path(FIGURE_DIR, "Figure6_Ecoplates.pdf"),
  plot = figure6,
  width = 22,
  height = 16,
  units = "cm",
  device = cairo_pdf
)

write_session_info(file.path(OUTPUT_DIR, "session_info_ecoplate_analysis.txt"))
