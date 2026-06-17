#roots
setwd("C:/Users/SandraCortes/OneDrive - Bactobio Ltd/Documents/Man/Phylo_M_filtered-20251013T190827Z-1-001/Phylo_M_filtered")


root_df <- read.csv("roots_dryweight.csv")

library(ggplot2)
library(ggpubr)
library(dplyr)

root_df <- root_df %>%
  mutate(
    Time = factor(Time, levels = c("Two weeks", "Four weeks")),
    Treatment = factor(Treatment, levels = c("NoHerb", "Herbivory"))
  )

p_root <- ggplot(
  root_df,
  aes(x = Treatment, y = Root, fill = Treatment)
) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.8,
    colour = "black",
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = Treatment),
    width = 0.08,
    size = 3,
    alpha = 0.9
  ) +
  stat_compare_means(
    comparisons = list(c("NoHerb", "Herbivory")),
    method = "wilcox.test",
    label = "p.format",
    label.y = 0.5,   # adjust for your data
  )+
  facet_wrap(~Time, nrow = 1) +
  scale_fill_manual(
    values = c(
      "NoHerb" = "#91bfdb",
      "Herbivory" = "#fc8d59"
    )
  ) +
  scale_colour_manual(
    values = c(
      "NoHerb" = "#91bfdb",
      "Herbivory" = "#fc8d59"
    )
  ) +
  labs(
    x = NULL,
    y = "Root dry weight (g)"
  ) +
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
  filename = "Supplementary_Figure_RootDryWeight.pdf",
  plot = p_root,
  width = 16,
  height = 10,
  units = "cm",
  device = cairo_pdf
)