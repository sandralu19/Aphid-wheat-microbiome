###############################################################################
# Microbiome analysis: 16S rRNA rhizosphere communities
# Manuscript: Aphid herbivory transiently restructures rhizosphere chemistry
#             and bacterial communities in wheat
# Author: Sandra Cortes-Patiño
# Tidied for GitHub/reproducibility
###############################################################################

# This script documents the microbiome workflow used in the manuscript:
#   - import ASV, taxonomy, metadata and tree files
#   - create phyloseq object
#   - remove chloroplast/mitochondrial ASVs
#   - rarefaction curves and alpha diversity
#   - prevalence filtering for beta diversity
#   - NMDS, PERMANOVA and dispersion testing
#   - partial dbRDA
#   - Venn diagrams for ASV turnover
#   - ANCOM-BC differential abundance at ASV level

# TODO: Confirm that the input file names below match the final archived files.
# TODO: If the raw sequencing reads are deposited in SRA/ENA, add accession in README.

source("scripts/00_project_setup.R")

packages <- c(
  "phyloseq", "readxl", "dplyr", "tidyr", "ggplot2", "vegan",
  "ggpubr", "rstatix", "patchwork", "VennDiagram", "grid", "gridExtra", "ANCOMBC",
  "mia", "miaViz", "SummarizedExperiment"
)

invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. Input/output paths
# -----------------------------------------------------------------------------

microbiome_dir <- file.path(DATA_DIR, "microbiome")
microbiome_out <- file.path(OUTPUT_DIR, "microbiome")
microbiome_fig <- file.path(FIGURE_DIR, "microbiome")

dir.create(microbiome_out, recursive = TRUE, showWarnings = FALSE)
dir.create(microbiome_fig, recursive = TRUE, showWarnings = FALSE)

asv_file <- file.path(microbiome_dir, "ASVs.xlsx")
tax_file <- file.path(microbiome_dir, "TAX_K.xlsx")
metadata_file <- file.path(microbiome_dir, "MET.xlsx")
tree_file <- file.path(microbiome_dir, "tree.nwk")

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

import_phyloseq_from_excel <- function(asv_file, tax_file, metadata_file, tree_file = NULL) {
  otu_mat <- readxl::read_excel(asv_file) |> as.data.frame()
  tax_mat <- readxl::read_excel(tax_file) |> as.data.frame()
  samples_df <- readxl::read_excel(metadata_file) |> as.data.frame()

  # First column is assumed to contain feature/sample identifiers.
  rownames(otu_mat) <- otu_mat[[1]]
  otu_mat <- otu_mat[, -1, drop = FALSE]
  otu_mat <- data.matrix(otu_mat)

  rownames(tax_mat) <- tax_mat[[1]]
  tax_mat <- tax_mat[, -1, drop = FALSE]
  tax_mat <- as.matrix(tax_mat)

  rownames(samples_df) <- samples_df[[1]]
  samples_df <- samples_df[, -1, drop = FALSE]

  physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax_mat),
    phyloseq::sample_data(samples_df)
  )

  if (!is.null(tree_file) && file.exists(tree_file)) {
    physeq <- phyloseq::merge_phyloseq(physeq, phyloseq::read_tree(tree_file, errorIfNULL = FALSE))
  }

  physeq
}

remove_non_bacterial_asvs <- function(physeq) {
  physeq |>
    phyloseq::subset_taxa(!(Family %in% c("f_Mitochondria", "Mitochondria"))) |>
    phyloseq::subset_taxa(!(Order %in% c("o_Chloroplast", "Chloroplast")))
}

filter_by_presence <- function(physeq, sampling_group, threshold = 50) {
  # Keeps taxa present above `threshold` percent of samples in at least one group.
  if (!phyloseq::taxa_are_rows(physeq)) {
    phyloseq::otu_table(physeq) <- t(phyloseq::otu_table(physeq))
  }

  group_vector <- phyloseq::sample_data(physeq)[[sampling_group]]
  groups <- unique(as.character(group_vector))
  cols_per_group <- lapply(groups, function(x) which(group_vector == x))
  names(cols_per_group) <- groups

  presence_matrix <- t(apply(phyloseq::otu_table(physeq), 1, function(x) {
    sapply(cols_per_group, function(y) 100 - (sum(x[y] == 0) / length(y) * 100))
  }))

  keep_taxa <- rownames(presence_matrix)[rowSums(presence_matrix > threshold) > 0]
  phyloseq::prune_taxa(keep_taxa, physeq)
}

save_pdf <- function(plot, filename, width = 16, height = 10) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "cm",
    device = cairo_pdf
  )
}

# -----------------------------------------------------------------------------
# 3. Import, clean and save phyloseq object
# -----------------------------------------------------------------------------

aphids <- import_phyloseq_from_excel(asv_file, tax_file, metadata_file, tree_file)
filtered16S <- remove_non_bacterial_asvs(aphids)

saveRDS(filtered16S, file.path(microbiome_out, "filtered16S_phyloseq.rds"))

# -----------------------------------------------------------------------------
# 4. Alpha rarefaction and alpha diversity
# -----------------------------------------------------------------------------

set.seed(123)
ps_rarefied <- phyloseq::rarefy_even_depth(
  filtered16S,
  sample.size = min(phyloseq::sample_sums(filtered16S)),
  replace = FALSE,
  rngseed = 123
)

# Rarefaction curves: before and after rarefaction.
pdf(file.path(microbiome_fig, "supplementary_rarefaction_before.pdf"), width = 6, height = 4)
otu_before <- t(as(phyloseq::otu_table(filtered16S), "matrix"))
vegan::rarecurve(otu_before, step = 20, sample = min(rowSums(otu_before)), col = "#998ec3", cex = 0.6)
dev.off()

pdf(file.path(microbiome_fig, "supplementary_rarefaction_after.pdf"), width = 6, height = 4)
otu_after <- t(as(phyloseq::otu_table(ps_rarefied), "matrix"))
vegan::rarecurve(otu_after, step = 20, sample = min(rowSums(otu_after)), col = "#998ec3", cex = 0.6)
dev.off()

richness <- phyloseq::estimate_richness(ps_rarefied, measures = c("Observed", "Shannon"))
richness$SampleID <- rownames(richness)
metadata_alpha <- phyloseq::sample_data(ps_rarefied) |> as.data.frame()
metadata_alpha$SampleID <- rownames(metadata_alpha)

alpha_df <- dplyr::left_join(richness, metadata_alpha, by = "SampleID") |>
  tidyr::pivot_longer(cols = c("Observed", "Shannon"), names_to = "Index", values_to = "Value")

# TODO: Check that Rep and Insect levels match the final metadata.
alpha_df$Rep <- factor(alpha_df$Rep, levels = c("Before", "Herb2", "NoH2", "Herb4", "NoH4", "Bulk"))

alpha_stats <- alpha_df |>
  dplyr::group_by(Index) |>
  rstatix::t_test(Value ~ Rep) |>
  dplyr::filter(p.adj < 0.05) |>
  rstatix::add_xy_position(x = "Rep")

alpha_plot <- ggplot(alpha_df, aes(x = Rep, y = Value, fill = Insect)) +
  geom_boxplot(outlier.shape = NA, colour = "black") +
  geom_jitter(width = 0.08, size = 2, alpha = 0.8) +
  facet_wrap(~Index, scales = "free_y") +
  ggpubr::stat_pvalue_manual(alpha_stats, label = "p.adj.signif", tip.length = 0.01) +
  scale_fill_manual(values = condition_colours, na.value = "grey70") +
  labs(x = NULL, y = NULL) +
  manuscript_theme(base_size = 12) +
  theme(legend.position = "right")

save_pdf(alpha_plot, file.path(microbiome_fig, "figure_alpha_diversity.pdf"), width = 18, height = 9)
write.csv(alpha_df, file.path(microbiome_out, "alpha_diversity_values.csv"), row.names = FALSE)
write.csv(alpha_stats, file.path(microbiome_out, "alpha_diversity_pairwise_tests.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Prevalence filtering and beta diversity
# -----------------------------------------------------------------------------

filtered_prevalent <- filter_by_presence(filtered16S, sampling_group = "Rep", threshold = 50)
saveRDS(filtered_prevalent, file.path(microbiome_out, "filtered16S_prevalence50_phyloseq.rds"))

relative_ps <- phyloseq::transform_sample_counts(filtered_prevalent, function(x) x / sum(x))
bray_dist <- phyloseq::distance(relative_ps, method = "bray")
metadata_beta <- phyloseq::sample_data(relative_ps) |> as.data.frame()

# NMDS with vegan, to report stress directly.
set.seed(123)
nmds <- vegan::metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE)
nmds_scores <- as.data.frame(vegan::scores(nmds, display = "sites"))
nmds_scores$SampleID <- rownames(nmds_scores)
metadata_beta$SampleID <- rownames(metadata_beta)
nmds_scores <- dplyr::left_join(nmds_scores, metadata_beta, by = "SampleID")

nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = Insect, shape = Time)) +
  stat_ellipse(aes(group = interaction(Insect, Time)), level = 0.95, linewidth = 0.7, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.95) +
  annotate("text", x = Inf, y = -Inf, label = paste0("Stress = ", round(nmds$stress, 3)), hjust = 1.1, vjust = -0.8, size = 4) +
  scale_colour_manual(values = condition_colours, na.value = "grey70") +
  coord_equal() +
  labs(x = "NMDS1", y = "NMDS2", colour = "Treatment", shape = "Timepoint") +
  manuscript_theme()

save_pdf(nmds_plot, file.path(microbiome_fig, "figure_nmds_bray.pdf"), width = 16, height = 10)

permanova_all <- vegan::adonis2(bray_dist ~ Insect * Time, data = metadata_beta, by = "terms", permutations = 999)
dispersion <- vegan::betadisper(bray_dist, metadata_beta$Insect)
dispersion_test <- vegan::permutest(dispersion)

write.csv(as.data.frame(permanova_all), file.path(microbiome_out, "permanova_bray_all.csv"))
capture.output(dispersion_test, file = file.path(microbiome_out, "betadisper_insect.txt"))

# -----------------------------------------------------------------------------
# 6. Partial dbRDA
# -----------------------------------------------------------------------------

# The manuscript focuses on weeks 2 and 4 for the herbivory effect.
# TODO: Confirm that Time values are exactly those used in your final metadata.
ps_rhizosphere <- phyloseq::subset_samples(relative_ps, SampleType %in% c("Rhizosphere"))
metadata_rda <- phyloseq::sample_data(ps_rhizosphere) |> as.data.frame()
bray_rda <- phyloseq::distance(ps_rhizosphere, method = "bray")

partial_dbrda <- vegan::capscale(bray_rda ~ Insect + Condition(Time), data = metadata_rda)
partial_dbrda_anova <- vegan::anova.cca(partial_dbrda, by = "terms", permutations = 999)

write.csv(as.data.frame(partial_dbrda_anova), file.path(microbiome_out, "partial_dbrda_anova.csv"))

# Basic dbRDA plot data.
rda_sites <- as.data.frame(vegan::scores(partial_dbrda, display = "sites"))
rda_sites$SampleID <- rownames(rda_sites)
metadata_rda$SampleID <- rownames(metadata_rda)
rda_sites <- dplyr::left_join(rda_sites, metadata_rda, by = "SampleID")

rda_plot <- ggplot(rda_sites, aes(x = CAP1, y = MDS1, colour = Insect)) +
  geom_point(size = 3) +
  scale_colour_manual(values = condition_colours, na.value = "grey70") +
  labs(x = "dbRDA1", y = "dbRDA2", colour = "Treatment") +
  manuscript_theme()

save_pdf(rda_plot, file.path(microbiome_fig, "figure_partial_dbrda.pdf"), width = 14, height = 10)

# -----------------------------------------------------------------------------
# 7. Venn diagrams for ASV turnover
# -----------------------------------------------------------------------------

get_present_asvs <- function(physeq_subset) {
  otu <- as.data.frame(phyloseq::otu_table(physeq_subset))
  if (!phyloseq::taxa_are_rows(physeq_subset)) otu <- t(otu)
  rownames(otu)[rowSums(otu > 0) > 0]
}

venn_groups <- list(
  RH = get_present_asvs(phyloseq::subset_samples(filtered_prevalent, Rep == "RH")),
  RN = get_present_asvs(phyloseq::subset_samples(filtered_prevalent, Rep == "RN")),
  TH = get_present_asvs(phyloseq::subset_samples(filtered_prevalent, Rep == "TH")),
  TN = get_present_asvs(phyloseq::subset_samples(filtered_prevalent, Rep == "TN"))
)

pdf(file.path(microbiome_fig, "figure_venn_week2_week4.pdf"), width = 8, height = 4)
venn_week2 <- VennDiagram::venn.diagram(
  x = list(Herbivory = venn_groups$RH, NoHerb = venn_groups$RN),
  filename = NULL,
  fill = c("#fc8d59", "#91bfdb"),
  alpha = 0.5,
  cex = 1.3,
  cat.cex = 1.1
)
venn_week4 <- VennDiagram::venn.diagram(
  x = list(Herbivory = venn_groups$TH, NoHerb = venn_groups$TN),
  filename = NULL,
  fill = c("tomato", "steelblue"),
  alpha = 0.5,
  cex = 1.3,
  cat.cex = 1.1
)
gridExtra::grid.arrange(venn_week2, venn_week4, ncol = 2)
dev.off()

# -----------------------------------------------------------------------------
# 8. Differential abundance with ANCOM-BC at ASV level
# -----------------------------------------------------------------------------

run_ancombc_asv <- function(physeq, timepoint = "Second") {
  tse <- mia::convertFromPhyloseq(physeq)
  tse_sub <- tse[, tse$Time %in% timepoint]
  tse_sub$Insect <- factor(tse_sub$Insect, levels = c("NoHerb", "Herbivory"))

  ANCOMBC::ancombc(
    data = tse_sub,
    formula = "Insect",
    p_adj_method = "holm",
    prv_cut = 0.10,
    lib_cut = 1000,
    group = "Insect",
    struc_zero = FALSE,
    neg_lb = TRUE,
    tol = 1e-5,
    max_iter = 100,
    conserve = TRUE,
    alpha = 0.05,
    global = TRUE,
    n_cl = 1,
    verbose = TRUE
  )
}

ancombc_second <- run_ancombc_asv(filtered_prevalent, timepoint = "Second")

lfc <- data.frame(ancombc_second$res$lfc[, -1] * ancombc_second$res$diff_abn[, -1], check.names = FALSE) |>
  dplyr::mutate(taxon_id = ancombc_second$res$diff_abn$taxon) |>
  dplyr::select(taxon_id, dplyr::everything())

se <- data.frame(ancombc_second$res$se[, -1] * ancombc_second$res$diff_abn[, -1], check.names = FALSE) |>
  dplyr::mutate(taxon_id = ancombc_second$res$diff_abn$taxon) |>
  dplyr::select(taxon_id, dplyr::everything())
colnames(se)[-1] <- paste0(colnames(se)[-1], "SE")

taxonomy <- SummarizedExperiment::rowData(mia::convertFromPhyloseq(filtered_prevalent)) |> as.data.frame()
taxonomy$taxon_id <- rownames(taxonomy)

da_asv <- lfc |>
  dplyr::left_join(se, by = "taxon_id") |>
  dplyr::left_join(taxonomy, by = "taxon_id") |>
  dplyr::filter(InsectHerbivory != 0) |>
  dplyr::arrange(InsectHerbivory) |>
  dplyr::mutate(
    direction = ifelse(InsectHerbivory > 0, "Enriched under herbivory", "Depleted under herbivory"),
    taxon_label = ifelse(!is.na(Genus), Genus, taxon_id),
    taxon_label = gsub("Burkholderia_Caballeronia_Paraburkholderia", "Burkholderia", taxon_label),
    taxon_label = factor(taxon_label, levels = unique(taxon_label))
  )

write.csv(da_asv, file.path(microbiome_out, "ancombc_second_week_asv_results.csv"), row.names = FALSE)

da_plot <- ggplot(da_asv, aes(x = taxon_label, y = InsectHerbivory, fill = Phylum)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = InsectHerbivory - InsectHerbivorySE, ymax = InsectHerbivory + InsectHerbivorySE), width = 0.2) +
  coord_flip() +
  labs(x = NULL, y = "Log fold change", fill = "Phylum") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

save_pdf(da_plot, file.path(microbiome_fig, "figure_ancombc_asv_second_week.pdf"), width = 16, height = 14)

# -----------------------------------------------------------------------------
# 9. Session information
# -----------------------------------------------------------------------------

write_session_info(file.path(microbiome_out, "session_info_microbiome.txt"))
