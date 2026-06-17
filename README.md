# Aphid herbivory and wheat rhizosphere microbiome analyses

This repository contains R scripts associated with the manuscript:

**Aphid herbivory transiently restructures rhizosphere chemistry and bacterial communities in wheat**

The scripts document the analysis and figure-generation workflows used for volatile organic compounds, rhizosphere metabolomics, microbiome composition, Procrustes integration analyses, EcoPlate substrate utilisation, and root dry weight.

## Repository structure

```text
.
├── README.md
├── data/
│   └── README_data.md
├── scripts/
│   ├── 00_project_setup.R
│   ├── 01_gc_voc_alignment_and_pca.R
│   ├── 02_leaf_voc_heatmap_selected_features.R
│   ├── 03_microbiome_analysis.R
│   ├── 04_procrustes_analysis.R
│   ├── 05_root_dry_weight.R
│   ├── 06_ecoplate_analysis.R
│   └── 07_rhizosphere_metabolite_heatmaps.R
├── figures/
└── outputs/
```

## Scripts

### `00_project_setup.R`
Shared paths, colours, plotting theme, and utility functions.

### `01_gc_voc_alignment_and_pca.R`
GC peak alignment and normalisation using `GCalignR`, followed by PCA and PERMANOVA for leaf and rhizosphere VOC datasets.

### `02_leaf_voc_heatmap_selected_features.R`
Heatmap of selected herbivory-associated leaf VOC features, including feature annotations and colour keys.

### `03_microbiome_analysis.R`
Microbiome workflow including phyloseq import, removal of chloroplast/mitochondrial ASVs, rarefaction curves, alpha diversity, Bray-Curtis NMDS, PERMANOVA, dbRDA, Venn diagrams, and ANCOM-BC differential abundance.

### `04_procrustes_analysis.R`
Procrustes analyses linking leaf VOC profiles with rhizosphere microbiome composition and rhizosphere metabolite profiles. Also generates the bootstrap Procrustes histogram.

### `05_root_dry_weight.R`
Wilcoxon tests and supplementary root dry weight figure.

### `06_ecoplate_analysis.R`
Biolog EcoPlate AWCD analysis, treatment x time model, post hoc grouping, substrate-level tests, and Figure 6 generation.

### `07_rhizosphere_metabolite_heatmaps.R`
Supplementary heatmaps for rhizosphere VOCs and LC-MS non-volatile metabolite features.

## Data availability


- Raw 16S rRNA sequencing reads: NCBI Sequence Read Archive under BioProject accession number PRJNA1390630.

- Raw LC-MS/GC-MS metabolomics files: 


See `data/README_data.md` for the expected local file structure.

## Reproducibility notes

These scripts were tidied from exploratory analysis scripts. They are intended to document the final workflow used in the manuscript rather than serve as a general-purpose R package.

Each analysis script writes a `session_info_*.txt` file to the `outputs/` folder to document R and package versions.

Several scripts contain `TODO` comments where final file names, sample orders, or repository DOI links should be confirmed before public release.

## Main R packages

- `GCalignR`
- `phyloseq`
- `vegan`
- `mia`
- `miaViz`
- `ANCOMBC`
- `FactoMineR`
- `factoextra`
- `ComplexHeatmap`
- `ggplot2`
- `rstatix`
- `ggpubr`
- `agricolae`
- `patchwork`

## Citation

#TODO ADD ZENODO DOI IF NEEDED

