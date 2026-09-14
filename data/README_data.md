# Expected data structure

This folder is a placeholder for processed input files used by the R scripts.

Large raw files should not be committed to GitHub. Raw sequencing and metabolomics data should be deposited in appropriate public repositories and linked in the manuscript Data Availability statement.

## Suggested local structure

```text
data/
├── gc_aboveground/
│   ├── metadata.csv
│   ├── Blank.xlsx
│   └── *.txt chromatogram peak tables
├── gc_belowground/
│   ├── metadata.csv
│   ├── Blank.xlsx
│   └── *.txt chromatogram peak tables
├── selected_leaf_vocs/
│   └── 100424_selectedpeaks.csv
 
├── microbiome/
│   ├── ASVs.xlsx
│   ├── TAX_K.xlsx
│   ├── MET.xlsx
│   └── tree.nwk
├── procrustes/
│   ├── pcoa_distances_procrusts.csv
│   ├── pca_distances_procrusts_above.csv
│   ├── pca_distances_procrusts_below_2weeks.csv
│   ├── pca_distances_procrusts_above2weeks.csv
│   ├── pca_distances_procrusts_metabolites.csv
│   └── meta_proc_all.csv
├── roots/
│   └── roots_dryweight.csv
├── ecoplates/
│   ├── 231023_summary.xlsx
│   └── Ecoplates_complete.xlsx
├── rhizosphere_vocs/
│   └── data_norm_VOCs_heat.csv
└── rhizosphere_metabolomics/
    └── data_norm_visual.csv
```

## Notes

- Raw FASTQ and raw mass spectrometry vendor files are deposited outside GitHub.
