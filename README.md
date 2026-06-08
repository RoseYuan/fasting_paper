# Effects of long-term fasting: longitudinal epigenetic responses in humans

This repository contains the analysis code for a longitudinal fasting study
combining clinical measurements, Illumina EPIC DNA methylation arrays,
differential methylation analysis, and epigenetic clock response analysis.

The workflow is written mainly in R and Quarto. It builds harmonized clinical
metadata, preprocesses methylation beta values, runs differential methylation
models across fasting timepoints, predicts multiple epigenetic clocks, and
generates manuscript figures and supplementary tables.

## Repository layout

| Path | Purpose |
| --- | --- |
| `0_merge_meta.R` | Merge sample-sheet, daily clinical metadata, three-timepoint metadata, and AGE/glycation data into project metadata tables. |
| `1_processing.qmd` | Process raw Illumina EPIC IDAT files with Meffil and Sesame; export QC summaries and beta matrices. |
| `1_visualization_clinical_params.qmd` | Plot longitudinal weight, BMI, and selected clinical parameters; export Figure 2 components and clinical summary tables. |
| `1_visualization_pca.qmd` | Run PCA on baseline clinical parameters and clinical change scores; assemble Figure 2 and supplementary PCA figures. |
| `1_pca_white_blood_cells.qmd` | Supplementary PCA analysis of white blood cell parameters. |
| `2_qc_and_prepare_filtered_betas.qmd` | Explore, normalize, filter, and batch-correct stringent beta matrices for differential methylation and classic clocks. |
| `2_prepare_complete_betas_for_clocks.qmd` | Prepare a less stringently filtered, imputed, normalized, batch-corrected beta matrix for newer epigenetic clocks. |
| `3_DM.qmd` | Longitudinal differential methylation analysis with limma, enrichment analysis, DMRcate, and Figure 3 outputs. |
| `3_ClassicClocks.R` | Predict Horvath, Hannum, SkinBlood, PhenoAge, and DunedinPACE values. |
| `3_PCClocks.R` | Predict principal-component epigenetic clocks and acceleration outputs. |
| `3_GrimAge.R` | Calculate DNAm GrimAge and GrimAge v2 outputs and components. |
| `3_OtherClocks.R` | Calculate SystemsAge, AdaptAge, DamAge, and CausAge. |
| `4_merge_clocks.R` | Merge all clock outputs into one sample-level clock table. |
| `5_clock_qc_and_response_analysis.qmd` | QC clock predictions and screen residual age acceleration associations with maximum weight loss. |
| `5_integrated_clock_analysis.qmd` | Integrated clock-response analysis, SRAA trajectories, clock PCA, clinical associations, and Figure 4 outputs. |
| `utils.R` | Shared metadata, plotting, modeling, duplicate-column, PCA, and clock-analysis helper functions. |
| `data/` | Generated or local analysis data, including merged metadata, beta matrices, model fits, and clock predictions. Ignored by Git. |
| `figures/` | Generated manuscript and supplementary figures/tables. Ignored by Git. |

## Analysis overview

The project has three main analysis branches:

1. Metadata and clinical visualization
   - `0_merge_meta.R` creates `data/meta_3timepoints.csv` and
     `data/meta_daily.csv`.
   - `1_visualization_clinical_params.qmd`,
     `1_visualization_pca.qmd`, and `1_pca_white_blood_cells.qmd`
     summarize clinical trajectories and clinical PCA structure.

2. Methylation preprocessing and differential methylation
   - `1_processing.qmd` performs raw IDAT QC and Sesame preprocessing.
   - `2_qc_and_prepare_filtered_betas.qmd` creates stringent beta matrices:
     `data/cleaned_beta_s.RDS`, `data/cleaned_beta_s_corrected.RDS`,
     `data/cleaned_beta_p.RDS`, and `data/cleaned_beta_p_corrected.RDS`.
   - `3_DM.qmd` uses uncorrected filtered betas and models scan date as a
     covariate in longitudinal limma models.

3. Epigenetic clocks and integrated response analysis
   - `2_prepare_complete_betas_for_clocks.qmd` creates
     `data/betas_complete_processed.RDS` for newer clock models that need
     broader CpG coverage.
   - `3_ClassicClocks.R`, `3_PCClocks.R`, `3_GrimAge.R`, and
     `3_OtherClocks.R` generate clock-specific output tables.
   - `4_merge_clocks.R` creates `data/AllClocks_Merged.csv`.
   - `5_clock_qc_and_response_analysis.qmd` and
     `5_integrated_clock_analysis.qmd` analyze clock trajectories, residual age
     acceleration, standardized residual age acceleration, weight-loss
     associations, clinical associations, and Figure 4 outputs.

## Notes on beta matrix usage

Two beta-processing strategies are used intentionally:

- Differential methylation and classic clock analyses use the more stringent
  filtered matrices from `2_qc_and_prepare_filtered_betas.qmd`.
- Newer clock analyses use the complete processed matrix from
  `2_prepare_complete_betas_for_clocks.qmd`, with missing values imputed to
  preserve clock CpG coverage.

For differential methylation, `3_DM.qmd` uses uncorrected filtered beta values
and controls scan-date batch effects inside the statistical design rather than
using pre-corrected beta values.

## License

This project is licensed under the GNU General Public License v3.0. See
`LICENSE` for details.
