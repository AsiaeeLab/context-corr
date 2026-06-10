# Methods factsheet (verified from code + objects)

This note records *what the code in this repo actually does* for preprocessing and analysis, to keep manuscript Methods text accurate.

## TCGA (FireBrowse) expression inputs

Source code: `ML-TCGA-code/10-fetchRNASEQ.Rmd`, `ML-TCGA-code/20-fetchMirSEQ.Rmd`

- **mRNA-seq (FireBrowse)**: `RSEM_genes_normalized` (RSEM normalized counts, normalized to the third quartile; downloaded per cohort from FireBrowse StandardData).
- **mRNA transform used for correlation**: `log2(10 + ncounts)` (so raw count 0 maps to `log2(10) = 3.321928...`).
- **miRNA-seq (FireBrowse)**: miRSeq RPM (`RPM.txt` files in the archive; downloaded per cohort).
- **miRNA transform used for correlation**: `log2(1 + ncounts)` (RPM of 0 stays 0).

Matched dataset used in analyses:

- File: `/home/amir/datasets/mirTCGA/clean/matchedNoGBM.rda`
- Objects: `mRNAdata` (20289×8890), `miRdata` (1146×8890), `cancerType` (31 cohorts), `sampleType`, `cohortColors`.

Practical consequence of the mRNA transform:

- Many mRNA values are exactly `log2(10)` (floor corresponding to zero raw count), so “undetectable in a cohort” manifests as **zero variance** within that cohort (all samples at the floor), yielding an undefined within-cohort correlation that we treat as missing.

## GTEx pilot (UCSC Xena)

- File: `/home/amir/datasets/mirTCGA/clean/xena-gtex.Rda`
- Object: `gtex$mrna` (19069×7792 data.frame)
- Units: verified by range/minimum to be **`log2(TPM + 0.001)`** (minimum ~= `log2(0.001) = -9.966`).
- Tissue context label: parsed from Xena sample IDs as a 4-digit code (e.g. `GTEX.*.0011.*`).

## New “major revision” robustness analyses (this branch)

All scripts are runnable from repo root and write outputs under `results/` and `doc/bioinformatics/figures/`:

- `scripts/robustness_signflip.R`
  - Effect-size–conditioned mixed-sign and Simpson reversal rates (eps grid + thresholds).
  - Outputs: `results/robustness_signflip_summary.csv`, `doc/bioinformatics/figures/Figure_RobustnessSignFlip.pdf`.
- `scripts/heterogeneity_I2_FDR.R`
  - Cochran’s Q / I² heterogeneity plus BH FDR within the tested domain.
  - Outputs: `results/heterogeneity_summary.csv`, `doc/bioinformatics/figures/Figure_Heterogeneity_I2.pdf`.
- `scripts/spearman_sensitivity.R`
  - Pearson vs Spearman concordance and Simpson reversal rate under Spearman on a subset.
  - Outputs: `results/spearman_sensitivity.csv`, `doc/bioinformatics/figures/Figure_PearsonVsSpearman.pdf`.
- `scripts/figure_global_distribution.R`
  - Main Figure 1: sampled global correlation distribution + miR-122 tissue-driven examples.
  - Output: `doc/bioinformatics/figures/Figure_GlobalDistribution.pdf`.

Minimal reusable API for context-aware correlation summaries:

- `R/contextcorr.R`

