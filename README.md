# context-corr

Reproducibility release for **"Context-dependent correlations mislead transcriptomic network inference in bulk and single-cell data"** (Asiaee, Bombina, McGee, Reed, Abrams, Abruzzo, Coombes; *Genome Biology* Brief Report).

A pooled correlation between molecular features mixes within-context and between-context structure. We quantify how often this Simpson-style reversal occurs across TCGA (8,890 tumors, 23M miRNA–mRNA pairs), GTEx tissues, PBMC scRNA-seq, and CITE-seq, and provide a compact reporting checklist that distinguishes genuine within-context association from compositional artifact.

## Layout

```
00-paths.R            # local-path bootstrap, expects ~/Paths/context.json
R/                    # core helpers: contextcorr, sc_context_pipeline, sc_datasets
scripts/              # one driver per paper figure + the four heavy pipelines
data-prep/            # KRC's Rmds that pull TCGA from FireBrowse
reports/              # analytical Rmds with their rendered HTMLs
figures/              # the exact PDFs/PNGs that appear in the paper
results/              # small CSV summaries + the few caches small enough for git
data/README.md        # how to obtain each raw dataset
```

## Quick reproduction — what runs out of the box

These scripts read pre-computed caches included in this release and regenerate the corresponding paper figures with no raw data:

| Figure | Script | Inputs (in `results/`) |
|---|---|---|
| Fig 2C (CITE-seq joint global-vs-resid) | `scripts/gb_fig2C_citeseq_globalvsresid.R` | `cache/sc_200.../pair_stats_Joint_RNA_ADT.rds` |
| Fig 2D (CITE-seq context-definition comparison) | `scripts/gb_fig2D_citeseq_contextdef.R` | `cache/sc_200.../pair_stats_*.rds` |
| Fig S1 (mixed-sign, reversal, I² panels A–D) | `scripts/regen_figS1_styling.R` | `robustness_signflip_summary.csv`, `cache/heterogeneity_I2_FDR.rds` |

From a fresh R session at the repo root:

```r
source("scripts/regen_figS1_styling.R")
source("scripts/gb_fig2C_citeseq_globalvsresid.R")
source("scripts/gb_fig2D_citeseq_contextdef.R")
```

The PDFs are written into `figures/`.

## Full reproduction — figures that need raw data

The remaining figures depend on raw TCGA, GTEx, or single-cell inputs, which are too large to ship:

| Figure | Script | Needs |
|---|---|---|
| Fig 1A, 1B, 1C–F (HTRA3 / miR-122) | `gb_fig1_histogram_and_mir122.R`, `gb_fig1_htra3_panels.R` | TCGA matched mRNA + miRNA |
| Fig 2A (GTEx CEACAM3) | `gb_fig3_panels.R` | GTEx UCSC Xena matrix |
| Fig 2B (PBMC LYZ–FTH1) | `gb_fig3_panels.R` | 10x PBMC 3k |
| Fig S2 (let-7a targets heatmap) | `reports/403-visualization.Rmd` | TCGA matched mRNA + miRNA |
| Fig S3 (BRCA PAM50 / IHC) | `gb_fig4_panels.R` or `figure_brca_subtypes.R` | TCGA BRCA + the heavy BRCA caches (`brca_subtype_PAM50.rds`, `brca_subtype_IHC.rds`) which `scripts/brca_subtype_analysis.R` builds from raw |

To obtain the raw data, see [`data/README.md`](data/README.md). The PDFs as rendered for the paper are already in `figures/` for inspection.

## Heavy pipelines that produce the cached summaries

If you want to regenerate the caches from raw TCGA / single-cell inputs:

```r
# TCGA mixed-sign / Simpson-reversal rates (feeds Fig S1 A,B)
source("scripts/robustness_signflip.R")

# Cochran's Q / I² distribution (feeds Fig S1 C,D)
source("scripts/heterogeneity_I2_FDR.R")

# Pearson vs Spearman robustness (supplementary note)
source("scripts/spearman_sensitivity.R")

# Single-cell pipelines
source("scripts/sc_100_pbmc_coarse_celltypes.R")  # PBMC coarse cell types
source("scripts/sc_110_pbmc_tcell_subtypes.R")    # T-cell-subtype reversal (supp. note)
source("scripts/sc_200_citeseq_protein_rna_contexts.R")  # CITE-seq context-def comparison

# BRCA subtype analysis
source("scripts/brca_subtype_analysis.R")  # builds the heavy BRCA caches
```

Each pipeline caches its intermediate results under `results/cache/` so repeat runs with identical parameters are fast.

## Path configuration

`00-paths.R` reads `$HOME/Paths/context.json`, expected to look like:

```json
{
  "paths": {
    "mirtcga":  "/abs/path/to/your/tcga/clean",
    "scratch":  "/abs/path/to/scratch",
    "results":  "results"
  }
}
```

`mirtcga` is the directory that holds `matchedNoGBM.Rda`, `xena-gtex.Rda`, and `cohortColors.Rda` once you've run the data-prep step.

## Reports

`reports/` carries the analytical narrative as R Markdown plus the rendered HTML. The HTMLs are the canonical record of the numbers as they were used in the paper; the Rmds let you re-run on a fresh dataset.

## Citation

Asiaee A, Bombina P, McGee RL II, Reed J, Abrams ZB, Abruzzo LV, Coombes KR. Context-dependent correlations mislead transcriptomic network inference in bulk and single-cell data. *Genome Biology* (2026).

## License

MIT — see [LICENSE](LICENSE).
