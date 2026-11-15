# context-corr

Reproducibility release for **"Context-dependent correlation reversals in transcriptomic network inference across bulk and single-cell data"** (Asiaee, Bombina, McGee, Reed, Abrams, Abruzzo, Coombes; *BMC Genomics*, in revision).

A pooled correlation between molecular features mixes within-context and between-context structure. We quantify how often pooled and context-conditional directions disagree across TCGA (8,890 tumors, 23M miRNA–mRNA pairs), GTEx tissues, PBMC scRNA-seq, and CITE-seq, and provide a compact checklist for reporting the pooled coefficient alongside within-context associations and context-mean structure.

## Layout

```
00-paths.R            # dataset paths: CONTEXT_PATHS or ~/Paths/context.json
R/                    # core helpers: contextcorr, sc_context_pipeline, sc_datasets
examples/             # worked_example.R (the reporting checklist on simulated data), benchmark_runtime.R
tests/                # run_tests.R, plain-R checks of the reporting functions
scripts/              # figure drivers and original/revision analysis pipelines
data-prep/            # KRC's Rmds that pull TCGA from FireBrowse
reports/              # analytical Rmds with their rendered HTMLs
figures/              # the exact PDFs/PNGs that appear in the paper
results/              # processed result tables and selected reusable summaries
data/README.md        # how to obtain each raw dataset
```

## Installation

This repository is not an R package (there is no `DESCRIPTION` file), so `remotes::install_github()` and `devtools::load_all()` do not apply. The reporting interface is one base-R file, `R/contextcorr.R`, with no dependencies beyond `stats`. Either clone and `source()` it:

```sh
git clone https://github.com/AsiaeeLab/context-corr.git
cd context-corr
```

```r
source("R/contextcorr.R")
```

or source it directly from GitHub without cloning:

```r
source("https://raw.githubusercontent.com/AsiaeeLab/context-corr/main/R/contextcorr.R")
```

Tested on R 4.6.1. The figure scripts and single-cell pipelines have their own dependencies (Seurat, ggplot2, rjson, and others, listed at the top of each script); the reporting functions need none of them.

## Reporting checklist interface (`R/contextcorr.R`)

Definitions, matching the paper:

- **Sign tolerance** `eps = 0.05`: a context votes +1 if `r_c >= eps`, -1 if `r_c <= -eps`, 0 otherwise. A pair is **mixed-sign** if both votes occur.
- **Simpson reversal**: `sign(r_global)` differs from the majority context vote (`simpson_flag`, `majority = "count"`).
- **Residualized correlation**: correlation after subtracting per-context means from both variables.
- **Fisher-z interval** for `r_c`: `z_c +/- qnorm(1 - alpha/2) / sqrt(n_c - 3)` mapped back to the r scale (the 1.96 interval at `alpha = 0.05`). Contexts with `n_c < 4` have no interval.
- **Heterogeneity**: Cochran's Q on Fisher z with weights `w = n_c - 3`, `I2 = max(0, (Q - df) / Q)`.
- **Random-effects pooling**: DerSimonian-Laird on Fisher z, `tau2 = max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))`, `w_RE = 1 / (1/w + tau2)`, `r_RE = tanh(sum(w_RE z) / sum(w_RE))`.

| Function | Returns |
|---|---|
| `context_correlations(expr_x, expr_y, context, method = c("pearson", "spearman"))` | `r_global`, `r_by_context`, `n_by_context`, `k`, `df` |
| `heterogeneity_Q_I2(r_by_context, n_by_context)` | `k`, `df`, `Q`, `p`, `I2` |
| `random_effects_pool(r_by_context, n_by_context)` | `k`, `tau2`, `z_RE`, `r_RE` |
| `mean_residualized_cor(expr_x, expr_y, context, method = c("pearson", "spearman"))` | the residualized correlation |
| `simpson_flag(r_global, r_by_context, eps = 0.05, majority = c("count", "median"))` | `flag`, signs, vote counts |
| `contextcorr_fisher_z_ci(r, n, alpha = 0.05)` | data.frame: `r`, `n`, `z`, `se`, `ci_lo`, `ci_hi`, `p`, `excludes_zero` |
| `classify_pair(r_global, r_by_context, n_by_context, eps = 0.05, switch_thresh = 0.3, alpha = 0.05)` | list: `label`, `reason`, per-context table, counts, `simpson` |
| `residualize_matrix(expr, context)` | features x samples matrix with per-context means removed |
| `random_effects_adjacency(expr, context, features = NULL, min_n = 4, method = c("pearson", "spearman"))` | list: `r_RE` (p x p, diagonal 1), `I2`, `tau2`, `k`, `n_by_context` |
| `context_report(expr_x, expr_y, context, eps = 0.05, method = c("pearson", "spearman"), min_n = 4, switch_thresh = 0.3, alpha = 0.05)` | the five checklist items plus `simpson` and `classification` |
| `print_context_report(report, digits = 3)` | console summary of a `context_report()` result |

`context_report()` returns, for one pair: (1) `r_global`; (2) `r_by_context` and `n_by_context`; (3) `heterogeneity` (`Q`, `df`, `p`, `I2`); (4) `random_effects` (`r_RE`, `tau2`); (5) `r_resid`; plus `simpson = simpson_flag(...)` and `classification = classify_pair(...)`. Contexts with fewer than `min_n` samples keep their `n` but have `r` set to `NA` and are excluded from every pooled or conditional estimand, mirroring the paper's retained-context analysis population; `n` is the resulting pooled sample size. A pair-specific perfect correlation (`|r_c| = 1`) is also returned as `NA` because its Fisher transform is infinite, but its observations remain in the retained pooled and residualized population.

`residualize_matrix()` and `random_effects_adjacency()$r_RE` are the two drop-in inputs for network tools: pass the residualized matrix to any correlation-based method (WGCNA, ARACNe, GENIE3), or use `r_RE` in place of a pooled correlation matrix. `random_effects_adjacency()` computes one `cor()` per context on the whole feature block, then pools each of the `p(p-1)/2` pairs with `random_effects_pool()`, so it is O(p^2 k) in time and stores a `k x p(p-1)/2` matrix of within-context correlations.

### Pair classification

`classify_pair()` uses only informative contexts (finite `|r_c| < 1`, `n_c >= 4`), requires at least two informative contexts for a substantive label, and tests the labels in this order, the first match winning:

1. `context_switching`: at least one context with `r_c >= +switch_thresh` and another with `r_c <= -switch_thresh`, both with Fisher-z `p < alpha`.
2. `compositional`: Simpson reversal by majority vote (pooled sign disagrees with the eps-tolerant majority of context votes).
3. `globally_robust` (the machine-readable token for “globally directionally robust” in the manuscript): `|r_global| >= eps`, at least one context has a Fisher-z CI excluding 0, every such context has the pooled sign, and no reversal.
4. `indeterminate`: none of the above (no informative context, no context individually distinguishable from 0, `|r_global| < eps`, or a significant dissenting context that does not reach `switch_thresh`); the `reason` element says which.

`compositional` is an operational token for pooled-versus-majority disagreement. It does not by itself establish that the between-context covariance term caused the disagreement; use `r_resid` and the covariance decomposition for that mechanism check.

Precedence: `context_switching` is checked before `compositional` because a pair with strong opposite-sign contexts often also reverses the pooled sign relative to the vote majority, and the switching label is the more specific description. A switching pair can never be `globally_robust`, since it has significant contexts of both signs; `compositional` precedes `globally_robust` because the latter requires no reversal.

## Worked example

`examples/worked_example.R` (run with `Rscript examples/worked_example.R` from the repo root) simulates 4 contexts x 150 samples and one pair of each type. The generator draws `x = mu_x[c] + e1` and `y = mu_y[c] + r[c] e1 + sqrt(1 - r[c]^2) e2` per context:

```r
source("R/contextcorr.R")
set.seed(20260903)
context <- factor(rep(paste0("C", 1:4), each = 150))

sim_pair <- function(r_within, mu_x, mu_y, context) {
  x <- y <- numeric(length(context))
  for (i in seq_along(levels(context))) {
    idx <- which(context == levels(context)[i])
    e1 <- rnorm(length(idx)); e2 <- rnorm(length(idx))
    x[idx] <- mu_x[i] + e1
    y[idx] <- mu_y[i] + r_within[i] * e1 + sqrt(1 - r_within[i]^2) * e2
  }
  list(x = x, y = y)
}

# r = -0.4 inside every context, but both means rise together across contexts
comp <- sim_pair(rep(-0.4, 4), mu_x = c(0, 1.5, 3, 4.5), mu_y = c(0, 1.5, 3, 4.5), context)
# r = +0.5 everywhere, no mean shifts
robust <- sim_pair(rep(0.5, 4), rep(0, 4), rep(0, 4), context)
# r = +0.5 in two contexts, -0.5 in the other two
switching <- sim_pair(c(0.5, 0.5, -0.5, -0.5), rep(0, 4), rep(0, 4), context)

rpt <- context_report(comp$x, comp$y, context, eps = 0.05, method = "pearson")
print_context_report(rpt)
```

Output for the built-in reversal:

```
Context-aware correlation report (pearson; eps = 0.05, min_n = 4)
1. Pooled correlation: r_global = 0.638 (n = 600)
2. Within-context correlations (k = 4 contexts used):
 context   n      r  ci_lo  ci_hi       p vote
      C1 150 -0.363 -0.494 -0.215 4.0e-06   -1
      C2 150 -0.390 -0.518 -0.245 6.1e-07   -1
      C3 150 -0.426 -0.549 -0.285 3.4e-08   -1
      C4 150 -0.456 -0.574 -0.319 2.3e-09   -1
3. Heterogeneity: Q = 1.070 on 3 df, p = 0.78, I2 = 0.0%
4. Random-effects pooled correlation: r_RE = -0.409 (tau2 = 0.000, k = 4)
5. Mean-residualized correlation: r_resid = -0.409
Simpson reversal flag: TRUE (pooled sign 1, majority sign -1; votes +0 / -4)
Classification: compositional
  pooled sign (+) disagrees with the majority of context votes (+0 / -4 at eps = 0.05)
```

The other pairs come out as `globally_robust` (r_global = 0.479, every context 0.39 to 0.51, I2 = 0%), `context_switching` (r_global = -0.015, contexts +0.60, +0.49, -0.57, -0.53, I2 = 98.7%, tau2 = 0.50), and a pure-noise pair as `indeterminate` ("no context has a Fisher-z CI excluding 0").

The script then builds a 20 x 600 toy matrix in which genes g1 to g5 share a latent factor while their context means move in opposite directions (g1 to g3 up, g4 and g5 down), and compares three adjacency choices:

```r
res <- residualize_matrix(expr, context)          # per-context means now ~1e-16
adj <- random_effects_adjacency(expr, context)    # adj$r_RE is p x p, symmetric, diagonal 1
cor(t(expr))["g1", "g4"]   # pooled:              -0.57
adj$r_RE["g1", "g4"]       # random-effects r_RE: +0.65
cor(t(res))["g1", "g4"]    # residualized:        +0.65
adj$I2["g1", "g4"]         # 0 (the four contexts agree)
```

`context_report(expr["g1", ], expr["g4", ], context)` labels the pair `compositional`.

## Tests

```sh
Rscript tests/run_tests.R
```

Plain R with `stopifnot()`; no testthat needed. The 19 checks cover: the four toy pairs receive their intended labels and `classify_pair()` precedence and edge cases behave as documented; `random_effects_pool()` reproduces the hand-coded DerSimonian-Laird formulas (heterogeneous case with `tau2 > 0`, homogeneous case with `tau2 = 0`, exclusion of `NA`, `n < 4`, and `|r| = 1` contexts) and `heterogeneity_Q_I2()` reproduces Q, df, p, I2; `contextcorr_fisher_z_ci()` matches `tanh(atanh(r) +/- 1.96 / sqrt(n - 3))`; `residualize_matrix()` gives per-context means below 1e-12, keeps dimnames, agrees with `mean_residualized_cor()`, drops NA-context samples, and zeroes singleton contexts; `random_effects_adjacency()` is symmetric with unit diagonal, has a symmetric `I2`, honors `min_n` and `features`, and every entry equals `random_effects_pool()` on the directly computed per-context correlations; `context_report()` is internally consistent with the component functions. The script exits with status 1 on any failure.

## Runtime

`examples/benchmark_runtime.R` times the functions on synthetic data at the paper's scale (31 contexts, 8,890 samples, context sizes 37 to 1,302). R 4.6.1 on x86_64 Linux (Intel Xeon w9-3575X), single-threaded:

| Call | Elapsed |
|---|---|
| `context_report()` for one pair, first call | 0.10 s |
| `context_report()` for 100 pairs in a loop | 0.55 s (5.5 ms per pair) |
| `residualize_matrix()` on 200 features x 8,890 samples | 0.03 s |
| `random_effects_adjacency()` on 200 features x 8,890 samples (19,900 pairs) | 1.0 s |

`random_effects_adjacency()` scales as `p^2`: about 200 pairs per millisecond, so 2,000 features (2M pairs) take on the order of two minutes.

## Reproducing the revised analyses

Run commands from the repository root. The reporting functions, tests and worked
example require only base R. Large-data analyses require the input files described
in [data/README.md](data/README.md) and additional packages listed in their scripts.
They do not run as part of the tests.

The current results are in `results/revision/`; [its README](results/revision/README.md)
identifies the result tables for each analysis. The earlier HTML reports in
`reports/` retain the original analyses and are not the source of the revised
population-specific rates or current figure numbering.

| Analysis | Script | Main processed outputs |
|---|---|---|
| Common reporting and preprocessing tables | `scripts/revision/uniform_table.R` | `uniform_table/table1_metrics.csv`, `table2_preprocessing.csv`, `bootstrap_cis_by_setting.csv` |
| TCGA voting, FDR, bootstrap and measurability checks | `scripts/revision/tcga_robustness.R` | `tcga_robustness/reversal_rules.csv`, `fdr_gated_by_cohort.csv`, `bootstrap_cis.csv` |
| Null calibration and downsampling | `scripts/revision/null_simulation.R` | `null_simulation/null_rates.csv`, `downsampling_tcga.csv` |
| Fixed-population T-cell subdivision and sparse-data sensitivity | `scripts/revision/singlecell_tcell.R` | `singlecell_tcell/tcell_controlled.csv`, `sparse_check.csv` |
| Eight-donor PBMC and pseudobulk analysis | `scripts/revision/multidonor.R` | `multidonor/cell_level_rates.csv`, `per_donor_rates.csv`, `pseudobulk_rates.csv` |
| Cell-composition adjustment | `scripts/revision/composition.R` | `composition/adjusted_rates.csv`, `reversal_transitions.csv` |
| WGCNA/GENIE3 comparisons and exploratory STRING/TRRUST benchmarks | `scripts/revision/downstream.R` | `downstream/edge_overlap.csv`, `module_agreement.csv`, `benchmark_precision.csv` |
| TCGA sample-type sensitivity and Fisher-$z$ variance checks | `scripts/revision/tcga_checks.R`, `tcga_skcm_participants.R`, `null_tail.R` | `sampletype_fisher_checks/tcga_sampletype_metrics.csv`, `fisher_B2_permutation_heterogeneity.csv`, `fisher_B3_bootstrap_variance_ratio.csv`, `fisher_B4_spearman_heterogeneity.csv` |

### Small reproductions using included outputs

These commands use processed results included in this repository and do not need
raw expression data or path configuration:

```sh
Rscript tests/run_tests.R
Rscript examples/worked_example.R
Rscript scripts/bmc_fig3_panels.R
Rscript scripts/revision/regenerate_downstream_figure.R
Rscript scripts/revision/regenerate_multidonor_figure.R
```

The first figure command writes the four TCGA panels into `figures/`; the last two
write their supplementary figures into the corresponding `results/revision/`
subdirectory. Those supplementary figures also appear in `figures/` as the copies
used by the manuscript.

### Full analyses

Configure paths first. The TCGA analyses require `matchedNoGBM.rda` (or `.Rda` where
supported), while the single-cell analyses require the exact Seurat inputs.
Run the original preparation scripts before the analyses that consume their caches:

```sh
Rscript scripts/robustness_signflip.R
Rscript scripts/heterogeneity_I2_FDR.R
Rscript scripts/brca_subtype_analysis.R
Rscript scripts/sc_100_pbmc_coarse_celltypes.R
Rscript scripts/sc_110_pbmc_tcell_subtypes.R
Rscript scripts/sc_200_citeseq_protein_rna_contexts.R
Rscript scripts/revision/uniform_table.R
Rscript scripts/revision/tcga_robustness.R
Rscript scripts/revision/null_simulation.R
Rscript scripts/revision/singlecell_tcell.R
Rscript scripts/revision/multidonor.R
Rscript scripts/revision/composition.R
Rscript scripts/revision/downstream.R
Rscript scripts/revision/tcga_checks.R
Rscript scripts/revision/tcga_skcm_participants.R
Rscript scripts/spearman_sensitivity.R
```

The single-cell pipelines use Seurat/SeuratObject and Matrix. The downstream
comparison also requires data.table, WGCNA, dynamicTreeCut, GENIE3, doRNG and the
parallel backends checked by that script. Composition adjustment requires the
`estimate` package. Figure scripts may require ggplot2. See `sessionInfo.txt` and
`runtime.txt` alongside the saved results for the environments used in the study.
Several large-data scripts use forked parallel workers and were run on Linux.
Rerunning a full analysis overwrites its processed outputs; use a separate clone
if you want to preserve the supplied results. Quick modes are computational checks,
not reproductions of the full-data tables.

For the updated GTEx and single-cell scatter panels, run:

```sh
Rscript scripts/revision/fig_gtex_resid.R
Rscript scripts/revision/regenerate_uniform_singlecell_figures.R
```

These require raw inputs and the single-cell preprocessing outputs above and write
to `figures/`. Current main-figure numbering is: 1, pooled distribution and miR-122;
2, HTRA3 example; 3, TCGA heterogeneity and validated targets; 4, GTEx; 5, PBMC and
CITE-seq; 6, BRCA partitions. Some filenames retain their original numbering.

## Path configuration

Copy `data/context-paths.example.json` to a local JSON file and set its `mirtcga`,
`scratch`, and `results` paths. Keep `results` equal to `results` when following the
commands above because the revision scripts read and write repository-relative
result directories. Install `rjson` for the configuration reader.

```sh
CONTEXT_PATHS=/absolute/path/to/context-paths.json Rscript scripts/revision/tcga_robustness.R
```

If `CONTEXT_PATHS` is unset, `00-paths.R` reads `~/Paths/context.json` for compatibility
with the original scripts. The reporting functions, tests, and worked example do
not require a configuration file. Put the GTEx v8 sample attributes at
`<mirtcga>/GTEx_SampleAttributes_v8.txt`, or set `GTEX_ATTRIBUTES` to that file.

## Citation

Asiaee A, Bombina P, McGee RL II, Reed J, Abrams ZB, Abruzzo LV, Coombes KR.
Context-dependent correlation reversals in transcriptomic network inference across
bulk and single-cell data. Manuscript under revision at *BMC Genomics* (2026).
See `CITATION.cff` for machine-readable citation information.

## License

MIT; see [LICENSE](LICENSE).
