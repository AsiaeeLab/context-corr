# Processed results for the revised manuscript

These tables contain the completed analyses reported in the revised manuscript.
The analysis scripts are in `scripts/revision/`; seeds, thresholds, filtering and
bootstrap settings are specified there. Fractions and percentages follow their
column names; do not interchange them when comparing tables.

| Directory | Contents |
|---|---|
| `uniform_table/` | The common reporting table, preprocessing counts, retained context sizes, and pair-bootstrap intervals across six settings. `table1_metrics.csv` and `table2_preprocessing.csv` supply the main manuscript tables. |
| `tcga_robustness/` | Cohort sizes, alternative voting rules, per-cohort FDR thresholds and reversal survival, bootstrap intervals/stability, and undefined-correlation checks. |
| `null_simulation/` | Balanced and unequal-context null scenarios, empirical mean-shift summaries, and TCGA downsampling results. |
| `singlecell_tcell/` | Fixed-population subdivision results, within-context aggregation and rank-correlation sensitivity, and input metadata. `<scratch>` in `pbmc_donor_check.csv` denotes the configured scratch directory; the original local path has been generalized. |
| `multidonor/` | Cell-level, per-donor and pseudobulk results from GSE96583 batch 2, including retained-cell counts and pairwise summaries needed to reproduce Figure S5. |
| `composition/` | Public TCGA composition covariates, adjusted correlation summaries, reversal transitions, and mean-shift attenuation. |
| `downstream/` | Network ranking overlap, module agreement, pair classifications, and exploratory STRING/TRRUST reference comparisons. |
| `sampletype_fisher_checks/` | TCGA metrics recomputed on primary tumors only (and other sample-type variants), per-cohort stability of within-cohort correlations, and checks of the Fisher-$z$ variance by within-cohort permutation, within-cohort bootstrap, and Spearman correlation. |

`reproduction_checks.csv`, `gtex_convention_grid.csv` and
`headline_reproduction.csv` retain explicit comparisons between alternative
analysis conventions. Their check columns are not substitutes for the final
common-population estimates in `uniform_table/table1_metrics.csv`.

`multidonor/cell_pair_summary.rds` and `pseudobulk_pair_summary.rds` are the
`summary_pairs` tables extracted without numerical changes from the full saved
cell-type and donor-by-cell-type analyses. They allow the supplementary figure
to be regenerated without shipping expression matrices or intermediate caches.
After rerunning `multidonor.R`, these summaries are also exported by that script.

The included `sessionInfo.txt` and `runtime.txt` files describe the original runs;
they are provenance records, not claims that all full analyses were rerun for the
repository update. The unit tests and worked example are independent of those
large-data runs. Raw downloads and regenerable intermediate objects are excluded.
