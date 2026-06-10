# Obtaining the raw data

None of the raw datasets are shipped with this release. Each table below points to its public source and tells you where to place the file so the scripts find it (paths are resolved through `~/Paths/context.json`).

## TCGA matched miRNA + mRNA (the heaviest input)

The paper uses 8,890 samples across 31 cancer cohorts from FireBrowse standardized data `stddata__2016_01_28` (mRNA: upper-quartile-normalized RSEM counts; miRNA: RPM).

The cleanest path is to run KRC's data-prep Rmds in order:

```r
# Working dir: data-prep/
rmarkdown::render("01-cohorts.Rmd")
rmarkdown::render("10-fetchRNASEQ.Rmd")
rmarkdown::render("20-fetchMirSEQ.Rmd")
```

These download per-cohort archives from FireBrowse, parse them into per-gene/per-miRNA matrices, and finally build the matched object `matchedNoGBM.Rda` (objects: `mRNAdata` 20289×8890, `miRdata` 1146×8890, `cancerType`, `sampleType`, `cohortColors`).

Place the resulting `matchedNoGBM.Rda` (and `cohortColors.Rda`) under whatever directory you point `paths$mirtcga` at in `~/Paths/context.json`.

Expected disk: ~5 GB for the FireBrowse archives during fetch; the matched object is much smaller.

## GTEx (UCSC Xena)

The Xena "TCGA TARGET GTEx" expression matrix (values in `log2(TPM + 0.001)`) is at:

- <https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count>

Tissue-of-origin labels come from the GTEx v8 sample annotations (the `SMTSD` field, joined on sample identifier). Save the resulting object as `xena-gtex.Rda` (containing `gtex$mrna`, a 19069×7792 data frame) under `paths$mirtcga`.

## PBMC 3k scRNA-seq

Auto-downloaded from `SeuratData`. No manual step:

```r
install.packages("SeuratData", repos = "https://satijalab.r-universe.dev")
SeuratData::InstallData("pbmc3k")
```

Then `R/sc_datasets.R::load_pbmc_scrna()` will find it.

## CITE-seq (PBMC, RNA + ADT)

The paper analyses 7,798 cells after QC from a PBMC CITE-seq dataset (Stoeckius et al., 2017). The release loader (`R/sc_datasets.R::load_pbmc_citeseq`) looks first for a local Seurat object at `<scratch>/sc_inputs/pbmc_citeseq.rds`, then falls back to the `SeuratData` candidates `bmcite` and `pbmcmultimodal`.

If you do not already have the exact dataset, the simplest fallback is:

```r
SeuratData::InstallData("bmcite")     # or
SeuratData::InstallData("pbmcmultimodal")
```

These are not numerically identical to the paper's input, but they reproduce the qualitative protein–RNA-context pattern.

## miRTarBase v10

Used for the validated-target heatmap (Fig S2):

- <https://mirtarbase.cuhk.edu.cn/>

Download the v10 release; one Excel/CSV per organism.
