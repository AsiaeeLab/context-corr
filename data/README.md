# Obtaining the raw data

Raw expression datasets are not shipped with this release. Paths are resolved through `CONTEXT_PATHS` or `~/Paths/context.json`; see `context-paths.example.json`. Saved result tables and selected pairwise summaries are included under `results/`.

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

Tissue-of-origin labels come from the GTEx v8 sample annotations (the `SMTSD` field, joined on sample identifier). Save the expression object as `xena-gtex.Rda` (containing `gtex$mrna`, a 19069×7792 data frame) under `paths$mirtcga`. Save `GTEx_SampleAttributes_v8.txt` in the same directory, or set `GTEX_ATTRIBUTES` to its full path. The revised analysis uses the 20 retained SMTSD tissues for both feature selection and correlation estimation.

## PBMC 3k scRNA-seq

The paper's input is the 10x Genomics PBMC 3k dataset. The loader accepts a Seurat object at `<scratch>/sc_inputs/pbmc_scrna.rds`. It can also use an installed SeuratData PBMC 3k object:

```r
install.packages("SeuratData", repos = "https://satijalab.r-universe.dev")
SeuratData::InstallData("pbmc3k")
```

The analysis scripts call `load_pbmc_scrna(auto_download = FALSE)`, so install the dataset or prepare the local object before running them. For exact numerical reproduction, retain the paper's raw PBMC 3k cells and apply the quality-control settings in the scripts.

## CITE-seq (PBMC, RNA + ADT)

The paper uses the **10x Genomics PBMC 10k with TotalSeq-B antibodies**, dataset `pbmc_10k_protein_v3`, with 7,798 cells after QC. Stoeckius et al. (2017) is the reference for the CITE-seq method, not the accession for this input. Prepare a Seurat object containing the RNA and ADT count assays and save it at `<scratch>/sc_inputs/pbmc_citeseq.rds`. The loader looks there first, then checks installed SeuratData candidates `bmcite` and `pbmcmultimodal`.

For testing the pipeline on alternative data only:

```r
SeuratData::InstallData("bmcite")     # or
SeuratData::InstallData("pbmcmultimodal")
```

These alternatives are not the paper's input and must not be used to claim reproduction of its numerical results.

## Multi-donor PBMC

`scripts/revision/multidonor.R` downloads GEO **GSE96583 batch 2**, control GSM2560248 and IFN-beta-stimulated GSM2560249, plus the accompanying donor/cell-type metadata. Downloads and processed Seurat objects are cached locally. The processed tables and pairwise summaries are included in `results/revision/multidonor/`.

## Composition and network references

`scripts/revision/composition.R` obtains the public TCGA PanCanAtlas leukocyte fractions, CIBERSORT fractions and ABSOLUTE purity estimates from GDC and calculates ESTIMATE scores. The script records the exact source URLs and covariate definitions. It requires the `estimate` package.

`scripts/revision/downstream.R` downloads STRING v12.0 human physical links and TRRUST v2 human TF-target relationships. The downloaded databases remain in local caches; the resulting benchmark tables are included. These are exploratory reference comparisons, not a perturbational validation dataset.

## miRTarBase v10

Used for the validated-target heatmap (Fig S2):

- <https://mirtarbase.cuhk.edu.cn/>

Download the v10 release; one Excel/CSV per organism.
