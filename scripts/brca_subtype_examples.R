#!/usr/bin/env Rscript

# Pick exemplar Simpson-reversal pairs for the BRCA subtype paper figure.
# Reads the cached results from brca_subtype_analysis.R, ranks reversals by
# (a) effect size of |r_global|, and (b) magnitude of within-context spread,
# and writes the top candidates plus the per-sample data needed to plot them.

suppressPackageStartupMessages({
  source("00-paths.R")
  source(file.path("R", "contextcorr.R"))
})

f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
matched_file <- if (file.exists(f1)) f1 else f2
load(matched_file)

clin <- read.delim(file.path(paths$clean, "BRCA_clinicalMatrix.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)

brca_idx <- which(cancerType == "BRCA" & sampleType == "tumor")
brca_ids15 <- gsub("\\.", "-", substr(colnames(mRNAdata)[brca_idx], 1, 15))
m <- match(brca_ids15, clin$sampleID)

pam50_raw <- clin$PAM50Call_RNAseq[m]
pam50 <- ifelse(pam50_raw %in% c("Basal", "Her2", "LumA", "LumB", "Normal"),
                pam50_raw, NA_character_)
er <- ifelse(clin$ER_Status_nature2012[m] %in% c("Positive","Negative"),
             clin$ER_Status_nature2012[m], NA_character_)
pr <- ifelse(clin$PR_Status_nature2012[m] %in% c("Positive","Negative"),
             clin$PR_Status_nature2012[m], NA_character_)
h2 <- ifelse(clin$HER2_Final_Status_nature2012[m] %in% c("Positive","Negative"),
             clin$HER2_Final_Status_nature2012[m], NA_character_)
ihc <- rep(NA_character_, length(er))
ok_ihc <- !is.na(er) & !is.na(pr) & !is.na(h2)
ihc[ok_ihc & h2 == "Positive"] <- "HER2+"
ihc[ok_ihc & h2 == "Negative" & (er == "Positive" | pr == "Positive")] <- "HR+"
ihc[ok_ihc & h2 == "Negative" & er == "Negative" & pr == "Negative"] <- "TNBC"

pick_examples <- function(label_name, label_vec, top_n = 30) {
  cache <- readRDS(file.path("results", "cache",
                             sprintf("brca_subtype_%s.rds", label_name)))
  rev_mask <- cache$reversal_mask
  rg <- cache$r_global
  ca <- cache$cor_array  # G x M x L
  spread <- apply(ca, 1:2, function(x) {
    if (sum(is.finite(x)) < 2) return(NA_real_)
    diff(range(x, na.rm = TRUE))
  })
  # Score: |r_global| * spread, only on reversals
  score <- abs(rg) * spread
  score[!rev_mask] <- NA_real_

  ord <- order(score, decreasing = TRUE, na.last = NA)
  G <- nrow(rg); M <- ncol(rg)
  if (length(ord) == 0) return(NULL)
  ord <- ord[seq_len(min(top_n, length(ord)))]
  cc <- ((ord - 1) %/% G) + 1L
  rr <- ((ord - 1) %% G) + 1L

  top <- data.frame(
    context = label_name,
    rank = seq_along(ord),
    gene = rownames(rg)[rr],
    mir = colnames(rg)[cc],
    r_global = rg[cbind(rr, cc)],
    spread = spread[cbind(rr, cc)],
    score = score[cbind(rr, cc)],
    stringsAsFactors = FALSE
  )
  # Per-context r values
  L <- dim(ca)[3]
  for (l in seq_len(L)) {
    nm <- paste0("r_", dimnames(ca)[[3]][l])
    top[[nm]] <- ca[cbind(rr, cc, l)]
  }
  top
}

ex_pam50 <- pick_examples("PAM50", pam50)
ex_ihc <- pick_examples("IHC", ihc)

write.csv(ex_pam50, file.path("results", "brca_subtype_examples_PAM50.csv"),
          row.names = FALSE)
write.csv(ex_ihc, file.path("results", "brca_subtype_examples_IHC.csv"),
          row.names = FALSE)

# Build per-sample data for the very top reversal under each context, so the
# Rmd can produce scatter panels without re-doing the heavy work.
make_pair_payload <- function(label_name, label_vec, gene, mir) {
  keep <- which(!is.na(label_vec))
  data.frame(
    sample = colnames(mRNAdata)[brca_idx][keep],
    context = label_vec[keep],
    mRNA = mRNAdata[gene, brca_idx][keep],
    miR = miRdata[mir, brca_idx][keep],
    stringsAsFactors = FALSE
  )
}

if (!is.null(ex_pam50) && nrow(ex_pam50) >= 1) {
  pl <- make_pair_payload("PAM50", pam50, ex_pam50$gene[1], ex_pam50$mir[1])
  saveRDS(list(gene = ex_pam50$gene[1], mir = ex_pam50$mir[1], data = pl,
               summary = ex_pam50[1, ]),
          file.path("results", "cache", "brca_top_reversal_PAM50.rds"))
}
if (!is.null(ex_ihc) && nrow(ex_ihc) >= 1) {
  pl <- make_pair_payload("IHC", ihc, ex_ihc$gene[1], ex_ihc$mir[1])
  saveRDS(list(gene = ex_ihc$gene[1], mir = ex_ihc$mir[1], data = pl,
               summary = ex_ihc[1, ]),
          file.path("results", "cache", "brca_top_reversal_IHC.rds"))
}

message("Wrote brca_subtype_examples_*.csv and top reversal cache.")
