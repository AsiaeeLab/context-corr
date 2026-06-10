# Figure 2C: CITE-seq protein-RNA global vs residualized correlation (joint context).
# Each point is a protein-RNA pair; orange = sign flips between global and
# mean-residualized correlation (37.9%). Rendered at a compact physical size with
# enlarged axis-title fonts so it visually matches panels A/B when height-locked.
suppressMessages(library(ggplot2))

cache <- "results/cache/sc_200_citeseq_protein_rna_contexts/pair_stats_Joint_RNA_ADT.rds"
out   <- "figures/Figure_CITEseq_GlobalVsResid_panel.pdf"

sp <- readRDS(cache)$summary_pairs
ok <- is.finite(sp$r_global) & is.finite(sp$r_resid)
df <- data.frame(rg = sp$r_global[ok], rr = sp$r_resid[ok])
flip_rate <- mean(sign(df$rr) != sign(df$rg))
df$flip <- factor(ifelse(sign(df$rr) != sign(df$rg),
                         sprintf("sign flip: %.1f%%", 100 * flip_rate), "no flip"),
                  levels = c(sprintf("sign flip: %.1f%%", 100 * flip_rate), "no flip"))
df <- df[order(df$flip == levels(df$flip)[1]), ]  # draw flips on top
message(sprintf("CITE-seq joint global-vs-resid: n=%d, sign-flip=%.1f%%", nrow(df), 100 * flip_rate))

p <- ggplot(df, aes(rg, rr, color = flip)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_point(size = 0.45, alpha = 0.5) +
  scale_color_manual(values = setNames(c("#E8601C", "grey70"), levels(df$flip))) +
  labs(x = expression(r[global]), y = expression(r[residualized])) +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.title = element_blank(), legend.key = element_blank(),
        legend.background = element_blank(), legend.text = element_text(size = 7.5),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0, "pt"),
        axis.title = element_text(size = 13), axis.text = element_text(size = 9),
        panel.grid.minor = element_blank())

ggsave(out, p, width = 2.7, height = 2.7)
message("Wrote ", out)
