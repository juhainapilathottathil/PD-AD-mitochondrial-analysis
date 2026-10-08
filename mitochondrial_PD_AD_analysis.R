# =============================================================================
# Exploratory comparison of mitochondrial gene expression changes in
#   Parkinson's disease (GSE7621, substantia nigra) and
#   Alzheimer's disease (GSE36980, hippocampus)
#
# How to run: put this file in the repository folder, open it in RStudio,
# then Session > Set Working Directory > To Source File Location.
# Needs an internet connection (GEO download) and the MitoCarta3.0 file from
# the Broad Institute (see README).
# =============================================================================

library(GEOquery)
library(Biobase)
library(limma)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Looks for the MitoCarta3.0 file in this folder first, then in Downloads
mitocarta_file <- c("Human.MitoCarta3.0.xls", "~/Downloads/Human.MitoCarta3.0.xls")
mitocarta_file <- mitocarta_file[file.exists(mitocarta_file)][1]
if (is.na(mitocarta_file)) {
  stop("Download Human.MitoCarta3.0.xls from the Broad Institute (see README) ",
       "and save it in this folder or in Downloads.")
}

# ---- 1. Mitochondrial gene list (MitoCarta3.0) ------------------------------
mito_genes   <- read_excel(mitocarta_file, sheet = "A Human MitoCarta3.0")
mito_symbols <- mito_genes$Symbol

# ---- 2. Parkinson's disease: GSE7621 ----------------------------------------
gse_pd  <- getGEO("GSE7621", GSEMatrix = TRUE)
pd_data <- gse_pd[[1]]

pd_titles <- pData(pd_data)$title
stopifnot(all(grepl("normal", pd_titles[1:9])), all(grepl("PD", pd_titles[10:25])))
group_pd <- factor(c(rep("control", 9), rep("PD", 16)), levels = c("control", "PD"))

# PD values are on a linear intensity scale (max ~34,000); limma expects log2
ex <- exprs(pd_data)
if (max(ex, na.rm = TRUE) > 100) exprs(pd_data) <- log2(pmax(ex, 1))

fit_pd     <- eBayes(lmFit(exprs(pd_data), model.matrix(~group_pd)))
pd_results <- topTable(fit_pd, coef = 2, number = Inf)
pd_results$ID <- rownames(pd_results)

# Probes with several symbols (e.g. "A /// B") are not split and will not match
probe_to_gene_pd <- fData(pd_data)[, c("ID", "Gene Symbol")]
colnames(probe_to_gene_pd) <- c("ID", "GeneSymbol")
pd_annot        <- merge(pd_results, probe_to_gene_pd, by = "ID")
pd_mito_results <- pd_annot[pd_annot$GeneSymbol %in% mito_symbols, ]
pd_mito_results <- pd_mito_results[order(pd_mito_results$P.Value), ]

# ---- 3. Alzheimer's disease: GSE36980 (hippocampus only) --------------------
gse_ad  <- getGEO("GSE36980", GSEMatrix = TRUE)
ad_data <- gse_ad[[1]]

ad_hi    <- ad_data[, grepl("HI", pData(ad_data)$title)]
group_ad <- factor(ifelse(grepl("non-AD", pData(ad_hi)$title), "control", "AD"),
                   levels = c("control", "AD"))
stopifnot(table(group_ad)[["control"]] == 10, table(group_ad)[["AD"]] == 8)

# Gene symbol = 2nd "//"-separated field of the first "///"-separated entry
extract_gene_symbol <- function(x) {
  if (is.na(x) || x == "---") return(NA_character_)
  first_entry <- strsplit(x, " /// ")[[1]][1]
  fields <- strsplit(first_entry, " // ")[[1]]
  if (length(fields) >= 2) return(trimws(fields[2]))
  NA_character_
}
ad_symbols       <- sapply(fData(ad_hi)$gene_assignment, extract_gene_symbol)
probe_to_gene_ad <- data.frame(ID = fData(ad_hi)$ID, GeneSymbol = unname(ad_symbols))

fit_ad     <- eBayes(lmFit(exprs(ad_hi), model.matrix(~group_ad)))
ad_results <- topTable(fit_ad, coef = 2, number = Inf)
ad_results$ID <- rownames(ad_results)

ad_annot        <- merge(ad_results, probe_to_gene_ad, by = "ID")
ad_mito_results <- ad_annot[ad_annot$GeneSymbol %in% mito_symbols, ]
ad_mito_results <- ad_mito_results[order(ad_mito_results$P.Value), ]

# ---- 4. Gene lists and GO enrichment ----------------------------------------
# Background ("universe") = all mitochondrial genes tested in that dataset.
# Lists use nominal p < 0.05 (exploratory; see README for limitations).
pd_universe <- unique(pd_mito_results$GeneSymbol)
ad_universe <- unique(ad_mito_results$GeneSymbol)
pd_sig <- unique(pd_mito_results$GeneSymbol[pd_mito_results$P.Value < 0.05])
ad_sig <- unique(ad_mito_results$GeneSymbol[ad_mito_results$P.Value < 0.05])

run_enrich <- function(genes, universe) {
  enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
           keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
}
enrich_pd <- run_enrich(pd_sig, pd_universe)
enrich_ad <- run_enrich(ad_sig, ad_universe)

# ---- 5. Cross-disease comparison ---------------------------------------------
common_universe <- intersect(pd_universe, ad_universe)
pd_s    <- intersect(pd_sig, common_universe)
ad_s    <- intersect(ad_sig, common_universe)
overlap <- intersect(pd_s, ad_s)
N       <- length(common_universe)
expected_overlap <- length(pd_s) * length(ad_s) / N
p_overlap <- phyper(length(overlap) - 1, length(pd_s), N - length(pd_s),
                    length(ad_s), lower.tail = FALSE)

# Gene-level fold-changes (probes for the same gene are averaged)
pd_fc  <- tapply(pd_mito_results$logFC, pd_mito_results$GeneSymbol, mean)
ad_fc  <- tapply(ad_mito_results$logFC, ad_mito_results$GeneSymbol, mean)
common <- intersect(names(pd_fc), names(ad_fc))
rho_mito <- cor(pd_fc[common], ad_fc[common], method = "spearman")

# Control: same correlation for random sets of non-mitochondrial genes
pd_all   <- tapply(pd_annot$logFC, pd_annot$GeneSymbol, mean)
ad_all   <- tapply(ad_annot$logFC, ad_annot$GeneSymbol, mean)
non_mito <- setdiff(intersect(names(pd_all), names(ad_all)), mito_symbols)
set.seed(1)
null_rho <- replicate(2000, {
  g <- sample(non_mito, length(common))
  cor(pd_all[g], ad_all[g], method = "spearman")
})
p_perm <- (sum(null_rho >= rho_mito) + 1) / (length(null_rho) + 1)

# Respiratory chain genes: every human gene annotated to GO:0022904 (respiratory
# electron transport chain, including child terms) in the GO annotation database.
# The set is defined independently of both datasets, so testing it in AD and in PD
# is not circular. Are these genes down more often than other mitochondrial genes?
resp_go <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0022904",
                                 keytype = "GOALL", columns = "SYMBOL")
resp <- intersect(unique(resp_go$SYMBOL), common)
rest <- setdiff(common, resp)
dir_table <- function(fc) {
  matrix(c(sum(fc[resp] < 0), sum(fc[resp] > 0), sum(fc[rest] < 0), sum(fc[rest] > 0)),
         nrow = 2, byrow = TRUE,
         dimnames = list(c("respiratory chain", "other mitochondrial"), c("down", "up")))
}
dir_pd <- dir_table(pd_fc);  fisher_pd <- fisher.test(dir_pd)
dir_ad <- dir_table(ad_fc);  fisher_ad <- fisher.test(dir_ad)
# Also compare the fold-changes themselves (not just their sign)
wilcox_pd <- wilcox.test(pd_fc[resp], pd_fc[rest])
wilcox_ad <- wilcox.test(ad_fc[resp], ad_fc[rest])

# ---- 6. Save results ----------------------------------------------------------
sink("results/summary_stats.txt")
cat("Mitochondrial genes tested: PD", length(pd_universe), "| AD", length(ad_universe),
    "| both", N, "\n")
cat("Nominal p < 0.05 (within shared universe): PD", length(pd_s), "| AD", length(ad_s), "\n")
cat("Best genome-wide adj.P.Val among mito genes: PD", signif(min(pd_mito_results$adj.P.Val), 3),
    "| AD", signif(min(ad_mito_results$adj.P.Val), 3), "\n")
cat("PD GO terms enriched:", nrow(as.data.frame(enrich_pd)),
    "| AD GO terms enriched:", nrow(as.data.frame(enrich_ad)), "\n\n")
cat("Overlap of nominally significant genes:", length(overlap), "observed vs",
    round(expected_overlap, 1), "expected; hypergeometric p =", signif(p_overlap, 3), "\n\n")
cat("Spearman rho of PD vs AD fold-changes (mitochondrial genes, n =", length(common), "):",
    round(rho_mito, 3), "\n")
cat("Random non-mitochondrial gene sets: mean", round(mean(null_rho), 3),
    "| max", round(max(null_rho), 3), "| permutation p =", signif(p_perm, 3), "\n\n")
cat("Respiratory chain genes (GO:0022904 annotation) tested in both datasets:", length(resp), "\n")
cat("\nPD direction:\n"); print(dir_pd); print(fisher_pd)
cat("Median logFC: respiratory chain", round(median(pd_fc[resp]), 3),
    "| other mitochondrial", round(median(pd_fc[rest]), 3), "\n")
print(wilcox_pd)
cat("\nAD direction:\n"); print(dir_ad); print(fisher_ad)
cat("Median logFC: respiratory chain", round(median(ad_fc[resp]), 3),
    "| other mitochondrial", round(median(ad_fc[rest]), 3), "\n")
print(wilcox_ad)
sink()

write.csv(pd_mito_results, "results/PD_mitochondrial_DE_results.csv", row.names = FALSE)
write.csv(ad_mito_results, "results/AD_mitochondrial_DE_results.csv", row.names = FALSE)
write.csv(as.data.frame(enrich_ad), "results/AD_GO_enrichment.csv", row.names = FALSE)
write.csv(data.frame(GeneSymbol = common,
                     PD_logFC = as.numeric(pd_fc[common]),
                     AD_logFC = as.numeric(ad_fc[common]),
                     respiratory_chain = common %in% resp),
          "results/PD_AD_gene_level_fold_changes.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

# ---- 7. Figures ----------------------------------------------------------------
p1 <- dotplot(enrich_ad, showCategory = 10) +
  ggtitle("AD hippocampus: GO enrichment among mitochondrial genes")
ggsave("figures/AD_GO_enrichment_dotplot.png", p1, width = 8, height = 6)

sc <- data.frame(PD = as.numeric(pd_fc[common]), AD = as.numeric(ad_fc[common]),
                 respiratory_chain = common %in% resp)
p2 <- ggplot(sc, aes(PD, AD, colour = respiratory_chain)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(alpha = 0.6) +
  labs(x = "PD log2 fold-change (substantia nigra)",
       y = "AD log2 fold-change (hippocampus)",
       title = "Mitochondrial genes: fold-changes in PD vs AD") +
  theme_minimal()
ggsave("figures/PD_AD_fold_change_scatter.png", p2, width = 7, height = 6)

p3 <- ggplot(data.frame(rho = null_rho), aes(rho)) +
  geom_histogram(bins = 40, fill = "grey70") +
  geom_vline(xintercept = rho_mito, colour = "red") +
  labs(x = "Spearman rho, PD vs AD fold-changes",
       y = "Random non-mitochondrial gene sets",
       title = "Observed mitochondrial concordance (red) vs random gene sets") +
  theme_minimal()
ggsave("figures/permutation_null_vs_observed.png", p3, width = 7, height = 5)
