# Mitochondrial gene expression in Parkinson's and Alzheimer's disease brain: an exploratory comparison

An exploratory undergraduate analysis of public microarray data. It has not been peer reviewed, and the results are hypothesis-generating, not conclusions.

## Question

Do Parkinson's disease (PD) and Alzheimer's disease (AD) brain tissue show similar changes in mitochondrial gene expression, and are the changes concentrated in particular mitochondrial processes?

## Data

| Disease | GEO accession | Tissue | Samples |
|---|---|---|---|
| PD | GSE7621 | Substantia nigra | 9 control, 16 PD |
| AD | GSE36980 | Hippocampus (subset of a three-region dataset) | 10 control, 8 AD |

Mitochondrial genes: MitoCarta3.0 (1,136 human genes), downloaded from the Broad Institute.

## Methods

1. PD expression values were on a linear intensity scale and were log2-transformed. AD values were already on a log scale.
2. Differential expression (disease vs control) with `limma`, run separately per dataset.
3. Probes were mapped to gene symbols and restricted to MitoCarta3.0 genes. Gene-level fold-changes average the probes for the same gene.
4. GO biological process enrichment (`clusterProfiler`, Benjamini-Hochberg) on mitochondrial genes with nominal p < 0.05, using all mitochondrial genes tested in that dataset as the background.
5. Cross-disease comparison: hypergeometric overlap test, Spearman correlation of fold-changes, a permutation control against random non-mitochondrial gene sets (2,000 draws, seed 1), and Fisher's exact tests on the direction of change in respiratory chain genes.

## Results

- **Individual genes.** 986 mitochondrial genes were tested in both datasets. No mitochondrial gene in PD passes a genome-wide FDR of 0.05 (best adjusted p = 0.13). The AD hits are modest as well.
- **Enrichment.** In AD, respiratory chain and oxidative phosphorylation terms are enriched among changed mitochondrial genes (adjusted p roughly 6e-07 to 2e-06 for the top terms, fold enrichment about 1.6 to 2.0), and 47 of the 49 respiratory chain genes tested in both datasets are down. In PD, no GO term is enriched.
- **Overlap.** 66 genes are nominally significant in both diseases, against 74.4 expected by chance (hypergeometric p = 0.93). There is no excess of shared individual genes.
- **Concordance.** Across all 986 genes, PD and AD fold-changes are positively correlated (Spearman rho = 0.30). Random sets of non-mitochondrial genes give rho of about 0.19 on average (maximum 0.27 in 2,000 draws; permutation p < 0.0005).
- **Respiratory chain direction.** 38 of the 49 respiratory chain genes are down in PD, a larger share than among the other mitochondrial genes (odds ratio 2.6, 95% CI 1.3 to 5.7, p = 0.005). The gene set was defined from the AD enrichment, so the PD test is an out-of-sample check. The equivalent AD test is circular and is reported only for completeness.

Figures are in `figures/`, tables in `results/`, and summary statistics in `results/summary_stats.txt`.

## Limitations

- Small sample sizes, with different array platforms and different brain regions in the two datasets.
- Postmortem tissue: neuron loss and changes in cell composition can shift mitochondrial gene expression in both diseases without disease-specific mitochondrial pathology. The non-mitochondrial correlation (about 0.19) suggests this contributes.
- Gene lists use a nominal p < 0.05 cutoff because almost nothing passes genome-wide FDR. This is exploratory.
- Mitochondrial genes are co-regulated, so the tests above treat correlated genes as independent and are likely optimistic.
- Probes annotated with several gene symbols were not split.
- The results are correlational and make no causal claim.

## Possible next steps

- Replicate the concordance in independent PD and AD expression datasets.
- Protein-protein interaction network analysis (STRING) and hub genes.
- Cell-type deconvolution to separate cell composition effects from per-cell expression changes.

## Reproducing the analysis

1. Download `Human.MitoCarta3.0.xls` from the Broad Institute MitoCarta3.0 page and set the path at the top of the script.
2. Open `mitochondrial_PD_AD_analysis.R` in RStudio and set the working directory to the script's folder.
3. Run the script. It downloads both GEO datasets and writes `results/` and `figures/`.

R 4.6.1 with the packages loaded at the top of the script. See `results/session_info.txt` after a run for exact versions.
