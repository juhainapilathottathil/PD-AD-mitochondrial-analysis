# Cross-disease mitochondrial transcriptional concordance in Parkinson's and Alzheimer's disease

An exploratory undergraduate analysis of public brain microarray datasets investigating whether mitochondrial gene-expression changes show similar patterns in Parkinson's disease (PD) and Alzheimer's disease (AD).

This project is hypothesis-generating and has not been peer reviewed. The analyses are intended to explore patterns in public transcriptomic data and do not establish causal mechanisms.

## Research question

Do Parkinson's disease and Alzheimer's disease brain tissue show similar mitochondrial gene-expression changes, and are these changes concentrated in particular mitochondrial processes?

## Biological motivation

Mitochondria are essential for ATP production, redox regulation, metabolism, calcium handling, and cellular stress responses. Neurons have particularly high energetic demands, making mitochondrial biology relevant to neurodegenerative disease.

Parkinson's disease and Alzheimer's disease have distinct pathological features, but both have been associated with mitochondrial dysfunction, altered energy metabolism, oxidative stress, and changes in cellular homeostasis.

This project asks whether mitochondrial transcriptional changes show measurable cross-disease concordance in public brain expression datasets.

## Data

| Disease | GEO accession | Brain region | Samples |
|---|---|---|---|
| Parkinson's disease | GSE7621 | Substantia nigra | 9 control, 16 PD |
| Alzheimer's disease | GSE36980 | Hippocampus subset | 10 control, 8 AD |

Mitochondrial genes were defined using **MitoCarta3.0**, a curated catalogue of human mitochondrial genes downloaded from the Broad Institute.

After platform annotation and filtering, the analysis tested:

- 996 mitochondrial genes in the PD dataset
- 1,019 mitochondrial genes in the AD dataset
- 986 mitochondrial genes represented in both datasets

## Analysis workflow

```text
Public GEO microarray data
          ↓
Expression preprocessing
          ↓
Differential expression with limma
          ↓
Probe → gene annotation
          ↓
MitoCarta3.0 filtering
          ↓
Gene-level fold changes
          ↓
┌───────────────────────────────────────┐
│ Functional enrichment                 │
│ Cross-disease gene overlap            │
│ Fold-change concordance                │
│ Permutation control                    │
│ Respiratory-chain direction analysis  │
└───────────────────────────────────────┘