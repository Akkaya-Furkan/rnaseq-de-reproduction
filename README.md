# RNA-seq differential expression — reproducing a study

**Question.** Can the top differentially expressed genes reported by a public GEO study be reproduced?

## Data
A public GEO RNA-seq dataset (≥ 3 samples per group) with a published DE list.

## Approach
1. Quantification with salmon (or the published count matrix)
2. QC: library size, PCA, sample correlation
3. DESeq2 / pyDESeq2 with LFC shrinkage, padj < 0.05
4. Overlap with the published top genes (Jaccard, rank correlation)

## Deliverables
DE table, overlap report, volcano and overlap figures, write-up.

## Limitations
- The original pipeline may not be fully specified; differences are reported as method differences.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.
