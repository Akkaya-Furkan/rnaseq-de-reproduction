# RNA-seq differential expression — reproducing a published study

Reproduction of the dexamethasone response of human airway smooth muscle cells
(GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778), Himes et al. 2014)
with pyDESeq2, compared against the authors' own differential expression table.

## Question

Can the differentially expressed genes reported for GSE52778 (dexamethasone vs untreated) be reproduced
with an independent pipeline, and how close are the gene sets and gene rankings?

## Data

| Item | Source |
|---|---|
| Counts | [recount3](https://rna.recount.bio) gene-level sums, Gencode v26, project SRP033351 (= GSE52778), 16 runs |
| Published DE result | GEO supplementary file `GSE52778_Dex_vs_Untreated_gene_exp.diff.gz` (Cuffdiff) |
| Design | 4 airway smooth muscle donors × {untreated, dexamethasone, albuterol, albuterol + dexamethasone} |

Primary analysis: the 8 untreated and dexamethasone-only samples (4 vs 4, paired by donor).
`scripts/download_data.py` fetches everything and writes SHA-256 checksums; raw data are not committed.

## Approach

1. **Counts.** recount3 sums are base-pair coverage; they were divided by each sample's mean mapped read length and
   rounded to pseudo read counts. Ensembl ids were mapped to symbols; where several ids share a symbol, the id with
   the largest total count was kept. Genes with fewer than 10 total counts over the 8 samples were dropped (26,150 kept).
2. **QC.** Library size, PCA of log2(normalised + 1), sample-to-sample Spearman correlation
   (`results/figures/qc_library_pca_correlation.png`).
3. **Differential expression.** pyDESeq2 0.5.4, design `~ cell_line + dex`, Wald test, Benjamini–Hochberg
   (padj < 0.05), apeglm-type LFC shrinkage for plotting; unshrunk log2FC and Wald statistics are used for comparisons.
4. **Comparison with the publication.** The Cuffdiff `log2(fold_change)` column is log2(Untreated/Dex); its sign
   was flipped to Dex/Untreated (checked against the `value_1`/`value_2` columns). Only genes with Cuffdiff status
   `OK` that are also tested by pyDESeq2 form the common universe (11,975 genes). Reported: Jaccard index of the
   significant sets, overlap of the top-N genes ranked by |signed score| (published score = sign(log2FC) × −log10 p,
   p floored at 2.2e-16 where Cuffdiff reports 0), Spearman correlation of scores and of log2FC, direction concordance.
5. **Checks.** Donor-label sign-flip checks, an unpaired model, an |log2FC| > 1 filter, and a 16-sample model.
6. An independent review (Codex) of leakage, statistics and number-code consistency was run; its findings are
   incorporated (permutation wording, figure labels, a chance-expectation edge case, exploratory status of the
   label-corrected 16-sample model).

## Results

All numbers come from `results/metrics.json` (produced by `scripts/run_analysis.py`).

| | Published (Cuffdiff) | pyDESeq2 |
|---|---|---|
| Genes tested | 13,513 (status OK) | 21,080 with padj |
| Significant at 5% FDR | 316 (302 in common universe) | 4,158 (3,275 in common universe) |

* **Recovery.** 301 of the 302 published significant genes in the common universe are also significant in
  pyDESeq2 (99.7%), all 301 in the same direction. The only miss is RASD1 (published q = 0.016).
  14 published significant genes are not tested by pyDESeq2 (symbol mismatch or filtered).
* **Jaccard of the significant sets: 0.092**, because pyDESeq2 calls ~11× more genes: only 9.2% of its calls are in the
  published list. Median baseMean is 1,136 for shared genes versus 760 for pyDESeq2-only genes; the pyDESeq2-only
  genes have a smaller median |log2FC| (0.55 vs 1.67 for shared genes). Restricting both sets to |log2FC| > 1 raises
  the Jaccard to 0.39 (302 published, 760 pyDESeq2).
* **Top-N overlap** (Jaccard of the top-N by |score|; shared genes in brackets; expected by chance < 0.01 for N ≤ 50):
  N=10: 0.54 (7), N=25: 0.47 (16), N=50: 0.39 (28), N=100: 0.57 (73), N=200: 0.57 (145), N=316: 0.61 (239).
  At N=100 the chance expectation is 0.8 shared genes.
* **Rank agreement.** Spearman ρ = 0.963 for signed scores and 0.960 for log2FC (n = 11,975); the log2FC sign agrees for 92.5%
  of all common genes and 100% of the published significant genes.
* **Known dexamethasone genes** (`results/tables/known_dex_genes.csv`): DUSP1, KLF15, PER1, TSC22D3, FKBP5,
  CRISPLD2, ZBTB16, MT2A are strongly up in both (pyDESeq2 log2FC 2.2–7.2, padj ≤ 1.6e-19). ANGPTL4 is
  significant in pyDESeq2 (padj 3.0e-5, log2FC 0.86) but not in Cuffdiff (q = 0.98), an example of the power difference.
* **QC.** 8 libraries of 19.2–31.6 M assigned counts; PC1 (30%) separates dexamethasone from untreated, PC2 (21%)
  largely reflects donor.
* **Sign-flip checks.** Swapping dexamethasone/untreated labels in 2 of 4 donors (6 patterns) gave 0 significant genes in
  all 6. Swapping 1 or 3 donors keeps part of the true effect and gave 0–23 significant genes. These are sensitivity
  checks, not a calibrated null.
* **Unpaired model** (`~ dex`): 3,230 significant, Jaccard with the published set 0.114.

Figures: `results/figures/volcano_dex_vs_untreated.png`, `overlap_published_vs_pydeseq2.png`,
`qc_library_pca_correlation.png`. Tables: `results/tables/` (DE table `de_pydeseq2_dex_vs_untreated.csv`,
merged comparison, top-N overlap, genes unique to either side, QC, sensitivity analyses).

### Why the gene lists differ

In this comparison the Cuffdiff table is much less powerful (it flags 316 genes against 4,158) than
a paired negative-binomial model; whether Cuffdiff modelled donor pairing is not documented. The published list is a high-effect, high-expression subset of what pyDESeq2
finds. The effect sizes and signed scores are nearly the same (ρ ≈ 0.96), so the disagreement is mostly in where the
significance threshold falls rather than in what the data show. Other differences: counts (recount3 STAR/Gencode v26
coverage-derived versus the authors' alignment/annotation), gene symbols as the join key, and pseudo-counts from coverage.

## Limitations

* The original Cuffdiff run (version, annotation, exact sample grouping, pooling of replicates) is not documented in the
  deposited table; differences are reported as method differences. Whether the table pools only dexamethasone-alone
  samples or also the albuterol groups is not stated and was not determined.
* Counts are not the authors' counts: recount3 coverage sums scaled by mean mapped length are an approximation of read counts.
* The join is by gene symbol; 14 of 316 published significant genes were not matched or tested. Duplicate symbols keep one Ensembl id.
* **Label inconsistency in the albuterol arms.** For donors N052611 and N080611, the "Albuterol" sample has a
  dexamethasone-like expression profile and the "Albuterol + Dexamethasone" sample does not
  (`results/tables/dex_response_score_all16.csv`). A 16-sample additive model with the labels as deposited finds 0
  significant genes. Swapping those two pairs' labels, a data-driven and therefore circular correction, gives 8,805
  significant genes; this is exploratory, FDR is not interpretable, and it is not used for any headline result.
  The cause (sample mix-up in GEO/SRA metadata or in processing) was not investigated. The primary 8-sample analysis
  uses only untreated and dexamethasone-only samples, whose scores are consistent with their labels.
* Four donors give a small sample; dispersion estimation relies on 3 residual degrees of freedom. pyDESeq2 differs
  from R DESeq2 in details (no direct R run was made here, so agreement with R DESeq2 was not tested).
* Jaccard is sensitive to the very different set sizes, and the published score is a p-value-based surrogate
  (Cuffdiff p-values floor at 2.2e-16 and some test statistics are sentinels), so top-N ties at the top are broken by a tiny
  term and the exact top-10 order is approximate.
* No salmon quantification from FASTQ was done; the published/recount3 count route was used.

## Reproduce

```bash
python -m venv .venv && . .venv/bin/activate   # or: conda env create -f environment.yml
pip install -r requirements.txt
python scripts/download_data.py
python scripts/run_analysis.py                  # ~10 min on 12 cores, mostly the permutation/sensitivity fits
pytest
```

Tested with Python 3.14.8 on Linux.

## Layout

`src/rnaseq_de/` loading, statistics and plotting helpers · `scripts/` download and analysis ·
`tests/` unit tests (10) · `results/` tables, figures, `metrics.json`.
