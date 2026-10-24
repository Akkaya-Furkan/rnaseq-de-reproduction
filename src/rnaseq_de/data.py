"""Loading and harmonising the count matrix, sample table, gene symbols and the published DE table."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


def load_gene_symbols(gtf_gz: Path) -> pd.Series:
    """Map Ensembl gene id (with version) -> gene symbol from the recount3 GTF."""
    ids, names = [], []
    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            m = re.search(r'gene_name "([^"]+)"', f[8])
            ids.append(gid)
            names.append(m.group(1) if m else gid)
    return pd.Series(names, index=ids, name="symbol")


def coverage_to_counts(coverage: pd.DataFrame, mean_mapped_length: pd.Series) -> pd.DataFrame:
    """recount3 gene sums are base-pair coverage; divide by the mean mapped read length
    per sample and round to obtain pseudo read counts (genes x samples)."""
    scaled = coverage.div(mean_mapped_length.reindex(coverage.columns), axis=1)
    return scaled.round().astype(int)


def collapse_duplicate_symbols(counts: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    """Index counts by gene symbol; when several Ensembl ids share a symbol keep the
    one with the largest total count (no summing, so no artificial inflation)."""
    df = counts.copy()
    df["_symbol"] = symbols.reindex(counts.index).values
    df["_total"] = counts.sum(axis=1).values
    df = df.dropna(subset=["_symbol"]).sort_values("_total", ascending=False, kind="mergesort")
    df = df.drop_duplicates("_symbol").set_index("_symbol").drop(columns="_total")
    df.index.name = "symbol"
    return df


def load_sample_table(sra_meta_gz: Path) -> pd.DataFrame:
    """Sample table indexed by SRR run id with cell_line, treatment and GSM accession."""
    meta = pd.read_csv(sra_meta_gz, sep="\t")
    rows = []
    for _, r in meta.iterrows():
        attrs = dict(p.split(";;", 1) for p in r["sample_attributes"].split("|"))
        rows.append({"run": r["external_id"], "gsm": r["sample_name"],
                     "cell_line": attrs["cell line"], "treatment": attrs["treatment"]})
    return pd.DataFrame(rows).set_index("run")


def load_published(diff_gz: Path) -> pd.DataFrame:
    """Authors' Cuffdiff table (GEO supplementary)."""
    d = pd.read_csv(diff_gz, sep="\t").rename(columns={"log2(fold_change)": "log2fc"})
    return d[["gene", "locus", "status", "value_1", "value_2", "log2fc", "test_stat", "p_value", "q_value", "significant"]]
