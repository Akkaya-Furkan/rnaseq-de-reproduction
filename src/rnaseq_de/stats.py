"""Overlap and rank-agreement metrics."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def top_n(df: pd.DataFrame, score: str, n: int) -> list[str]:
    """Index labels of the n rows with the largest |score| (stable tie order)."""
    return df[score].abs().sort_values(ascending=False, kind="mergesort").index[:n].tolist()


def top_n_overlap(published: pd.DataFrame, mine: pd.DataFrame, score_pub: str, score_mine: str, ns) -> pd.DataFrame:
    """Overlap of the top-N genes by |score| in both tables, restricted to genes present in both."""
    common = published.index.intersection(mine.index)
    if len(common) == 0:
        raise ValueError("no genes in common")
    p, m = published.loc[common], mine.loc[common]
    out = []
    for n in ns:
        a, b = set(top_n(p, score_pub, n)), set(top_n(m, score_mine, n))
        out.append({"n": n, "shared": len(a & b), "jaccard": jaccard(a, b),
                    "expected_shared_by_chance": len(a) * len(b) / len(common)})
    return pd.DataFrame(out)


def rank_correlation(x: pd.Series, y: pd.Series) -> dict:
    common = x.index.intersection(y.index)
    x, y = x.loc[common], y.loc[common]
    ok = x.notna() & y.notna()
    rho, p = stats.spearmanr(x[ok], y[ok])
    return {"n_genes": int(ok.sum()), "spearman_rho": float(rho), "spearman_p": float(p)}


def direction_concordance(x: pd.Series, y: pd.Series) -> float:
    common = x.index.intersection(y.index)
    return float((np.sign(x.loc[common]) == np.sign(y.loc[common])).mean())
