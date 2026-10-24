"""Figures."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def volcano(res: pd.DataFrame, path, label_n: int = 12, alpha: float = 0.05) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    x = res["lfc_shrunk"]
    y = -np.log10(res["padj"].clip(lower=1e-300))
    sig = res["padj"] < alpha
    ax.scatter(x[~sig], y[~sig], s=5, c="#b0b0b0", label="padj ≥ 0.05", rasterized=True)
    ax.scatter(x[sig], y[sig], s=6, c="#c0392b", label=f"padj < 0.05 (n={int(sig.sum())})", rasterized=True)
    for g in res.loc[sig].assign(_a=lambda d: d["stat"].abs()).nlargest(label_n, "_a").index:
        ax.annotate(g, (x[g], y[g]), fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.axhline(-np.log10(alpha), ls=":", c="k", lw=0.6)
    ax.set_xlabel("log2 fold change, dexamethasone vs untreated (shrunk)")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title("pyDESeq2, paired by donor (GSE52778, 4 vs 4)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def overlap_figure(merged: pd.DataFrame, topn: pd.DataFrame, jacc_sig: float, path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axes[0]
    ax.scatter(merged["log2fc_pub"], merged["lfc_mine"], s=4, alpha=0.4, c="#34495e", rasterized=True)
    lim = float(np.nanpercentile(np.abs(merged[["log2fc_pub", "lfc_mine"]].to_numpy()), 99.5))
    ax.plot([-lim, lim], [-lim, lim], "r--", lw=0.8)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("published (Cuffdiff) log2FC"); ax.set_ylabel("pyDESeq2 unshrunk log2FC")
    ax.set_title("Effect sizes")
    ax = axes[1]
    ax.scatter(merged["stat_pub"], merged["stat_mine"], s=4, alpha=0.4, c="#34495e", rasterized=True)
    ax.set_xlabel("published signed score: sign(log2FC) x -log10(p), p floored at 2.2e-16"); ax.set_ylabel("pyDESeq2 Wald statistic")
    ax.set_title("Test statistics")
    ax = axes[2]
    n, e = topn["n"], topn["expected_shared_by_chance"]
    ax.plot(n, topn["jaccard"], "o-", c="#c0392b", label="Jaccard, top-N by |stat|")
    ax.plot(n, e / (2 * n - e), "k:", label="chance (approx.: E[K]/(2N-E[K]))")
    ax.axhline(jacc_sig, c="#2980b9", ls="--", label=f"significant sets: {jacc_sig:.2f}")
    ax.set_xscale("log"); ax.set_ylim(0, 1)
    ax.set_xlabel("N"); ax.set_ylabel("Jaccard index"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("Top-gene overlap")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def qc_figure(lib_sizes: pd.Series, logn: pd.DataFrame, labels: list[str], path) -> None:
    from sklearn.decomposition import PCA

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axes[0]
    ax.bar(range(len(lib_sizes)), lib_sizes.values / 1e6, color=["#c0392b" if "Dex" in l else "#7f8c8d" for l in labels])
    ax.set_xticks(range(len(lib_sizes))); ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("assigned counts (millions)"); ax.set_title("Library size (red = Dex)")
    ax = axes[1]
    pca = PCA(n_components=2).fit(logn.T.values)
    pc = pca.transform(logn.T.values)
    donors = sorted({l.split("_")[0] for l in labels})
    cmap = dict(zip(donors, ["#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]))
    for (x, y), l in zip(pc, labels):
        ax.scatter(x, y, c=cmap[l.split("_")[0]], marker="^" if "Dex" in l else "o", s=60)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.0f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.0f}%)")
    ax.set_title("PCA, log2(normalised+1); colour=donor, triangle=Dex", fontsize=9)
    ax = axes[2]
    corr = logn.corr(method="spearman").values
    im = ax.imshow(corr, vmin=corr.min(), vmax=1, cmap="viridis")
    ax.set_xticks(range(len(labels))); ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6); ax.set_yticklabels(labels, fontsize=6)
    fig.colorbar(im, ax=ax, label="Spearman")
    ax.set_title("Sample correlation")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
