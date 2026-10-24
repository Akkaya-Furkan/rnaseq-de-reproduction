#!/usr/bin/env python3
"""Reproduce the GSE52778 dexamethasone-vs-untreated DE result with pyDESeq2 and compare it
with the authors' published Cuffdiff table. Run scripts/download_data.py first."""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rnaseq_de import plots, stats  # noqa: E402
from rnaseq_de.data import (collapse_duplicate_symbols, coverage_to_counts, load_gene_symbols,  # noqa: E402
                            load_published, load_sample_table)

DATA, TAB, FIG = ROOT / "data", ROOT / "results" / "tables", ROOT / "results" / "figures"
ALPHA, MIN_TOTAL_COUNT, SEED = 0.05, 10, 0
P_FLOOR = 2.22045e-16  # Cuffdiff reports p = 0 below this value
inference = DefaultInference(n_cpus=4)


def run_deseq(counts: pd.DataFrame, meta: pd.DataFrame, design: str, factor: str, shrink: bool = True):
    """counts: genes x samples. Returns result table (and the fitted dds)."""
    dds = DeseqDataSet(counts=counts.T, metadata=meta, design=design, refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=[factor, "yes", "no"], alpha=ALPHA, inference=inference, quiet=True)
    ds.summary()
    res = ds.results_df.rename(columns={"log2FoldChange": "lfc", "lfcSE": "lfc_se", "pvalue": "pvalue"})
    if shrink:
        coeff = f"{factor}[T.yes]"
        ds.lfc_shrink(coeff=coeff)
        res["lfc_shrunk"] = ds.results_df["log2FoldChange"]
    return res, dds


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    np.random.seed(SEED)

    # ---- inputs --------------------------------------------------------------------
    cov = pd.read_csv(DATA / "gene_sums.tsv.gz", sep="\t", comment="#", index_col=0)
    qc = pd.read_csv(DATA / "recount_qc.tsv.gz", sep="\t").set_index("external_id")
    samples = load_sample_table(DATA / "sra_meta.tsv.gz")
    symbols = load_gene_symbols(DATA / "gene_annotation_G026.gtf.gz")
    pub = load_published(DATA / "GSE52778_Dex_vs_Untreated_gene_exp.diff.gz")

    counts_all = collapse_duplicate_symbols(
        coverage_to_counts(cov, qc["star.average_mapped_length"]), symbols)
    samples["dex"] = samples["treatment"].isin(["Dexamethasone", "Albuterol_Dexamethasone"]).map({True: "yes", False: "no"})
    samples["alb"] = samples["treatment"].isin(["Albuterol", "Albuterol_Dexamethasone"]).map({True: "yes", False: "no"})
    samples["label"] = samples["cell_line"] + "_" + samples["treatment"].replace(
        {"Untreated": "Ctrl", "Dexamethasone": "Dex", "Albuterol": "Alb", "Albuterol_Dexamethasone": "AlbDex"})
    samples = samples.loc[counts_all.columns]
    assert counts_all.shape[1] == 16 and samples["cell_line"].nunique() == 4

    # ---- primary analysis: untreated vs dexamethasone only, paired by donor ---------
    prim = samples[samples["treatment"].isin(["Untreated", "Dexamethasone"])]
    counts = counts_all[prim.index]
    counts = counts[counts.sum(axis=1) >= MIN_TOTAL_COUNT]
    meta = prim[["cell_line", "dex"]].copy()
    meta["dex"] = pd.Categorical(meta["dex"], categories=["no", "yes"])
    n_genes_tested = counts.shape[0]

    res, dds = run_deseq(counts, meta, "~cell_line + dex", "dex")
    res.index.name = "gene"
    res.sort_values("padj").to_csv(TAB / "de_pydeseq2_dex_vs_untreated.csv")

    # ---- QC ------------------------------------------------------------------------
    norm = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names).T
    logn = np.log2(norm + 1)
    lib = counts_all[prim.index].sum()
    qc_tab = pd.DataFrame({"label": prim["label"], "assigned_counts": lib,
                           "size_factor": dds.obs["size_factors"]})
    qc_tab.to_csv(TAB / "qc_samples.csv")
    plots.qc_figure(lib, logn, prim["label"].tolist(), FIG / "qc_library_pca_correlation.png")
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2).fit(logn.T.values)
    pc_ratio = pca.explained_variance_ratio_.tolist()

    # ---- comparison with the published Cuffdiff table ------------------------------
    p = pub[pub["status"] == "OK"].copy()
    # Cuffdiff's log2(fold_change) column is log2(value_2/value_1) = Untreated/Dex; flip to Dex/Untreated.
    p["lfc_dex"] = -p["log2fc"]
    ts = p["test_stat"].where(p["test_stat"].abs() < 1e100, np.nan).abs().fillna(0).clip(upper=100)
    p["score"] = np.sign(p["lfc_dex"]) * (-np.log10(p["p_value"].clip(lower=P_FLOOR)) + 1e-6 * ts)
    p["sig"] = p["significant"].eq("yes")
    p = p.drop_duplicates("gene").set_index("gene")

    m = res.dropna(subset=["padj"]).copy()
    common = p.index.intersection(m.index)
    merged = pd.DataFrame({
        "log2fc_pub": p.loc[common, "lfc_dex"], "stat_pub": p.loc[common, "score"], "q_pub": p.loc[common, "q_value"],
        "sig_pub": p.loc[common, "sig"], "lfc_mine": m.loc[common, "lfc"], "lfc_shrunk_mine": m.loc[common, "lfc_shrunk"],
        "stat_mine": m.loc[common, "stat"], "padj_mine": m.loc[common, "padj"]})
    merged["sig_mine"] = merged["padj_mine"] < ALPHA
    merged.to_csv(TAB / "comparison_merged.csv")

    ns = [10, 25, 50, 100, 200, 316]
    topn = stats.top_n_overlap(merged.rename(columns={"stat_pub": "s"}), merged.rename(columns={"stat_mine": "s"}),
                               "s", "s", ns)
    # significant sets within the common universe
    A, B = set(merged.index[merged.sig_pub]), set(merged.index[merged.sig_mine])
    # top-N-by-|stat| comparison above uses the same universe `common`
    sig_pub_all = set(p.index[p.sig])
    sig_mine_all = set(m.index[m["padj"] < ALPHA])
    sig_both_lfc = {g for g in A & B if np.sign(merged.loc[g, "log2fc_pub"]) == np.sign(merged.loc[g, "lfc_mine"])}
    topn.to_csv(TAB / "overlap_topN.csv", index=False)

    metrics = {
        "n_samples": int(counts.shape[1]),
        "genes_tested_pydeseq2_after_filter": int(n_genes_tested),
        "genes_with_padj_pydeseq2": int(res["padj"].notna().sum()),
        "pydeseq2_sig_padj_lt_0.05": int((res["padj"] < ALPHA).sum()),
        "pydeseq2_sig_padj_lt_0.05_abs_lfc_gt_1": int(((res["padj"] < ALPHA) & (res["lfc"].abs() > 1)).sum()),
        "pydeseq2_up": int(((res["padj"] < ALPHA) & (res["lfc"] > 0)).sum()),
        "pydeseq2_down": int(((res["padj"] < ALPHA) & (res["lfc"] < 0)).sum()),
        "published_genes_status_OK": int(len(p)),
        "published_sig_q_lt_0.05_total": int(len(sig_pub_all)),
        "common_tested_genes": int(len(common)),
        "published_sig_in_common": int(len(A)),
        "pydeseq2_sig_in_common": int(len(B)),
        "sig_intersection": int(len(A & B)),
        "sig_jaccard": stats.jaccard(A, B),
        "published_sig_recovered_fraction": len(A & B) / len(A),
        "pydeseq2_sig_also_published_fraction": len(A & B) / len(B),
        "sig_intersection_same_direction": int(len(sig_both_lfc)),
        "published_sig_not_in_common_universe": int(len(sig_pub_all - set(common))),
        "spearman_signed_score": stats.rank_correlation(merged["stat_pub"], merged["stat_mine"]),
        "spearman_log2fc": stats.rank_correlation(merged["log2fc_pub"], merged["lfc_mine"]),
        "direction_concordance_all_common": stats.direction_concordance(merged["log2fc_pub"], merged["lfc_mine"]),
        "direction_concordance_published_sig": stats.direction_concordance(
            merged.loc[list(A), "log2fc_pub"], merged.loc[list(A), "lfc_mine"]),
        "pca_explained_variance_pc1_pc2": pc_ratio,
        "library_size_min_max": [int(lib.min()), int(lib.max())],
    }

    plots.volcano(res.dropna(subset=["padj"]), FIG / "volcano_dex_vs_untreated.png")
    plots.overlap_figure(merged, topn, metrics["sig_jaccard"], FIG / "overlap_published_vs_pydeseq2.png")

    # positive-control genes known from the dexamethasone literature
    known = ["DUSP1", "KLF15", "PER1", "TSC22D3", "FKBP5", "CRISPLD2", "PDK4", "ZBTB16", "ANGPTL4", "MT2A"]
    kt = pd.DataFrame({"pydeseq2_lfc": res["lfc"].reindex(known), "pydeseq2_padj": res["padj"].reindex(known),
                       "published_lfc": p["lfc_dex"].reindex(known), "published_q": p["q_value"].reindex(known)})
    kt.to_csv(TAB / "known_dex_genes.csv")

    # ---- negative control: relabel dex within donors (7 non-trivial swaps) ---------
    perm_rows = []
    ctrl_ids = {cl: g.index[g["dex"] == "no"][0] for cl, g in prim.groupby("cell_line")}
    trt_ids = {cl: g.index[g["dex"] == "yes"][0] for cl, g in prim.groupby("cell_line")}
    donors = sorted(ctrl_ids)
    for flips in product([0, 1], repeat=4):
        if sum(flips) in (0, 4):  # identity and its complement reproduce the real contrast
            continue
        lab = {}
        for cl, f in zip(donors, flips):
            a, b = ctrl_ids[cl], trt_ids[cl]
            lab[a], lab[b] = ("yes", "no") if f else ("no", "yes")
        mp = meta.copy(); mp["dex"] = pd.Categorical([lab[i] for i in mp.index], categories=["no", "yes"])
        rp, _ = run_deseq(counts, mp, "~cell_line + dex", "dex", shrink=False)
        perm_rows.append({"flips": "".join(map(str, flips)), "donors_swapped": int(sum(flips)),
                          "n_sig_padj_lt_0.05": int((rp["padj"] < ALPHA).sum())})
    perm = pd.DataFrame(perm_rows); perm.to_csv(TAB / "permutation_null.csv", index=False)
    # Sign-flip sensitivity check, not a calibrated null: with 2 of 4 donors flipped the paired effect is expected to
    # cancel (6 patterns = 3 complementary pairs); 1 or 3 flipped donors keep part of the real signal.
    metrics["permutation_balanced_signflip_sig_counts"] = perm.loc[perm.donors_swapped == 2, "n_sig_padj_lt_0.05"].tolist()
    metrics["permutation_partial_swap_sig_counts"] = perm.loc[perm.donors_swapped != 2, "n_sig_padj_lt_0.05"].tolist()

    # ---- all 16 samples: label consistency check, then additive model ---------------
    # Dex-response score: mean mean-centred (not variance-scaled) log-expression of the 200 strongest primary-analysis genes
    # (baseMean > 50), signed by their primary log2FC. Alb / Alb+Dex samples were not used to pick the genes.
    top = res[res["baseMean"] > 50].sort_values("padj").head(200)
    lc = np.log2(counts_all.loc[top.index] / (counts_all.sum() / 1e6) + 1)
    lc = lc.sub(lc.mean(axis=1), axis=0)
    score = lc.mul(np.sign(top["lfc"]), axis=0).mean()
    chk = samples[["cell_line", "treatment", "dex"]].assign(dex_response_score=score.round(3))
    chk.to_csv(TAB / "dex_response_score_all16.csv")
    all_meta = samples[["cell_line", "dex", "alb"]].copy()
    swapped = []
    for cl, g in samples.groupby("cell_line"):
        a, ad = g.index[g.treatment == "Albuterol"][0], g.index[g.treatment == "Albuterol_Dexamethasone"][0]
        if score[a] > 0 > score[ad]:  # Alb looks dex-treated, Alb+Dex does not -> labels look swapped
            all_meta.loc[[a, ad], "dex"] = ["yes", "no"]
            swapped.append(cl)
    metrics["apparent_albuterol_label_swaps_donors"] = swapped
    call = counts_all[counts_all.sum(axis=1) >= MIN_TOTAL_COUNT]
    for tag, mt in (("as_labelled", samples[["cell_line", "dex", "alb"]].copy()), ("swap_corrected_posthoc_exploratory", all_meta)):
        mt = mt.copy()
        for c in ("dex", "alb"):
            mt[c] = pd.Categorical(mt[c], categories=["no", "yes"])
        r16, _ = run_deseq(call, mt, "~cell_line + alb + dex", "dex", shrink=False)
        r16.sort_values("padj").to_csv(TAB / f"de_pydeseq2_dex_16samples_{tag}.csv")
        s16 = set(r16.index[r16["padj"] < ALPHA]) & set(common)
        metrics[f"sensitivity_16samples_{tag}"] = {
            "sig": int((r16["padj"] < ALPHA).sum()), "sig_in_common": len(s16),
            "jaccard_with_published": stats.jaccard(A, s16), "intersection": len(A & s16)}

    # ---- sensitivity: unpaired model (no donor term) shows why pairing matters ------
    ru, _ = run_deseq(counts, meta[["dex"]], "~dex", "dex", shrink=False)
    ru.sort_values("padj").to_csv(TAB / "de_pydeseq2_dex_unpaired.csv")
    su = set(ru.index[ru["padj"] < ALPHA]) & set(common)
    metrics["sensitivity_unpaired"] = {"sig": int((ru["padj"] < ALPHA).sum()), "jaccard_with_published": stats.jaccard(A, su)}

    # ---- sensitivity: |log2FC| > 1 on both sides -----------------------------------
    A1 = {g for g in A if abs(merged.loc[g, "log2fc_pub"]) > 1}
    B1 = {g for g in B if abs(merged.loc[g, "lfc_mine"]) > 1}
    metrics["sensitivity_abs_lfc_gt_1"] = {"published": len(A1), "pydeseq2": len(B1), "jaccard": stats.jaccard(A1, B1)}

    # ---- discordant genes -----------------------------------------------------------
    only_pub = merged.loc[list(A - B)].sort_values("q_pub").head(20)
    only_mine = merged.loc[list(B - A)].sort_values("padj_mine").head(20)
    only_pub.to_csv(TAB / "top_published_only.csv"); only_mine.to_csv(TAB / "top_pydeseq2_only.csv")
    metrics["published_only_median_basemean_pydeseq2"] = float(res.loc[list(A - B), "baseMean"].median())
    metrics["shared_median_basemean_pydeseq2"] = float(res.loc[list(A & B), "baseMean"].median())
    metrics["pydeseq2_only_median_basemean"] = float(res.loc[list(B - A), "baseMean"].median())
    metrics["pydeseq2_only_median_abs_lfc"] = float(res.loc[list(B - A), "lfc"].abs().median())
    metrics["shared_median_abs_lfc"] = float(res.loc[list(A & B), "lfc"].abs().median())

    (ROOT / "results" / "metrics.json").write_text(json.dumps(metrics, indent=2))
    topn_s = topn.to_string(index=False, float_format=lambda v: f"{v:.3f}")
    (ROOT / "results" / "overlap_report.txt").write_text(
        "Top-N overlap (published Cuffdiff vs pyDESeq2), common tested genes\n" + topn_s + "\n\n" + json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2)); print(topn_s); print(kt)


if __name__ == "__main__":
    main()
