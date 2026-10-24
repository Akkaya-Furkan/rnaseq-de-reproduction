import gzip

import numpy as np
import pandas as pd
import pytest

from rnaseq_de import data, stats


def test_jaccard_basic():
    assert stats.jaccard({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)
    assert stats.jaccard({"a"}, {"a"}) == 1.0
    assert stats.jaccard({"a"}, {"b"}) == 0.0
    assert np.isnan(stats.jaccard(set(), set()))


def test_top_n_uses_absolute_value_and_is_stable():
    df = pd.DataFrame({"s": [1.0, -5.0, 3.0, -3.0]}, index=list("abcd"))
    assert stats.top_n(df, "s", 2) == ["b", "c"]  # ties keep input order
    assert stats.top_n(df, "s", 10) == ["b", "c", "d", "a"]


def test_top_n_overlap_identical_and_chance():
    idx = [f"g{i}" for i in range(100)]
    x = pd.DataFrame({"s": np.arange(100.0)}, index=idx)
    out = stats.top_n_overlap(x, x, "s", "s", [10])
    assert out.loc[0, "shared"] == 10 and out.loc[0, "jaccard"] == 1.0
    assert out.loc[0, "expected_shared_by_chance"] == pytest.approx(1.0)
    rev = pd.DataFrame({"s": np.arange(100.0)[::-1]}, index=idx)
    assert stats.top_n_overlap(x, rev, "s", "s", [10]).loc[0, "shared"] == 0


def test_top_n_overlap_restricts_to_common_genes():
    a = pd.DataFrame({"s": [3.0, 2.0, 1.0]}, index=["x", "y", "z"])
    b = pd.DataFrame({"s": [3.0, 2.0, 1.0]}, index=["y", "z", "w"])
    out = stats.top_n_overlap(a, b, "s", "s", [1])
    # common = {y, z}; top1 of a -> y, top1 of b -> y
    assert out.loc[0, "shared"] == 1


def test_rank_correlation_and_direction():
    x = pd.Series([1, 2, 3, 4.0], index=list("abcd"))
    assert stats.rank_correlation(x, x * 10)["spearman_rho"] == pytest.approx(1.0)
    assert stats.rank_correlation(x, -x)["spearman_rho"] == pytest.approx(-1.0)
    assert stats.direction_concordance(x - 2.5, pd.Series([-1, -1, 1, 1.0], index=list("abcd"))) == 1.0
    assert stats.direction_concordance(x - 2.5, -(x - 2.5)) == 0.0


def test_coverage_to_counts_scales_per_sample():
    cov = pd.DataFrame({"s1": [100, 250], "s2": [200, 500]}, index=["g1", "g2"])
    out = data.coverage_to_counts(cov, pd.Series({"s1": 100.0, "s2": 50.0}))
    assert out.to_dict("list") == {"s1": [1, 2], "s2": [4, 10]}  # 2.5 rounds to 2 (banker's)


def test_collapse_duplicate_symbols_keeps_highest_total_not_sum():
    counts = pd.DataFrame({"s": [5, 50, 7]}, index=["E1", "E2", "E3"])
    sym = pd.Series({"E1": "A", "E2": "A", "E3": "B"})
    out = data.collapse_duplicate_symbols(counts, sym)
    assert out.loc["A", "s"] == 50 and out.loc["B", "s"] == 7 and len(out) == 2


def test_load_gene_symbols(tmp_path):
    p = tmp_path / "a.gtf.gz"
    rows = ['1\tE\tgene\t1\t2\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "FOO";',
            '1\tE\texon\t1\t2\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "FOO";',
            '1\tE\tgene\t3\t4\t.\t+\t.\tgene_id "ENSG2.1";']
    with gzip.open(p, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    s = data.load_gene_symbols(p)
    assert s.to_dict() == {"ENSG1.1": "FOO", "ENSG2.1": "ENSG2.1"}


def test_load_sample_table(tmp_path):
    p = tmp_path / "m.tsv.gz"
    pd.DataFrame({"external_id": ["SRR1"], "sample_name": ["GSM1"],
                  "sample_attributes": ["cell line;;N1|treatment;;Dexamethasone|tissue;;x"]}).to_csv(p, sep="\t", index=False)
    t = data.load_sample_table(p)
    assert t.loc["SRR1"].to_dict() == {"gsm": "GSM1", "cell_line": "N1", "treatment": "Dexamethasone"}


def test_top_n_overlap_n_larger_than_universe():
    a = pd.DataFrame({"s": [1.0, 2.0]}, index=["x", "y"])
    out = stats.top_n_overlap(a, a, "s", "s", [10])
    assert out.loc[0, "shared"] == 2 and out.loc[0, "expected_shared_by_chance"] == pytest.approx(2.0)
