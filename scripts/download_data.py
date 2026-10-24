#!/usr/bin/env python3
"""Download all public inputs into data/ and record SHA-256 checksums.

Sources
  * recount3 (Gencode v26) gene-level coverage sums + QC for SRP033351 (= GSE52778)
  * GEO series matrix and the authors' Cuffdiff Dex-vs-Untreated table (published DE result)
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

RC = "https://recount-opendata.s3.amazonaws.com/recount3/release/human"
SRA = f"{RC}/data_sources/sra"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778"
FILES = {
    "gene_sums.tsv.gz": f"{SRA}/gene_sums/51/SRP033351/sra.gene_sums.SRP033351.G026.gz",
    "recount_qc.tsv.gz": f"{SRA}/metadata/51/SRP033351/sra.recount_qc.SRP033351.MD.gz",
    "sra_meta.tsv.gz": f"{SRA}/metadata/51/SRP033351/sra.sra.SRP033351.MD.gz",
    "gene_annotation_G026.gtf.gz": f"{RC}/annotations/gene_sums/human.gene_sums.G026.gtf.gz",
    "GSE52778_series_matrix.txt.gz": f"{GEO}/matrix/GSE52778_series_matrix.txt.gz",
    "GSE52778_Dex_vs_Untreated_gene_exp.diff.gz": f"{GEO}/suppl/GSE52778_Dex_vs_Untreated_gene_exp.diff.gz",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(out: str = "data") -> None:
    out_dir = Path(out)
    out_dir.mkdir(exist_ok=True)
    lines = []
    for name, url in FILES.items():
        dest = out_dir / name
        if not dest.exists():
            print(f"downloading {name}")
            urllib.request.urlretrieve(url, dest)
        lines.append(f"{sha256(dest)}  {name}  {dest.stat().st_size} bytes  {url}")
    (out_dir / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(*sys.argv[1:])
