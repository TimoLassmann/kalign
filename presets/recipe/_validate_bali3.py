"""Validation aid (not part of the pipeline): drive5 case marking vs BAliBASE 3 XML core blocks.

For the 218 BB cases: (1) compare the drive5 core mask (``core_columns``) with the XML column mask
(``<colsco-data>`` of the original BAliBASE_R1-5 release; the reference columns are the same), and
(2) score the SAME kept alignments both ways, drive5 mask on the drive5 reference
(``scoring_protein.score_case``) and XML mask on the original MSF reference (as bench_balibase.py).

    python -m recipe._validate_bali3 --manifest cases.tsv --aln-dir DIR --bb3 bb3_release --out OUT.tsv
"""
from __future__ import annotations

import argparse
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

from .scoring import read_manifest
import tempfile

from .scoring_protein import core_mask, dash_gaps, project, read_fasta, write_fasta


def xml_mask(xml):
    t = ET.parse(xml).getroot().find(".//column-score/colsco-data").text.split()
    return [1 if int(v) == 1 else 0 for v in t]


def main():
    import kalign

    from .scoring_protein import score_case

    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--aln-dir", required=True, help="dir with <case_id>.aln.fa")
    ap.add_argument("--bb3", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = [r for r in read_manifest(a.manifest) if r["benchmark"] == "bali3" and ".BBS" not in r["case_id"]]
    out = []
    for r in rows:
        cid = r["case_id"].split(".", 1)[1]
        rv = f"RV{cid[2:4]}"
        xml, msf = Path(a.bb3, rv, cid + ".xml"), Path(a.bb3, rv, cid + ".msf")
        aln = Path(a.aln_dir) / f"{r['case_id']}.aln.fa"
        xm = xml_mask(xml)
        dm, _ = core_mask(dash_gaps(read_fasta(r["reference"])), "core_columns")
        same_w = len(xm) == len(dm)
        both = sum(1 for x, d in zip(xm, dm) if x and d) if same_w else -1
        new = score_case(r, aln)
        # The drive5 case can hold sequences the release does not (e.g. BB40002: +APP1_SCHPO). The old
        # scorer aligned the release .tfa, so project our alignment onto the release names first.
        rel_names = [n for n, _ in read_fasta(msf.with_suffix(".tfa"))]
        d5_names = [n for n, _ in read_fasta(r["reference"])]
        with tempfile.TemporaryDirectory() as td:
            tf = Path(td) / "t.fasta"
            write_fasta(tf, project(read_fasta(aln), rel_names))
            old = kalign.compare_detailed(str(msf), str(tf), column_mask=xm)
        out.append({"case": cid, "rv": rv, "width_xml": len(xm), "width_d5": len(dm), "n_d5": len(d5_names), "n_release": len(rel_names),
                    "only_d5": ";".join(n for n in d5_names if n not in rel_names),
                    "only_release": ";".join(n for n in rel_names if n not in d5_names), "xml_core": sum(xm),
                    "d5_core": sum(dm), "both": both, "xml_only": (sum(xm) - both) if same_w else "",
                    "d5_only": (sum(dm) - both) if same_w else "",
                    "d5_sp": new["recall"], "d5_p": new["precision"], "d5_f1": new["f1"], "d5_tc": new["tc"],
                    "xml_sp": old["recall"], "xml_p": old["precision"], "xml_f1": old["f1"], "xml_tc": old["tc"]})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    n = len(out)
    print("cases", n, "identical masks", sum(1 for o in out if o["xml_only"] == 0 and o["d5_only"] == 0),
          "width mismatch", sum(1 for o in out if o["width_xml"] != o["width_d5"]))
    for k in ("sp", "p", "f1", "tc"):
        print(k, "d5 %.4f xml %.4f" % (sum(o["d5_" + k] for o in out) / n, sum(o["xml_" + k] for o in out) / n))


if __name__ == "__main__":
    main()
