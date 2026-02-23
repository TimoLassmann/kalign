#!/usr/bin/env python3
"""Compare gap statistics across alignment methods at varying divergence."""

import os
import subprocess
import kalign


def read_aln(path):
    seqs = []
    with open(path) as f:
        cur = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur:
                    seqs.append(cur)
                cur = ""
            else:
                cur += line
        if cur:
            seqs.append(cur)
    return seqs


def gap_opens_per_seq(aln_seqs):
    n = len(aln_seqs)
    alnlen = len(aln_seqs[0])
    opens = 0
    for s in aln_seqs:
        in_gap = False
        for j in range(alnlen):
            if s[j] == "-":
                if not in_gap:
                    opens += 1
                    in_gap = True
            else:
                in_gap = False
    return opens / n


cases = [
    "WAG_t16_d0.5_ir0.1_il2.0_r0",
    "WAG_t16_d1.0_ir0.1_il2.0_r0",
    "WAG_t16_d2.0_ir0.1_il2.0_r0",
    "WAG_t16_d4.0_ir0.1_il2.0_r0",
    "WAG_t16_d4.0_ir0.2_il2.0_r0",
    "WAG_t16_d4.0_ir0.2_il5.0_r0",
]

hdr = (
    f"{'case':40s} {'depth':>5s} "
    f"{'true':>6s} {'t_go':>6s} | "
    f"{'kal':>6s} {'k_go':>6s} "
    f"{'mft':>6s} {'m_go':>6s} "
    f"{'mus':>6s} {'mu_go':>6s} "
    f"{'clu':>6s} {'c_go':>6s}"
)
print(hdr)
print("-" * len(hdr))

for case in cases:
    base = f"benchmarks/data/{case}"
    fasta = f"{base}/unaligned.fasta"
    true_path = f"{base}/true_alignment.fasta"
    if not os.path.exists(fasta) or not os.path.exists(true_path):
        print(f"{case:40s} MISSING")
        continue

    depth = case.split("_d")[1].split("_")[0]

    true = read_aln(true_path)
    t_len = len(true[0])
    t_go = gap_opens_per_seq(true)

    # Kalign
    r = kalign.align_from_file(fasta, seq_type="protein")
    k_len = len(r.sequences[0])
    k_go = gap_opens_per_seq(list(r.sequences))

    # MAFFT
    res = subprocess.run(
        ["mafft", "--quiet", fasta], capture_output=True, text=True
    )
    with open("/tmp/mafft_out.fa", "w") as f:
        f.write(res.stdout)
    m_aln = read_aln("/tmp/mafft_out.fa")
    m_len = len(m_aln[0])
    m_go = gap_opens_per_seq(m_aln)

    # MUSCLE
    subprocess.run(
        ["muscle", "-align", fasta, "-output", "/tmp/muscle_out.fa", "-quiet"],
        capture_output=True,
    )
    mu_aln = read_aln("/tmp/muscle_out.fa")
    mu_len = len(mu_aln[0])
    mu_go = gap_opens_per_seq(mu_aln)

    # Clustal Omega
    subprocess.run(
        ["clustalo", "-i", fasta, "-o", "/tmp/clustal_out.fa", "--force"],
        capture_output=True,
    )
    c_aln = read_aln("/tmp/clustal_out.fa")
    c_len = len(c_aln[0])
    c_go = gap_opens_per_seq(c_aln)

    print(
        f"{case:40s} {depth:>5s} "
        f"{t_len:>6d} {t_go:>6.1f} | "
        f"{k_len:>6d} {k_go:>6.1f} "
        f"{m_len:>6d} {m_go:>6.1f} "
        f"{mu_len:>6d} {mu_go:>6.1f} "
        f"{c_len:>6d} {c_go:>6.1f}"
    )

print()
print("Columns: true=true_alnlen, t_go=true_gap_opens/seq")
print("         kal/mft/mus/clu = alignment length, _go = gap_opens/seq")
