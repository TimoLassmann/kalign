"""Group related cases (connected components) so that related cases never straddle CV folds.

    cd presets && python -m recipe.grouping --pool protein [--ws ~/tmp/kalign-refit] [--validate N]
    -> <ws>/data/groups_<pool>.tsv        case_id, group_id   (group ids G0001.. by size, then first case id)
       <ws>/data/groups_<pool>.edges.tsv  the links that merged two components (evidence)
       <ws>/data/groups_<pool>.stats.json counts, size distribution, benchmark spans, prefilter recall

Rule (plan section 1, step [2]). Two cases are LINKED if a sequence of one
  (a) is identical to a sequence of the other (hash), or
  (b) is contained in a sequence of the other, or
  (c) shares >= IDENTITY (0.5) identity over the length of the shorter sequence with it.
(b) is the special case identity = 1.0 of (c) and is found by the same machinery.
Groups = connected components (union-find) over cases of ALL benchmarks of a pool together.

Sequences compared: every sequence of the case's unaligned input AND of its reference (the input is a
superset of the reference except for tiny drive5 differences). Upper-case letters only, T=U for the
nucleotide pool. Sequences shorter than MIN_LEN are only used for the exact-hash rule (a short exact
match is still a shared sequence, but short fragments are not allowed to create 50 % links).

Identity (c): semi-global alignment (free end gaps on both sequences, linear gap cost; scores in
PARAMS), identical residues counted on the best-scoring path (ties: more identities), divided by the
length of the SHORTER sequence. X / N / non-ACGT never count as identical. The alignment is numba code
(optional; a slow pure-python fallback is used without numba, fine for tests only).

Speed: unique sequences only; exact hashes first; then a k-mer prefilter (shared distinct k-mers via an
inverted index, k and the minimum shared count in PARAMS, k-mers occurring in more than MAX_POST
sequences ignored) proposes pairs; pairs whose cases are already in the same component are dropped,
candidates are verified best-first in rounds, and components are re-checked between rounds so that
big families are not verified exhaustively. The prefilter is not exact: `--validate N` aligns all pairs
of N random sequences and reports the fraction of true links (identity >= IDENTITY) the prefilter
had proposed. Since groups are connected components, a missed pair only matters if it is the only
bridge between two components.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from .scoring import read_manifest
from .scoring_protein import read_fasta

try:
    import numba
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover - exercised on machines without numba
    HAVE_NUMBA = False

    def njit(*a, **k):
        if a and callable(a[0]):
            return a[0]
        return lambda f: f

    prange = range

IDENTITY = 0.5
ALPHABETS = {"protein": "ACDEFGHIKLMNPQRSTVWY", "nucleotide": "ACGT"}
PARAMS = {
    "protein": dict(k=5, min_shared=3, min_len=30, match=2, mismatch=-1, gap=-3),
    "nucleotide": dict(k=7, min_shared=4, min_len=30, match=2, mismatch=-3, gap=-5),
}
# Length-dependent identity needed to link (shorter length of the pair < upper bound -> threshold).
# Protein: a shuffle null (2,000 sequences, 2M pairs) never reached 0.38 even at 30-40 residues, so the
# plain 0.5 rule is used. Nucleotides: with 4 letters, short sequences reach 50 % by chance (shuffle null
# 2M pairs: shorter length 30-40 max 0.90, 40-60 max 0.83, 60-100 max 0.75, >=100 max 0.50 with 0 pairs
# >= 0.5) so the rule is tightened below 100 nt. See docs/step0/features.md.
IDENT_SCHEDULE = {"protein": [], "nucleotide": [(50, 0.95), (100, 0.80)]}
MAX_POST = 5000


# ----------------------------------------------------------------------------- sequences
def normalise(seq: str, pool: str) -> str:
    s = "".join(c for c in seq.upper() if "A" <= c <= "Z")
    return s.replace("U", "T") if pool == "nucleotide" else s


def encode(seq: str, pool: str) -> np.ndarray:
    """uint8 codes: 1.. for the pool alphabet, 0 for anything else (never matches)."""
    lut = np.zeros(256, dtype=np.uint8)
    for i, c in enumerate(ALPHABETS[pool]):
        lut[ord(c)] = i + 1
    return lut[np.frombuffer(seq.encode("latin-1"), dtype=np.uint8)]


def load_sequences(rows, pool, source="both"):
    """-> (uid_seqs list[str], uid_cases list[list[case index]]) over unique normalised sequences."""
    seen: dict[str, int] = {}
    cases: list[list[int]] = []
    seqs: list[str] = []
    for ci, row in enumerate(rows):
        mine = set()
        paths = {"both": (row["unaligned"], row["reference"]), "input": (row["unaligned"],),
                 "reference": (row["reference"],)}[source]
        for path in paths:
            for _, s in read_fasta(path):
                s = normalise(s, pool)
                if s:
                    mine.add(s)
        for s in sorted(mine):
            u = seen.get(s)
            if u is None:
                u = seen[s] = len(seqs)
                seqs.append(s)
                cases.append([])
            cases[u].append(ci)
    return seqs, cases


# ----------------------------------------------------------------------------- identity
@njit(cache=True, nogil=True)
def sg_identity(a, b, match, mismatch, gap):
    """Semi-global identity over the shorter sequence: (identical residues on best path) / min(len)."""
    if a.shape[0] > b.shape[0]:
        a, b = b, a
    la, lb = a.shape[0], b.shape[0]
    if la == 0:
        return 0.0
    S = np.zeros(lb + 1, dtype=np.int32)
    M = np.zeros(lb + 1, dtype=np.int32)
    best_s, best_m = -1000000, 0
    for i in range(1, la + 1):
        ai = a[i - 1]
        diag_s, diag_m = S[0], M[0]
        S[0] = 0
        M[0] = 0
        for j in range(1, lb + 1):
            up_s, up_m = S[j], M[j]
            eq = 1 if (ai != 0 and ai == b[j - 1]) else 0
            s1 = diag_s + (match if eq == 1 else mismatch)
            m1 = diag_m + eq
            s2 = up_s + gap
            m2 = up_m
            s3 = S[j - 1] + gap
            m3 = M[j - 1]
            bs, bm = s1, m1
            if s2 > bs or (s2 == bs and m2 > bm):
                bs, bm = s2, m2
            if s3 > bs or (s3 == bs and m3 > bm):
                bs, bm = s3, m3
            S[j] = bs
            M[j] = bm
            diag_s, diag_m = up_s, up_m
        if S[lb] > best_s or (S[lb] == best_s and M[lb] > best_m):
            best_s, best_m = S[lb], M[lb]
    for j in range(lb + 1):
        if S[j] > best_s or (S[j] == best_s and M[j] > best_m):
            best_s, best_m = S[j], M[j]
    return best_m / la


@njit(cache=True, nogil=True, parallel=True)
def verify_pairs(codes, offs, I, J, match, mismatch, gap):
    out = np.zeros(I.shape[0], dtype=np.float64)
    for p in prange(I.shape[0]):
        i, j = I[p], J[p]
        out[p] = sg_identity(codes[offs[i]:offs[i + 1]], codes[offs[j]:offs[j + 1]], match, mismatch, gap)
    return out


def required_identity(pool, shortlen, base=None):
    """Identity needed to link a pair whose shorter sequence has length shortlen (array or int)."""
    base = IDENTITY if base is None else base
    sl = np.asarray(shortlen)
    req = np.full(sl.shape, base, dtype=np.float64)
    for bound, thr in reversed(IDENT_SCHEDULE[pool]):
        req = np.where(sl < bound, max(thr, base), req)
    return req


# ----------------------------------------------------------------------------- union-find
@njit(cache=True)
def uf_find(parent, x):
    r = x
    while parent[r] != r:
        r = parent[r]
    while parent[x] != r:
        nx = parent[x]
        parent[x] = r
        x = nx
    return r


@njit(cache=True)
def uf_union(parent, a, b):
    ra, rb = uf_find(parent, a), uf_find(parent, b)
    if ra == rb:
        return False
    if ra < rb:
        parent[rb] = ra
    else:
        parent[ra] = rb
    return True


@njit(cache=True)
def uf_roots(parent):
    out = np.empty(parent.shape[0], dtype=np.int64)
    for i in range(parent.shape[0]):
        out[i] = uf_find(parent, i)
    return out


# ----------------------------------------------------------------------------- k-mer prefilter
@njit(cache=True)
def _kmers_one(c, k, base, out, pos):
    """Append k-mer integers of code array c (codes 1..base) to out starting at pos; skips k-mers with code 0."""
    run, val = 0, 0
    top = base ** (k - 1)
    for t in range(c.shape[0]):
        x = np.int64(c[t])
        if x == 0:
            run, val = 0, 0
            continue
        if run == k:
            val = val - (val // top) * top  # drop oldest digit
        else:
            run += 1
        val = val * base + (x - 1)
        if run == k:
            out[pos] = val
            pos += 1
    return pos


@njit(cache=True)
def all_kmers(codes, offs, k, base):
    n = offs.shape[0] - 1
    total = 0
    for u in range(n):
        total += max(0, offs[u + 1] - offs[u] - k + 1)
    out_k = np.empty(total, dtype=np.int64)
    out_u = np.empty(total, dtype=np.int64)
    pos = 0
    for u in range(n):
        start = pos
        pos = _kmers_one(codes[offs[u]:offs[u + 1]], k, base, out_k, pos)
        for t in range(start, pos):
            out_u[t] = u
    return out_k[:pos], out_u[:pos]


def build_index(codes, offs, k, base):
    """Distinct (uid, kmer) pairs as two CSR structures: kmer -> uids and uid -> kmers."""
    n = len(offs) - 1
    kk, uu = all_kmers(codes, offs, k, base)
    key = np.unique(kk * n + uu)  # distinct pairs, sorted by kmer then uid
    km, ui = key // n, key % n
    kvals, kstart = np.unique(km, return_index=True)
    post_start = np.append(kstart, len(km)).astype(np.int64)
    # uid -> kmer index (position in kvals)
    kidx = np.searchsorted(kvals, km)
    order = np.argsort(ui, kind="stable")
    u_sorted, k_sorted = ui[order], kidx[order].astype(np.int64)
    ustart = np.searchsorted(u_sorted, np.arange(n + 1)).astype(np.int64)
    return dict(post_start=post_start, post_uid=ui.astype(np.int64), uid_start=ustart, uid_kidx=k_sorted,
                nk=np.diff(ustart).astype(np.int64))


@njit(cache=True, nogil=True)
def cand_chunk(lo, hi, uid_start, uid_kidx, post_start, post_uid, label, nk, min_shared, max_post, n):
    counts = np.zeros(n, dtype=np.int32)
    touched = np.empty(n, dtype=np.int64)
    oi, oj, oc = [], [], []
    for i in range(lo, hi):
        nt = 0
        for q in range(uid_start[i], uid_start[i + 1]):
            kx = uid_kidx[q]
            a, b = post_start[kx], post_start[kx + 1]
            if b - a > max_post:
                continue
            for p in range(a, b):
                j = post_uid[p]
                if j <= i or label[j] == label[i]:
                    continue
                if counts[j] == 0:
                    touched[nt] = j
                    nt += 1
                counts[j] += 1
        for t in range(nt):
            j = touched[t]
            need = min(min_shared, max(1, min(nk[i], nk[j]) // 4))
            if counts[j] >= need:
                oi.append(i)
                oj.append(j)
                oc.append(counts[j])
            counts[j] = 0
    return np.array(oi, dtype=np.int64), np.array(oj, dtype=np.int64), np.array(oc, dtype=np.int64)


def candidates(index, label, params, n, threads):
    chunks = max(1, min(n, threads * 8))
    bounds = np.linspace(0, n, chunks + 1).astype(np.int64)

    def run(c):
        return cand_chunk(bounds[c], bounds[c + 1], index["uid_start"], index["uid_kidx"], index["post_start"],
                          index["post_uid"], label, index["nk"], params["min_shared"], MAX_POST, n)

    with ThreadPoolExecutor(threads) as ex:
        parts = list(ex.map(run, range(chunks)))
    I = np.concatenate([p[0] for p in parts]) if parts else np.zeros(0, np.int64)
    J = np.concatenate([p[1] for p in parts]) if parts else np.zeros(0, np.int64)
    C = np.concatenate([p[2] for p in parts]) if parts else np.zeros(0, np.int64)
    return I, J, C


# ----------------------------------------------------------------------------- main routine
def group_cases(rows, pool, threads=16, rnd=50000, validate=0, seed=1, log=print, source="both", identity=None):
    """-> (group label per case index (list of root case ids), edges, stats)."""
    global IDENTITY
    if identity is not None:
        IDENTITY = identity
    p = PARAMS[pool]
    base = len(ALPHABETS[pool])
    t0 = time.time()
    seqs, ucases = load_sequences(rows, pool, source)
    ncase, nuid = len(rows), len(seqs)
    log(f"[{time.time()-t0:6.1f}s] {ncase} cases, {nuid} unique sequences "
        f"({sum(len(s) for s in seqs)} residues)")
    parent = np.arange(ncase, dtype=np.int64)
    edges = []  # (case_a, case_b, kind, identity)
    # (a) exact: all cases sharing a sequence
    for u, cs in enumerate(ucases):
        for c in cs[1:]:
            if uf_union(parent, cs[0], c):
                edges.append((cs[0], c, "exact", 1.0))
    log(f"[{time.time()-t0:6.1f}s] exact links: {len(edges)} merges, {len(set(uf_roots(parent).tolist()))} components")
    # (b,c) prefilter + verify, only sequences long enough
    usable = np.array([len(s) >= p["min_len"] for s in seqs])
    enc = [encode(s, pool) for s in seqs]
    lens = np.array([len(e) for e in enc], dtype=np.int64)
    # sequences below min_len get zero length in the index so they never produce candidates
    codes_list = [e if usable[i] else e[:0] for i, e in enumerate(enc)]
    offs = np.zeros(nuid + 1, dtype=np.int64)
    offs[1:] = np.cumsum([len(e) for e in codes_list])
    codes = np.concatenate(codes_list).astype(np.uint8) if nuid else np.zeros(0, np.uint8)
    rep = np.array([cs[0] for cs in ucases], dtype=np.int64)
    index = build_index(codes, offs, p["k"], base)
    label = uf_roots(parent)[rep]
    I, J, C = candidates(index, label, p, nuid, threads)
    log(f"[{time.time()-t0:6.1f}s] prefilter: {len(I)} candidate pairs "
        f"({int((np.diff(index['post_start']) > MAX_POST).sum())} k-mers over MAX_POST ignored)")
    cand_keys = I * nuid + J
    order = np.lexsort((J, I, -C))
    I, J = I[order], J[order]
    verified = n_pass = 0
    for s in range(0, len(I), rnd):
        ii, jj = I[s:s + rnd], J[s:s + rnd]
        roots = uf_roots(parent)
        keep = roots[rep[ii]] != roots[rep[jj]]
        ii, jj = ii[keep], jj[keep]
        if len(ii) == 0:
            continue
        # verify on the original (unshortened) sequences
        ident = verify_pairs(allcodes(enc), allofs(enc), ii, jj, p["match"], p["mismatch"], p["gap"])
        verified += len(ii)
        req = required_identity(pool, np.minimum(lens[ii], lens[jj]))
        for a, b, idn, rq in zip(ii.tolist(), jj.tolist(), ident.tolist(), req.tolist()):
            if idn >= rq:
                n_pass += 1
                if uf_union(parent, rep[a], rep[b]):
                    edges.append((int(rep[a]), int(rep[b]), "identity" if idn < 0.9999 else "contained/identical", idn))
        log(f"[{time.time()-t0:6.1f}s] verified {verified} / {len(I)} candidates, {n_pass} links, "
            f"{len(set(uf_roots(parent).tolist()))} components")
    stats = dict(n_cases=ncase, n_unique_seqs=nuid, n_usable_seqs=int(usable.sum()), candidates=int(len(I)),
                 verified=int(verified), links=int(n_pass), params=p, identity_threshold=IDENTITY,
                 max_post=MAX_POST, numba=HAVE_NUMBA)
    if validate:
        stats["prefilter_recall"] = validate_prefilter(enc, usable, label, cand_keys, nuid, p, validate, seed, log, pool)
    return uf_roots(parent), edges, stats


_ALL_CACHE: dict = {}


def allcodes(enc):
    k = id(enc)
    if k not in _ALL_CACHE:
        _ALL_CACHE.clear()
        _ALL_CACHE[k] = (np.concatenate(enc).astype(np.uint8) if enc else np.zeros(0, np.uint8),
                         np.concatenate([[0], np.cumsum([len(e) for e in enc])]).astype(np.int64))
    return _ALL_CACHE[k][0]


def allofs(enc):
    allcodes(enc)
    return _ALL_CACHE[id(enc)][1]


def validate_prefilter(enc, usable, label, cand_keys, nuid, p, nsample, seed, log, pool):
    """All-pairs identity on a random sample of usable sequences; fraction of true links (in different
    initial components) that the prefilter had proposed."""
    rng = np.random.default_rng(seed)
    pool_ids = np.flatnonzero(usable)
    pick = np.sort(rng.choice(pool_ids, size=min(nsample, len(pool_ids)), replace=False))
    ii, jj = np.triu_indices(len(pick), 1)
    I, J = pick[ii], pick[jj]
    keep = label[I] != label[J]
    I, J = I[keep], J[keep]
    ident = verify_pairs(allcodes(enc), allofs(enc), I, J, p["match"], p["mismatch"], p["gap"])
    lens = np.array([len(e) for e in enc])
    true = ident >= required_identity(pool, np.minimum(lens[I], lens[J]))
    keys = I[true] * nuid + J[true]
    found = np.isin(keys, cand_keys)
    res = dict(sample=int(len(pick)), pairs=int(len(I)), true_links=int(true.sum()),
               proposed=int(found.sum()), recall=float(found.mean()) if len(keys) else None)
    # by identity band
    bands = {}
    idt = ident[true]
    for lo, hi in ((0.5, 0.6), (0.6, 0.8), (0.8, 1.01)):
        sel = (idt >= lo) & (idt < hi)
        bands[f"{lo}-{min(hi,1.0)}"] = dict(n=int(sel.sum()), recall=float(found[sel].mean()) if sel.any() else None)
    res["by_identity"] = bands
    log(f"prefilter validation: {res}")
    return res


# ----------------------------------------------------------------------------- output
def relabel(roots, case_ids):
    """Group ids G0001.. ordered by size (desc) then first case id; returns list parallel to case_ids."""
    members = defaultdict(list)
    for i, r in enumerate(roots.tolist()):
        members[r].append(i)
    order = sorted(members.values(), key=lambda m: (-len(m), min(case_ids[i] for i in m)))
    gid = {}
    for n, m in enumerate(order, 1):
        for i in m:
            gid[i] = f"G{n:04d}"
    return [gid[i] for i in range(len(case_ids))]


def summarise(rows, gids, edges):
    members = defaultdict(list)
    for r, g in zip(rows, gids):
        members[g].append(r)
    sizes = sorted((len(m) for m in members.values()), reverse=True)
    dist = Counter(sizes)
    bucket = lambda s: "1" if s == 1 else "2" if s == 2 else "3-5" if s <= 5 else "6-20" if s <= 20 else "21-100" if s <= 100 else ">100"
    bd = Counter(bucket(s) for s in sizes)
    largest = []
    for g, m in sorted(members.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:15]:
        largest.append(dict(group=g, n_cases=len(m), benchmarks=dict(Counter(r["benchmark"] for r in m))))
    multi_bench = [g for g, m in members.items() if len({r["benchmark"] for r in m}) > 1]
    pair_span = Counter()
    for g in multi_bench:
        pair_span[" + ".join(sorted({r["benchmark"] for r in members[g]}))] += 1
    pair_cases = Counter()
    for g in multi_bench:
        for b in {r["benchmark"] for r in members[g]}:
            pair_cases[b] += sum(1 for r in members[g] if r["benchmark"] == b)
    # cross-benchmark merging edges
    bench = {i: r["benchmark"] for i, r in enumerate(rows)}
    cross = Counter()
    for a, b, kind, idn in edges:
        if bench[a] != bench[b]:
            cross[" - ".join(sorted((bench[a], bench[b])))] += 1
    per_bench = {}
    for bname in sorted({r["benchmark"] for r in rows}):
        rs = [r for r in rows]
        gs = {g for r, g in zip(rows, gids) if r["benchmark"] == bname}
        alone = sum(1 for g in gs if all(x["benchmark"] == bname for x in members[g]))
        per_bench[bname] = dict(cases=sum(1 for r in rows if r["benchmark"] == bname), groups_touched=len(gs),
                                groups_only_this_benchmark=alone,
                                cases_in_multi_benchmark_groups=pair_cases.get(bname, 0))
    return dict(n_groups=len(members), n_singletons=bd["1"], size_buckets=dict(bd),
                size_histogram={str(k): v for k, v in sorted(dist.items())}, max_size=sizes[0],
                largest_groups=largest, multi_benchmark_groups=len(multi_bench),
                multi_benchmark_group_types=dict(pair_span), merging_edges_by_benchmark_pair=dict(cross),
                per_benchmark=per_bench)


# ----------------------------------------------------------------------------- final (paper) grouping
def merge_extra(roots, pairs, ids):
    """Union extra links on top of component roots. roots: array (root case index per case); pairs: list of
    (case index, case index, kind, identity). -> (new roots array, list of links that merged two components)."""
    parent = {int(r): int(r) for r in set(roots.tolist())}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    merged = []
    for a, b, kind, idn in pairs:
        ra, rb = find(int(roots[a])), find(int(roots[b]))
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
            merged.append((a, b, kind, idn))
    return np.array([find(int(r)) for r in roots.tolist()], dtype=np.int64), merged


def clan_pairs(rows, clans):
    """Links between Rfam cases whose families are in one clan (`#=GF CL`); star topology per clan."""
    by_clan = defaultdict(list)
    for i, r in enumerate(rows):
        fam = r.get("family") if r["benchmark"] == "rfam" else None
        if fam in clans:
            by_clan[clans[fam]].append(i)
    return [(m[0], x, "clan", 1.0) for m in by_clan.values() for x in m[1:]]


def containment_pairs(rows, pool, source="reference", min_len=None, k=12):
    """Exact substring containment between sequences of different cases (letters compared literally, so X / B /
    Z count as characters). The identity machinery (sg_identity) never counts X / B / Z as identical, so a
    sequence that is contained in another but holds an ambiguous residue scores < 1.0 there and is missed;
    this pass closes that gap (review 2026-10-05: 4 such cross-fold pairs, e.g. 1ubpC in two PREFAB cases).
    -> [(case index, case index, "contained-exact", 1.0)] (one representative case per sequence)."""
    min_len = PARAMS[pool]["min_len"] if min_len is None else min_len
    seqs, ucases = load_sequences(rows, pool, source)
    idx = defaultdict(set)
    for i, s in enumerate(seqs):
        if len(s) >= min_len:
            for j in range(len(s) - k + 1):
                idx[s[j:j + k]].add(i)
    out = []
    for i, a in enumerate(seqs):
        if len(a) < min_len:
            continue
        for j in idx.get(a[:k], ()):
            b = seqs[j]
            if j != i and len(b) >= len(a) and (len(b) > len(a) or j > i) and a in b:
                out.append((ucases[i][0], ucases[j][0], "contained-exact", 1.0))
    return out


def final_groups(rows, pool, threads=16, log=print, clans=None):
    """The groups used for the folds of the final fit.
    protein:     a REFERENCE sequence of one case is identical to / contained in a reference sequence of the other.
    nucleotide:  the step-0 rule (input + reference sequences, length-dependent identity) plus Rfam clan links."""
    if pool == "protein":
        roots, edges, stats = group_cases(rows, pool, threads=threads, log=log, source="reference", identity=1.0)
        cp = containment_pairs(rows, pool, "reference")
        roots, merged = merge_extra(roots, cp, None)
        edges = list(edges) + merged
        stats["containment_exact_pairs"] = len(cp)
        stats["containment_exact_merges"] = len(merged)
        log(f"exact containment: {len(cp)} pairs ({len(merged)} merged components); {len(set(roots.tolist()))} groups")
        stats["rule"] = ("reference sequences identical or contained (identity 1.0 over the shorter, min length 30; "
                         "plus exact substring containment, so ambiguous residues X/B/Z cannot hide a containment)")
    else:
        roots, edges, stats = group_cases(rows, pool, threads=threads, log=log, source="both", identity=0.5)
        stats["rule"] = "step-0 length-dependent identity rule on input+reference sequences, plus Rfam clans"
        if clans:
            cp = clan_pairs(rows, clans)
            roots, merged = merge_extra(roots, cp, None)
            edges = list(edges) + merged
            stats["clan_links"] = len(cp)
            stats["clan_merges"] = len(merged)
            log(f"clan links: {len(cp)} ({len(merged)} merged components); {len(set(roots.tolist()))} groups")
    return roots, edges, stats


def homolog_flags(rows, pool, fold_of, threads=16, log=print):
    """homolog_in_train covariate: True for a case if some sequence of it (input or reference) has a link by the
    step-0 identity rule (protein: >= 0.5 identity over the shorter; nucleotide: length-dependent) to a sequence
    of a case in ANOTHER fold. fold_of: list of fold ints parallel to rows. -> (np.bool_ array, info dict)."""
    p = PARAMS[pool]
    base = len(ALPHABETS[pool])
    seqs, ucases = load_sequences(rows, pool, "both")
    nuid = len(seqs)
    fold_of = np.asarray(fold_of, dtype=np.int64)
    umask = np.array([sum(1 << int(f) for f in set(fold_of[cs].tolist())) for cs in ucases], dtype=np.int64)
    flags = np.zeros(len(rows), dtype=bool)
    exact = 0
    for u, cs in enumerate(ucases):  # a sequence shared by cases of different folds
        m = int(umask[u])
        if m & (m - 1):
            for c in cs:
                flags[c] = True
            exact += 1
    log(f"homolog flags: {nuid} unique sequences, {exact} shared across folds, {int(flags.sum())} cases flagged (exact)")
    usable = np.array([len(x) >= p["min_len"] for x in seqs])
    enc = [encode(x, pool) for x in seqs]
    lens = np.array([len(e) for e in enc], dtype=np.int64)
    codes_list = [e if usable[i] else e[:0] for i, e in enumerate(enc)]
    offs = np.zeros(nuid + 1, dtype=np.int64)
    offs[1:] = np.cumsum([len(e) for e in codes_list])
    codes = np.concatenate(codes_list).astype(np.uint8)
    index = build_index(codes, offs, p["k"], base)
    # sequences of one fold only share a label (their pairs are same-fold, irrelevant); multi-fold sequences are unique
    single = (umask & (umask - 1)) == 0
    label = np.where(single, np.log2(np.maximum(umask, 1)).astype(np.int64), 1000 + np.arange(nuid))
    I, J, C = candidates(index, label, p, nuid, threads)
    order = np.lexsort((J, I, -C))
    I, J = I[order], J[order]
    log(f"homolog flags: {len(I)} candidate pairs")
    nun = np.array([len(cs) for cs in ucases], dtype=np.int64)  # unflagged case count per sequence (upper bound)
    ac, ao = allcodes(enc), allofs(enc)
    verified = links = 0
    for st in range(0, len(I), 50000):
        ii, jj = I[st:st + 50000], J[st:st + 50000]
        # skip pairs whose cases are all flagged already
        keep = np.array([not (all(flags[c] for c in ucases[a]) and all(flags[c] for c in ucases[b]))
                         for a, b in zip(ii.tolist(), jj.tolist())], dtype=bool)
        ii, jj = ii[keep], jj[keep]
        if not len(ii):
            continue
        ident = verify_pairs(ac, ao, ii, jj, p["match"], p["mismatch"], p["gap"])
        req = required_identity(pool, np.minimum(lens[ii], lens[jj]), base=0.5)
        verified += len(ii)
        for a, b, idn, rq in zip(ii.tolist(), jj.tolist(), ident.tolist(), req.tolist()):
            if idn < rq:
                continue
            links += 1
            for u, v in ((a, b), (b, a)):
                for c in ucases[u]:
                    if int(umask[v]) & ~(1 << int(fold_of[c])):
                        flags[c] = True
        if (st // 50000) % 20 == 0:
            log(f"  verified {verified}, links {links}, flagged {int(flags.sum())}")
    info = dict(unique_seqs=nuid, shared_across_folds=exact, candidates=int(len(I)), verified=int(verified),
                links=int(links), flagged=int(flags.sum()), prefilter=dict(p))
    return flags, info


def main_final(a):
    from .cases_nucleotide import rfam_clans
    if HAVE_NUMBA:
        numba.set_num_threads(a.threads)
    ws = Path(a.ws).expanduser()
    rows = read_manifest(ws / "data" / f"cases_{a.pool}.final.tsv")
    ids = [r["case_id"] for r in rows]
    clans = rfam_clans(ws / "data" / "raw" / "rfam" / "Rfam.seed") if a.pool == "nucleotide" else None
    for r in rows:  # family accession needed for clan links (nucleotide manifest column)
        r.setdefault("family", "")
    roots, edges, stats = final_groups(rows, a.pool, threads=a.threads, clans=clans)
    gids = relabel(roots, ids)
    d = ws / "data"
    with open(d / f"groups_{a.pool}.final.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["case_id", "group_id"])
        w.writerows(zip(ids, gids))
    with open(d / f"groups_{a.pool}.final.edges.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["case_a", "case_b", "benchmark_a", "benchmark_b", "kind", "identity"])
        for x, y, kind, idn in edges:
            w.writerow([ids[x], ids[y], rows[x]["benchmark"], rows[y]["benchmark"], kind, f"{idn:.4f}"])
    stats.update(summarise(rows, gids, edges))
    stats["largest_group_share"] = stats["max_size"] / len(rows)
    (d / f"groups_{a.pool}.final.stats.json").write_text(json.dumps(stats, indent=1))
    print(json.dumps({k: v for k, v in stats.items() if k not in ("params",)}, indent=1))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True, choices=["protein", "nucleotide"])
    ap.add_argument("--ws", default=os.environ.get("KALIGN_REFIT_WS", "~/tmp/kalign-refit"))
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--seqs", default="both", choices=["both", "input", "reference"],
                    help="which sequences of a case are compared (default both = input + reference; the official groups)")
    ap.add_argument("--identity", type=float, default=None, help="identity threshold (default 0.5; sensitivity runs only)")
    ap.add_argument("--out-dir", default=None, help="default <ws>/data; use a scratch dir for sensitivity runs")
    ap.add_argument("--validate", type=int, default=0, help="all-pairs prefilter recall check on N random sequences")
    ap.add_argument("--final", action="store_true",
                    help="final grouping of the paper fit: cases_<pool>.final.tsv -> groups_<pool>.final.tsv")
    a = ap.parse_args(argv)
    if a.final:
        return main_final(a)
    if HAVE_NUMBA:
        numba.set_num_threads(a.threads)
    ws = Path(a.ws).expanduser()
    rows = read_manifest(ws / "data" / f"cases_{a.pool}.tsv")
    ids = [r["case_id"] for r in rows]
    roots, edges, stats = group_cases(rows, a.pool, threads=a.threads, validate=a.validate, source=a.seqs, identity=a.identity)
    gids = relabel(roots, ids)
    d = Path(a.out_dir).expanduser() if a.out_dir else ws / "data"
    d.mkdir(parents=True, exist_ok=True)
    with open(d / f"groups_{a.pool}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["case_id", "group_id"])
        w.writerows(zip(ids, gids))
    with open(d / f"groups_{a.pool}.edges.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["case_a", "case_b", "benchmark_a", "benchmark_b", "kind", "identity"])
        for x, y, kind, idn in edges:
            w.writerow([ids[x], ids[y], rows[x]["benchmark"], rows[y]["benchmark"], kind, f"{idn:.4f}"])
    stats.update(summarise(rows, gids, edges))
    stats["sequence_source"] = a.seqs
    stats["identity_threshold"] = IDENTITY
    (d / f"groups_{a.pool}.stats.json").write_text(json.dumps(stats, indent=1))
    print(json.dumps({k: v for k, v in stats.items() if k not in ("params",)}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
