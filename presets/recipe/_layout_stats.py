"""Layout statistics of the drive5 bench archive (used for docs/step0/data.md)."""
import statistics as st
import sys
from pathlib import Path


def fa(p):
    out, name = [], None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            out.append([l[1:], []]); 
        elif l:
            out[-1][1].append(l)
    return [(n, "".join(s)) for n, s in out]


root = Path(sys.argv[1])
for b in sorted(p.name for p in root.iterdir() if p.is_dir()):
    d = root / b
    ids = (d / "info/ids.txt").read_text().split()
    infos = sorted(x.name for x in (d / "info").iterdir())
    ns, L, nup, nolen, gapless, noup, mism, inlens, lower, other = [], [], [], 0, 0, 0, 0, [], 0, 0
    miss = 0
    for i in ids:
        if not (d / "ref" / i).exists() or not (d / "in" / i).exists():
            miss += 1; continue
        r = fa(d / "ref" / i); q = fa(d / "in" / i)
        ns.append(len(r)); L.append(len(r[0][1])); inlens.append(len(q))
        if len({len(s) for _, s in r}) != 1: nolen += 1
        cols = [c for _, s in r for c in s]
        up = sum(c.isupper() for c in cols); lo = sum(c.islower() for c in cols)
        if up == 0: noup += 1
        if lo: lower += 1
        if not any("-" in s or "." in s for _, s in r): gapless += 1
        if {n for n, _ in r} != {n for n, _ in q}: mism += 1
        if any("-" in s for _, s in q): other += 1
    print(f"{b:11s} ids={len(ids)} missing={miss} info={infos} nseq(ref) min/med/max={min(ns)}/{st.median(ns)}/{max(ns)} "
          f"alnlen min/med/max={min(L)}/{st.median(L)}/{max(L)} in-vs-ref seq count differs={sum(a!=b for a,b in zip(ns,inlens))} "
          f"names differ={mism} ragged={nolen} no-upper={noup} has-lower={lower} gapless-ref={gapless} in-has-gaps={other}")
