"""Build presets/results/v2/report.html (andamentum report dialect) from the committed result files.

    uv run python presets/results/v2/report_build/build.py

Every number on the page is read from a file under presets/results/ (tables) or written in body.html next to the
file it comes from; every chart is a rendered SVG from recipe.report_figures, inlined unchanged.
"""

import csv
import html
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
RES = V2.parent
SKILL = Path("~/code/claude-plugins/document-tools/skills/html-report/assets").expanduser()
CORE_N = {"protein": {f: 840 for f in ("0", "1", "2", "3", "4", "final")},
          "nucleotide": {"0": 2201, "1": 1938, "2": 1934, "3": 1931, "4": 1948, "final": 2488}}  # PLAN.md / report 5b


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def table(head, rows, label_col=True):
    h = "".join(f"<th>{html.escape(c)}</th>" for c in head)
    b = ""
    for r in rows:
        cells = "".join((f"<th>{c}</th>" if (i == 0 and label_col) else f"<td>{c}</td>") for i, c in enumerate(r))
        b += f"<tr>{cells}</tr>"
    return f'<div class="scroll"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def bold_new(name):
    return f"<strong>{name}</strong>" if name.startswith("Kalign new") else name


# ------------------------------------------------------------------ tables from the A2 TABLES.md (one program per row)
def md_tables(path):
    out, cur, title = {}, None, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("### "):
            title = line[4:].strip()
            cur = out.setdefault(title, [])
        elif line.startswith("|") and cur is not None and not set(line) <= set("|-"):
            cur.append([c.strip() for c in line.strip("|").split("|")])
    return out


T = md_tables(V2 / "addendum_a2" / "TABLES.md")


def md_to_html(title, keep=None):
    rows = T[title]
    head, body = rows[0], rows[1:]
    if keep:
        body = [r for r in body if r[0] in keep]
    return table(head, [[bold_new(r[0])] + r[1:] for r in body])


KEEP_MAIN = ["Kalign 3.1", "Kalign 3.3.5", "Kalign 3.6.1 recall", "Kalign 3.7.0 fast", "Kalign 3.7.0 accurate",
             "Kalign new fast", "Kalign new accurate", "Clustal Omega", "MAFFT L-INS-i", "MUSCLE 5"]

tables = {
    "overall_protein": md_to_html("protein: held out, benchmark-weighted (n = 2916)"),
    "overall_nucleotide": md_to_html("nucleotide: held out, benchmark-weighted (n = 3225)"),
    "bali_sp": md_to_html("BAliBASE 3, held out: SP", KEEP_MAIN),
    "bali_tc": md_to_html("BAliBASE 3, held out: TC", KEEP_MAIN),
    "bali_cpu": md_to_html("BAliBASE 3, held out: CPU per case (s)", KEEP_MAIN),
}

# ------------------------------------------------------------------ speed between releases (Addendum A1)
speed = list(csv.DictReader(open(RES / "v1_1" / "figures" / "fig_speed.csv")))
tables["speed"] = table(["pool and mode", "Kalign 3.6.1 CPU per case (s)", "Kalign 3.7.0 CPU per case (s)",
                         "3.6.1 / 3.7.0"],
                        [[f"{r['pool']} {r['mode']}", f"{float(r['cpu_361']):.3f}", f"{float(r['cpu_370']):.3f}",
                          f"{float(r['speedup']):.1f}"] for r in speed])
sp = [float(r["speedup"]) for r in speed]

# ------------------------------------------------------------------ optimism (v2)
opt_rows = []
for pool in ("protein", "nucleotide"):
    for r in tsv(V2 / pool / "optimism.tsv"):
        if r["mode"] in ("fast", "accurate"):
            f = lambda k: f"{float(r[k]):+.3f} [{float(r[k + '_lo']):+.3f}, {float(r[k + '_hi']):+.3f}]"  # noqa: E731
            opt_rows.append([f"{pool} {r['mode']}", r["weighting"], f("opt_sp"), f("opt_tc")])
tables["optimism"] = table(["pool and mode", "weighting", "optimism SP [95% CI]", "optimism TC [95% CI]"], opt_rows)

# ------------------------------------------------------------------ per-fold picks (v2)
fold_rows = []
for pool in ("protein", "nucleotide"):
    for fold in ("0", "1", "2", "3", "4", "final"):
        s = json.loads((V2 / pool / f"selection_{'final' if fold == 'final' else 'fold' + fold}.json").read_text())
        n = CORE_N[pool][fold]
        cells = [f"{pool} {('final fit' if fold == 'final' else 'fold ' + fold)}"]
        for m in ("fast", "accurate"):
            e = s["selection"][m]
            runs = e["params"]["n_runs"]
            cells += [f"{e['q']:.3f}", f"{e['core']['cpu'] / n:.3f}", f"{runs}"]
        fold_rows.append(cells)
tables["folds"] = table(["fit", "fast (SP+TC)/2", "fast CPU (s/case)", "fast runs", "accurate (SP+TC)/2",
                         "accurate CPU (s/case)", "accurate runs"], fold_rows)

# ------------------------------------------------------------------ laptop reproduction
rep = []
for line in (V2 / "repro_laptop" / "compare_workstation.txt").read_text().splitlines()[1:6]:
    p = line.split()
    rep.append([" ".join(p[:-8]), p[-8], p[-7], p[-6], p[-5], p[-4], p[-1]])
tables["repro"] = table(["preset", "cases", "SP laptop", "SP workstation", "TC laptop", "TC workstation",
                         "cases within 0.001 SP"], rep)

# ------------------------------------------------------------------ parameters
P = json.loads((V2 / "presets_v2.json").read_text())["presets"]
MAT = {1: "PFASUM43", 2: "PFASUM60", 3: "CorBLOSUM66", 10: "Kimura 1 PAM", 11: "Kimura 20 PAM", 12: "Kimura 200 PAM"}
REF = {0: "none", 1: "all", 2: "confident", 3: "inline"}
FILL = {0: "off", 1: "realign", 2: "project", 3: "attach"}


def params_html(pool, mode):
    p = P[pool][mode]["params"]
    rows = [[str(k + 1), MAT[p["run_types"][k]], f"{p['run_gpo'][k]:.3f}", f"{p['run_gpe'][k]:.3f}",
             f"{p['run_tgpe'][k]:.3f}", f"{p['run_vsm_amax'][k]:.3f}", REF[p["run_refine"][k]]]
            for k in range(p["n_runs"])]
    t = table(["run", "matrix", "gap open", "gap extend", "terminal gap extend", "VSM amax", "refinement"], rows)
    shared = (f"sequence weights {p['seq_weights']:.3f}; consistency anchors {p['consistency']} "
              f"(weight {p['consistency_weight']:.3f}); realign {p['realign']}; minimum support {p['min_support']}; "
              f"fill-in {FILL[p['fill_in']]}" + (f" (anchor fraction {p['fill_anchor_frac']:.3f})" if p["fill_in"] else ""))
    return t + f"<p>Shared settings: {shared}. Parameter key <code>{P[pool][mode]['key'][:12]}</code>.</p>"


for pool in ("protein", "nucleotide"):
    for m in ("fast", "accurate"):
        tables[f"params_{pool}_{m}"] = params_html(pool, m)

vals = {"SPEEDUP_RANGE": f"{min(sp):.1f} to {max(sp):.1f} times"}

# ------------------------------------------------------------------ assemble
FIGS = {"fig_tradeoff_v2_protein": V2 / "addendum_a2" / "figures", "fig_tradeoff_v2_nucleotide": V2 / "addendum_a2" / "figures",
        "fig_bali_categories": V2 / "addendum_a2" / "figures", "fig_mode_choice_protein": RES / "v1_1" / "figures",
        "fig_mode_choice_nucleotide": RES / "v1_1" / "figures", "fig_speed": RES / "v1_1" / "figures",
        "fig_convergence": RES / "v1_1" / "figures"}
body = (HERE / "body.html").read_text()
body = re.sub(r"\{SVG:(\w+)\}", lambda m: (FIGS[m.group(1)] / f"{m.group(1)}.svg").read_text(), body)
body = body.replace("{SCHEMATIC}", (HERE / "schematic.svg").read_text())
body = re.sub(r"\{TABLE:(\w+)\}", lambda m: tables[m.group(1)], body)
for k, v in vals.items():
    body = body.replace("{" + k + "}", v)
body = body.replace("{PLACEHOLDERS}", str(body.count("[SOURCE NEEDED")))
left = re.findall(r"\{[A-Z_]+:?\w*\}", body)
if left:
    sys.exit(f"unfilled: {left}")
tpl = (SKILL / "template.html").read_text()
head = tpl[: tpl.index("<article")].replace("/* paste assets/report.css here, unchanged */",
                                             (SKILL / "report.css").read_text())
head = head.replace("<title>TITLE</title>", "<title>Kalign Preset Re-fit v2</title>")
out = head + '<article class="doc">\n' + body + "\n</article>\n</body>\n</html>\n"
(V2 / "report.html").write_text(out)
print("written", V2 / "report.html", len(out))
