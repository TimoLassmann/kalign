"""Assemble report.html: template + report.css + body.html, filling figures, tables and the A1 numbers."""
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
SKILL = Path("~/code/claude-plugins/document-tools/skills/html-report/assets").expanduser()
FIGS = HERE.parent.parent / "figures"  # presets/results/v1_1/figures
A1 = HERE.parent


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t" if str(p).endswith(".tsv") else ","))


speed = rd(FIGS / "fig_speed.csv")
sp = [float(r["speedup"]) for r in speed]
n = {r["pool"]: int(r["n"]) for r in speed}
fmt_x = lambda v: f"{v:.1f}"
vals = {
    "SPEEDUP_RANGE": f"{fmt_x(min(sp))} to {fmt_x(max(sp))} times",
    "N_ALL": f"{n['protein'] + n['nucleotide']:,}",
    "N_PROT": f"{n['protein']:,}",
    "N_NUC": f"{n['nucleotide']:,}",
}
err = sum(int(r["errors_361"]) + int(r["errors_370"]) for r in speed)
vals["SPEED_ERRORS"] = ("No alignment failed in either build." if err == 0 else
                        f"{err} alignments failed and keep their CPU in the means.")

# speed table
rows = "".join(
    f"<tr><th>{r['pool']} {r['mode']}</th><td>{float(r['cpu_361']):.3f}</td><td>{float(r['cpu_370']):.3f}</td>"
    f"<td>{float(r['speedup']):.1f}</td></tr>" for r in speed)
tables = {"speed": "<table><thead><tr><th>pool and mode</th><th>3.6.1 CPU per case (s)</th>"
                   "<th>3.7.0 CPU per case (s)</th><th>3.6.1 / 3.7.0</th></tr></thead><tbody>" + rows + "</tbody></table>"}

# accuracy table from the tradeoff CSVs
acc_rows = []
arm_name = {"3.6.1": "Kalign 3.6.1", "3.7.0": "Kalign 3.7.0", "re-fit": "Kalign re-fit"}
for pool in ("protein", "nucleotide"):
    for r in rd(FIGS / f"fig_tradeoff_{pool}.csv"):
        name = r["mode"] if r["arm"] == "external" else f"{arm_name[r['arm']]} {r['mode']}"
        acc_rows.append(f"<tr><th>{pool}</th><td>{name}</td><td>{float(r['sp']):.3f}</td><td>{float(r['tc']):.3f}</td>"
                        f"<td>{float(r['precision']):.3f}</td><td>{float(r['cpu']):.3f}</td></tr>")
tables["accuracy"] = ("<table><thead><tr><th>pool</th><th>program and mode</th><th>SP</th><th>TC</th><th>precision</th>"
                      "<th>CPU per case (s)</th></tr></thead><tbody>" + "".join(acc_rows) + "</tbody></table>")

extra = {}
if (HERE / "values.txt").exists():
    for line in (HERE / "values.txt").read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            extra[k.strip()] = v.strip()
vals.update(extra)

body = (HERE / "body.html").read_text()
body = re.sub(r"\{SVG:(\w+)\}", lambda m: (FIGS / f"{m.group(1)}.svg").read_text(), body)
body = body.replace("{SCHEMATIC}", (HERE / "schematic.svg").read_text())
body = re.sub(r"\{TABLE:(\w+)\}", lambda m: tables[m.group(1)], body)
body = body.replace("{PLACEHOLDERS}", str(body.count("[SOURCE NEEDED")))
missing = sorted(set(re.findall(r"\{([A-Z0-9_]+)\}", body)) - set(vals))
for k, v in vals.items():
    body = body.replace("{" + k + "}", v)
if missing:
    sys.exit(f"unfilled: {missing}")

tpl = (SKILL / "template.html").read_text()
css = (SKILL / "report.css").read_text()
head = tpl[: tpl.index("<article")]
head = head.replace("/* paste assets/report.css here, unchanged */", css).replace("<title>TITLE</title>",
                                                                                 "<title>Kalign Preset Re-fit</title>")
out = head + '<article class="doc">\n' + body + "\n</article>\n</body>\n</html>\n"
dest = HERE.parent.parent / "report.html"
dest.write_text(out)
print("written", dest, len(out))
