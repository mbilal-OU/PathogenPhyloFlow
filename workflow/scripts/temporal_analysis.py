import csv
import json
import subprocess
from pathlib import Path

from Bio import Phylo
from pathogenphyloflow.metrics import decimal_year, linear_regression


settings = dict(snakemake.params.settings)
mode = settings["mode"]
out_screen = Path(snakemake.output.screen)
out_points = Path(snakemake.output.points)
out_status = Path(snakemake.output.status)
log_path = Path(snakemake.log[0])
for path in [out_screen.parent, log_path.parent]:
    path.mkdir(parents=True, exist_ok=True)

with open(snakemake.input.samples, newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
dates = {}
for row in rows:
    value = row.get("date", "").strip()
    if value:
        dates[row["sample"].strip()] = value

tree = Phylo.read(snakemake.input.tree, "newick")
try:
    tree.root_at_midpoint()
except Exception:
    pass

points = []
for tip in tree.get_terminals():
    if tip.name not in dates:
        continue
    try:
        x = decimal_year(dates[tip.name])
    except Exception:
        continue
    distance = tree.distance(tree.root, tip)
    points.append((tip.name, dates[tip.name], x, distance))

with open(out_points, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample", "date", "decimal_year", "root_to_tip_distance"])
    writer.writerows(points)

if len(points) >= 2:
    x = [p[2] for p in points]
    y = [p[3] for p in points]
    regression = linear_regression(x, y)
    span = max(x) - min(x)
else:
    regression = {"slope": 0.0, "intercept": 0.0, "r2": 0.0}
    span = 0.0

criteria = {
    "enough_dated_samples": len(points) >= int(settings["min_dated_samples"]),
    "enough_sampling_span": span >= float(settings["min_year_span"]),
    "positive_slope": regression["slope"] > 0,
    "r2_threshold_met": regression["r2"] >= float(settings["min_r2"]),
}
passed = all(criteria.values())
run_treetime = mode == "on" or (mode == "auto" and passed)

status = "screen_only"
treetime_returncode = None
if mode == "off":
    status = "disabled"
    run_treetime = False
elif run_treetime:
    dates_csv = out_screen.parent / "treetime_dates.csv"
    with open(dates_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "date"])
        for sample, raw_date, _, _ in points:
            writer.writerow([sample, raw_date])
    treetime_dir = out_screen.parent / "treetime"
    treetime_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "treetime",
        "--tree",
        str(snakemake.input.tree),
        "--aln",
        str(snakemake.input.alignment),
        "--dates",
        str(dates_csv),
        "--outdir",
        str(treetime_dir),
        "--reroot",
        "least-squares",
        "--coalescent",
        str(settings["coalescent"]),
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True)
    log_path.write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
    treetime_returncode = completed.returncode
    if completed.returncode != 0:
        status = "treetime_failed"
        raise RuntimeError(f"TreeTime failed; see {log_path}")
    status = "treetime_completed"
elif mode == "auto" and not passed:
    status = "auto_skipped_screen_not_met"
else:
    log_path.write_text("TreeTime was not requested for this run.\n", encoding="utf-8")

summary = {
    "mode": mode,
    "dated_samples": len(points),
    "sampling_span_years": span,
    "slope": regression["slope"],
    "intercept": regression["intercept"],
    "r2": regression["r2"],
    "criteria": criteria,
    "screen_passed": passed,
    "treetime_status": status,
    "treetime_returncode": treetime_returncode,
    "warning": "Root-to-tip regression is a diagnostic and is not proof of temporal signal.",
}
out_screen.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
out_status.write_text(status + "\n", encoding="utf-8")
