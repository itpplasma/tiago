"""Recovery test of tiago_joint_fit on the public VMEC++ Solov'ev tokamak.

Synthetic ne, Te and Ti point values and one interferometer chord are
predicted from known profiles (predict_only), then fitted from a prior whose
mean is 40 % too high. The MAP profiles must recover the known ones.
Usage: joint_fit.py TIAGO_JOINT_FIT SOLOVEV_JSON WORKDIR
"""
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

exe, solovev, work = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
work.mkdir(parents=True, exist_ok=True)
inp = json.loads(solovev.read_text())
inp["ns_array"], inp["ftol_array"], inp["niter_array"] = [11, 31], [1e-11, 1e-12], [2000, 5000]
(work / "solovev.json").write_text(json.dumps(inp))


def truth(s):
    return 5e19 * (1 - 0.6 * s) + 1e19, 2000 * (1 - s) + 100, 1500 * (1 - s) + 100


def table(path, scale):
    rows = ["s,ne,Te,Ti,Zeff"]
    for k in range(21):
        s = k / 20
        ne, te, ti = truth(s)
        rows.append(f"{s},{ne * scale},{te * scale},{ti * scale},1.5")
    path.write_text("\n".join(rows) + "\n")


def probes(path, values=None):
    rows = ["group,quantity,R,phi,Z,value,sigma"]
    pts = [(1, q, 4.05 + 0.1 * k) for k in range(9) for q in (1, 2)] + [(2, 3, 4.1 + 0.1 * k) for k in range(8)]
    for j, (g, q, r) in enumerate(pts):
        v = values[j] if values else 1.0
        rows.append(f"{g},{q},{r},0.0,0.0,{v},{0.02 * v}")
    path.write_text("\n".join(rows) + "\n")


def chords(path, value=1.0):
    path.write_text("x0,y0,z0,dx,dy,dz,length,value,sigma\n"
                    f"2.5,0.0,0.2,1.0,0.0,0.0,3.0,{value},{0.01 * value}\n")


def namelist(path, out, prior, predict):
    path.write_text(f"""&joint_fit
  vmec_input = '{work / "solovev.json"}'
  output_dir = '{work / out}'
  probes_file = '{work / "probes.csv"}'
  chords_file = '{work / "chords.csv"}'
  prior_table = '{work / prior}'
  names = 'PHIEDGE'
  prior = 1.0
  prior_sigma = 0.01
  breakpoints = 0, 0.25, 0.5, 0.75, 1.0
  z_imp = 6.0
  sigma_f = 1.0, 1.0, 1.0
  l_core = 0.4, l_edge = 0.2, s0 = 0.7, sw = 0.2
  laplace_draws = 50
  predict_only = {".true." if predict else ".false."}
/
""")


def run(nml):
    subprocess.run([exe, str(nml)], check=True, cwd=work)


table(work / "truth.csv", 1.0)
table(work / "prior.csv", 1.4)
probes(work / "probes.csv")
chords(work / "chords.csv")
namelist(work / "predict.nml", "predict", "truth.csv", True)
run(work / "predict.nml")
with open(work / "predict/predictions.csv") as f:
    model = [float(r["model"]) for r in csv.DictReader(f)]
probes(work / "probes.csv", model[:-1])
chords(work / "chords.csv", model[-1])
namelist(work / "fit.nml", "fit", "prior.csv", False)
run(work / "fit.nml")

failures = []
for name in ("fit/profiles.csv", "fit/draws/latents.csv", "fit/draws/profiles.csv"):
    with open(work / name) as f:
        if any(None in r for r in csv.DictReader(f)):
            failures.append((name, "more values than header columns"))
with open(work / "fit/profiles.csv") as f:
    for r in csv.DictReader(f):
        s = float(r["s"])
        if r["quantity"] == "p" or not 0.05 <= s <= 0.9:
            continue
        exact = dict(zip(("ne", "Te", "Ti"), truth(s)))[r["quantity"]]
        sd = 0.5 * (float(r["q84"]) - float(r["q16"]))
        if abs(float(r["map"]) - exact) > max(3 * sd, 0.03 * exact):
            failures.append((r["quantity"], s, float(r["map"]), exact, sd))
summary = dict(line.split(" = ") for line in (work / "fit/summary.txt").read_text().splitlines())
chi2 = float(summary["chi2_map"])
if not math.isfinite(chi2) or chi2 > 5.0:
    failures.append(("chi2_map", chi2))
for f in failures:
    print("FAIL", f)
print(f"chi2_map = {chi2:.3e}, {len(failures)} failures")
sys.exit(1 if failures else 0)
