"""Joint mode of tiago_axisym_fit with a provider stand-in that returns the
FreeGS TestTokamakSensor G-EQDSK for every request.

Usage: axisym_joint_cli.py TIAGO_AXISYM_FIT GEQDSK MACHINE STUB WORKDIR

Predict mode writes n_e and T_e at 12 outboard points from known profiles
(the prior table). Fitted with 2 % noise from a prior table 30 % too high,
the MAP profiles must recover the known ones within three Laplace sd or 2 %
at s = 0.05 ... 0.2, where the probes are and the sd is below 5 %, and the
profile tables must be written.
"""
import csv
import subprocess
import sys
from pathlib import Path

import math

exe, geqdsk, machine, stub, work = sys.argv[1:6]
work = Path(work)
work.mkdir(parents=True, exist_ok=True)


def ne(s):
    return 1.0e19 * (1 - 0.6 * s * s)


def te(s):
    return 300.0 * (1 - 0.8 * s * s)


def table(path, factor):
    with open(path, "w") as fh:
        fh.write("s,ne,te,ti,zeff\n")
        for k in range(21):
            s = k / 20
            fh.write(f"{s},{factor * ne(s)},{factor * te(s)},{factor * te(s)},1.5\n")


table(work / "truth.csv", 1.0)
table(work / "prior.csv", 1.3)
labels = [line.split()[1] for line in open(machine)
          if line.split() and line.split()[0] in ("flux", "pickup", "rogowski")]
with open(work / "measurements.csv", "w") as fh:
    fh.write("label,value,sigma\n")
    for label in labels:
        fh.write(f"{label},0,1\n")


def probes(path, values=None):
    with open(path, "w") as fh:
        fh.write("group,quantity,R,phi,Z,value,sigma\n")
        k = 0
        for q in (1, 2):
            for i in range(12):
                v = values[k] if values else 1.0
                fh.write(f"1,{q},{1.15 + 0.30 * i / 11},0,0.05,{v},{0.02 * abs(v)}\n")
                k += 1


def namelist(path, predict, prior):
    path.write_text(f"""&axisym_fit
 machine = '{machine}', provider = 'sh {stub} {geqdsk}', workdir = '{work / "solve"}',
 measurements = '{work / "measurements.csv"}', output = '{work / ("predict" if predict else "fit")}',
 coil_current = 0, 0, 0, 0, plasma_current = 2e5,
 parameters = 'PLASMA_CURRENT', parameter_scale = 2e5, prior_values = 2e5, prior_sigma = 2e3,
 probes_file = '{work / "probes.csv"}', prior_table = '{prior}',
 breakpoints = 0, 0.3, 0.6, 0.85, 0.95, 1.0, z_imp = 6.0, sigma_f = 1.0, 1.0, 1.0,
 model_rel = 0.0, pressure_psin = 0.3, 0.6, pressure_sigma = 1e9, 1e9,
 {"predict_geqdsk = '" + geqdsk + "'," if predict else ""}
/
""")


probes(work / "probes.csv")
namelist(work / "predict.nml", True, work / "truth.csv")
subprocess.run([exe, str(work / "predict.nml")], check=True, capture_output=True)
pred = [float(r["predicted"]) for r in csv.DictReader(open(work / "predict_profile_signals.csv"))]
noise = [math.sin(3.1 * k) for k in range(len(pred))]
probes(work / "probes.csv", [p * (1 + 0.02 * n) for p, n in zip(pred, noise)])
namelist(work / "fit.nml", False, work / "prior.csv")
subprocess.run([exe, str(work / "fit.nml")], check=True, capture_output=True)

failures = 0
for r in csv.DictReader(open(work / "fit_profiles.csv")):
    s = float(r["s"])
    if not 0.05 <= s <= 0.2:
        continue
    for key, f in (("ne", ne), ("te", te)):
        err = abs(float(r[key]) - f(s))
        if float(r[key + "_sd"]) > 0.05 * f(s) or err > max(3 * float(r[key + "_sd"]), 0.02 * f(s)):
            print(f"{key} at s = {s}: {float(r[key]):.4g} against {f(s):.4g}, sd {float(r[key + '_sd']):.3g}")
            failures += 1
names = [r["name"] for r in csv.DictReader(open(work / "fit_parameters.csv"))]
if not any(n.startswith("LATENT(") for n in names):
    print("no profile latents in fit_parameters.csv")
    failures += 1
if len(list(csv.DictReader(open(work / "fit_profile_signals.csv")))) != 24:
    print("fit_profile_signals.csv does not hold the 24 probe signals")
    failures += 1
print(f"{failures} failures")
sys.exit(1 if failures else 0)
