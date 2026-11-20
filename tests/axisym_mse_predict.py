"""MSE rows of tiago_axisym_fit against FreeGS's field on TestTokamakSensor.

Usage: axisym_mse_predict.py TIAGO_AXISYM_FIT GEQDSK MACHINE FREEGS_CSV WORKDIR

One channel at each field point of the FreeGS fixture, at toroidal angle
phi_k = 0.3 k, viewed radially (n = e_R, reference +Z) by a horizontal beam
v = -cos(a) e_R + sin(a) e_phi. For this geometry (docs/mse.md)
tan(gamma) = v_R B_Z / (v_R B_phi - v_phi B_R) with FreeGS's B_R, B_phi, B_Z.
"""
import csv
import math
import subprocess
import sys
from pathlib import Path

exe, geqdsk, machine, reference, work = sys.argv[1:6]
work = Path(work)
work.mkdir(parents=True, exist_ok=True)

field, currents = {}, []
for row in csv.reader(open(reference)):
    kind, label = row[0], row[1]
    if kind in ("br", "bphi", "bz"):
        field.setdefault(label, {"r": float(row[2]), "z": float(row[3])})[kind] = float(row[4])
    elif kind == "current" and label != "plasma":
        currents.append(float(row[4]))

lines, expected = [], {}
for k, (label, f) in enumerate(sorted(field.items())):
    phi, a = 0.3 * k, 0.4 + 0.05 * k
    er = (math.cos(phi), math.sin(phi), 0.0)
    ep = (-math.sin(phi), math.cos(phi), 0.0)
    p = (f["r"] * er[0], f["r"] * er[1], f["z"])
    vr, vp = -math.cos(a), math.sin(a)
    v = tuple(vr * er[i] + vp * ep[i] for i in range(3))
    beam = tuple(p[i] - v[i] for i in range(3))
    sight = tuple(p[i] + 0.5 * er[i] for i in range(3))
    values = [*beam, *v, *sight, *er, 0.0, 0.0, 1.0, 1e-3]
    lines.append(f"MSE_{label} " + " ".join(f"{x:.17e}" for x in values))
    expected[f"MSE_{label}"] = math.atan(vr * f["bz"] / (vr * f["bphi"] - vp * f["br"]))
(work / "mse.txt").write_text(f"{len(lines)}\n" + "\n".join(lines) + "\n")

labels = [line.split()[1] for line in open(machine)
          if line.split() and line.split()[0] in ("flux", "pickup", "rogowski")]
with open(work / "measurements.csv", "w") as fh:
    fh.write("label,value,sigma\n")
    for label in labels + list(expected):
        fh.write(f"{label},0,1\n")
(work / "predict.nml").write_text(f"""&axisym_fit
 machine = '{machine}', provider = 'none', workdir = '{work / "solve"}',
 measurements = '{work / "measurements.csv"}', output = '{work / "predict"}',
 coil_current = {", ".join(f"{c:.17e}" for c in currents)},
 mse_channels = '{work / "mse.txt"}', predict_geqdsk = '{geqdsk}',
/
""")
subprocess.run([exe, str(work / "predict.nml")], check=True)

worst = 0.0
seen = 0
for row in csv.DictReader(open(work / "predict_signals.csv")):
    if row["label"] in expected:
        d = math.remainder(float(row["predicted"]) - expected[row["label"]], math.pi)
        worst = max(worst, abs(d))
        seen += 1
print(f"{seen} MSE channels, largest angle difference {worst:.3e} rad")
sys.exit(0 if seen == len(expected) and worst < 5e-6 else 1)
