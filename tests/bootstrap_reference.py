"""tiago_bootstrap against TORAX's Sauter bootstrap current on the FreeGS
TestTokamakSensor equilibrium (tests/freegs/generate_bootstrap_reference.py).

Usage: bootstrap_reference.py TIAGO_BOOTSTRAP GEQDSK KINETIC REFERENCE TIAGO_AXISYM_FIT
                              MACHINE WORKDIR

Trapped fraction, electron collisionality and <j.B>_bs at nine labels, with
collisionalities from 0.03 to 100. Measured differences are below 5e-4, 4e-3
and 1.5e-3, mostly from the interpolation of q, so the tolerances are 2e-3,
1e-2 and 5e-3. Also checks that the bootstrap current flows along the plasma
current (same sign as the equilibrium's <j.B>), and the equilibrium's
<j.B> = F p' + F' <B^2> / mu0 of tiago_bootstrap and of the jdotb rows of
tiago_axisym_fit against the reference's contour averages, within 1e-3 of the
largest |<j.B>|, since F p' and F' <B^2> / mu0 cancel near the edge.
"""
import csv
import subprocess
import sys
from pathlib import Path

exe, geqdsk, kinetic, reference, fit, machine, work = sys.argv[1:8]
ref = list(csv.DictReader(open(reference)))
run = subprocess.run([exe, geqdsk, kinetic, "0.3"] + [r["psin"] for r in ref],
                     check=True, capture_output=True, text=True)
got = list(csv.DictReader(run.stderr.splitlines()))
rows = list(csv.DictReader(run.stdout.splitlines()))
failures = 0
for g, r, row in zip(got, ref, rows):
    for key, rkey, tol in (("f_trap", "f_trap", 2e-3), ("nu_e", "nu_e", 1e-2), ("jdotb_bs", "jdotb", 5e-3)):
        d = float(g[key]) / float(r[rkey]) - 1
        if abs(d) > tol:
            print(f"psi_N {float(r['psin']):.2f} {key}: relative difference {d:.2e}")
            failures += 1
    if float(g["jdotb_bs"]) * float(g["jdotb_eq"]) <= 0:
        print(f"psi_N {float(r['psin']):.2f}: bootstrap current against the plasma current")
        failures += 1
    if abs(float(row["sigma"]) - 0.3 * abs(float(row["value"]))) > 1e-9 * abs(float(row["value"])):
        print("sigma row does not follow relative_sigma")
        failures += 1
scale = max(abs(float(r["jdotb_eq"])) for r in ref)
worst = max(abs(float(g["jdotb_eq"]) - float(r["jdotb_eq"])) for g, r in zip(got, ref)) / scale
if worst > 1e-3:
    print(f"<j.B>_eq of tiago_bootstrap: largest difference {worst:.2e} of the largest |<j.B>|")
    failures += 1

# the same labels as jdotb rows of the free-boundary fit, predicted from the G-EQDSK
work = Path(work)
work.mkdir(parents=True, exist_ok=True)
labels = [line.split()[1] for line in open(machine)
          if line.split() and line.split()[0] in ("flux", "pickup", "rogowski")]
coils = [line for line in open(machine) if line.split() and line.split()[0] == "coil"]
groups = {line.split()[2] for line in coils}
with open(work / "measurements.csv", "w") as fh:
    fh.write("label,value,sigma\n")
    for label in labels + [f"jdotb@{float(r['psin'])}" for r in ref]:
        fh.write(f"{label},0,1\n")
(work / "predict.nml").write_text(f"""&axisym_fit
 machine = '{machine}', provider = 'none', workdir = '{work / "solve"}',
 measurements = '{work / "measurements.csv"}', output = '{work / "predict"}',
 coil_current = {", ".join("0" for _ in groups)}, predict_geqdsk = '{geqdsk}',
/
""")
subprocess.run([fit, str(work / "predict.nml")], check=True, capture_output=True)
rows = [r for r in csv.DictReader(open(work / "predict_signals.csv")) if r["label"].startswith("jdotb@")]
worst = max(abs(float(p["predicted"]) - float(r["jdotb_eq"])) for p, r in zip(rows, ref)) / scale
if len(rows) != len(ref) or worst > 1e-3:
    print(f"jdotb rows of tiago_axisym_fit: {len(rows)} rows, largest difference {worst:.2e}")
    failures += 1
print(f"{len(got)} labels, {failures} failures")
sys.exit(1 if failures or len(got) != len(ref) else 0)
