"""Recovery of knot profiles with tiago_axisym_fit on FreeGS's TestTokamakSensor.

Usage: axisym_fit.py TIAGO_AXISYM_FIT PYTHON PROVIDER MACHINE WORKDIR

Plants p' and F F' knots, solves with the provider, predicts the sensors and
three pressure rows with `predict_geqdsk`, then fits, starting every solve
from the planted plasma flux, plasma current, coil
currents (measured to 1 %) and four knot values from a perturbed start, and
draws three Laplace samples. Checks recovery within two local standard
deviations and the edge j_phi table.
"""
import csv
import subprocess
import sys
from pathlib import Path

exe, python, provider, machine, work = sys.argv[1:6]
work = Path(work)
work.mkdir(parents=True, exist_ok=True)
currents = [1.53020064498019026e5, 6.27160501916148351e4, -9.80609577221965155e4,
            -5.75578364420235448e4]
knots = [0.0, 0.5, 0.9, 1.0]
# Zero current at psi_N = 1: otherwise the provider's grid nodes entering and
# leaving the plasma make the signals jump (docs/axisymmetric-free-boundary.md).
pprime = [1.0e4, 7.0e3, 3.0e3, 0.0]
ffprime = [0.02, 0.012, 0.006, 0.0]
truth = {"PLASMA_CURRENT": 2.0e5, "PPRIME(2)": pprime[1], "PPRIME(3)": pprime[2],
         "FFPRIME(2)": ffprime[1], "FFPRIME(3)": ffprime[2]}
truth.update({f"CURRENT({g + 1})": c for g, c in enumerate(currents)})
command = f"{python} -I {provider}"


def lst(values):
    return ", ".join(f"{v:.17e}" for v in values)


def namelist(path, predict, start):
    names = list(truth)
    scale = [abs(truth[n]) for n in names]
    prior = [truth[n] if n.startswith("CURRENT") else 0.0 for n in names]
    sigma = [0.01 * abs(truth[n]) if n.startswith("CURRENT") else 0.0 for n in names]
    path.write_text(f"""&axisym_fit
 machine = '{machine}', provider = '{command}', workdir = '{work / "solve"}',
 measurements = '{work / "measurements.csv"}', output = '{work / ("predict" if predict else "fit")}',
 nx = 65, ny = 65, fvac = 2.0, tolerance = 3e-6, coil_current = {lst(currents)}, plasma_current = 2e5,
 profile = 'knots', knot_psin = {lst(knots)}, pprime_knot = {lst(pprime)},
 ffprime_knot = {lst(ffprime)},
 parameters = {", ".join(repr(n) for n in names)},
 start_values = {lst(start[n] for n in names)}, parameter_scale = {lst(scale)},
 prior_values = {lst(prior)}, prior_sigma = {lst(sigma)},
 max_iterations = 20, draws = 3, edge_psin = 0.9, 0.95,
 initial = '{work / "truth.geqdsk.psi.npy"}',
 {"predict_geqdsk = '" + str(work / "truth.geqdsk") + "'," if predict else ""}
/
""")


# Planted equilibrium straight from the provider.
request = work / "truth.txt"
request.write_text("\n".join([
    f"machine = {machine}", "nx = 65", "ny = 65", "rmin = 0.1", "rmax = 2.0",
    "zmin = -1.0", "zmax = 1.0", "tolerance = 1e-7", "plasma_current = 2e5",
    "pressure_axis = 0", "fvac = 2.0", "profile = knots",
    "psin_knots = " + ",".join(map(str, knots)),
    "pprime_knots = " + ",".join(map(str, pprime)),
    "ffprime_knots = " + ",".join(map(str, ffprime))]
    + [f"current({g + 1}) = {c}" for g, c in enumerate(currents)]) + "\n")
subprocess.run(command.split() + [str(request), str(work / "truth.geqdsk")], check=True)

labels = [line.split()[1] for line in open(machine)
          if line.split() and line.split()[0] in ("flux", "pickup", "rogowski")]
rows = [("pressure@0.2", 1.0), ("pressure@0.5", 1.0), ("pressure@0.8", 1.0)]
with open(work / "measurements.csv", "w") as fh:
    fh.write("label,value,sigma\n")
    for label in labels:
        fh.write(f"{label},0,1\n")
    for label, sigma in rows:
        fh.write(f"{label},0,{sigma}\n")
namelist(work / "predict.nml", True, truth)
subprocess.run([exe, str(work / "predict.nml")], check=True)
predicted = list(csv.DictReader(open(work / "predict_signals.csv")))
sigma = {"FL": 1e-4, "BP": 1e-4, "Rog1": 50.0}
with open(work / "measurements.csv", "w") as fh:
    fh.write("label,value,sigma\n")
    for row in predicted:
        label = row["label"]
        s = 0.01 * abs(float(row["predicted"])) if "@" in label else sigma.get(label[:2], sigma.get(label))
        fh.write(f"{label},{row['predicted']},{s}\n")

start = dict(truth)
start["PLASMA_CURRENT"] *= 0.99
for k in ("PPRIME(2)", "PPRIME(3)", "FFPRIME(2)", "FFPRIME(3)"):
    start[k] *= 1.05
for g in range(4):
    start[f"CURRENT({g + 1})"] *= 1.003
namelist(work / "fit.nml", False, start)
subprocess.run([exe, str(work / "fit.nml")], check=True)

failures = 0
for row in csv.DictReader(open(work / "fit_parameters.csv")):
    name, value, sd = row["name"], float(row["value"]), float(row["sd"])
    error = value - truth[name]
    print(f"{name:16s} {value:14.6e} {truth[name]:14.6e} {sd:12.4e}")
    if abs(error) > 2 * sd + 1e-6 * abs(truth[name]):
        print(f"{name} outside two standard deviations")
        failures += 1
edge = list(csv.DictReader(open(work / "fit_edge.csv")))
samples = {int(r["sample"]) for r in edge}
print(f"edge table: {len(edge)} rows from samples {sorted(samples)}")
if 0 not in samples or len(samples) < 3:
    print("edge table misses the MAP or the draws")
    failures += 1
sys.exit(1 if failures else 0)
