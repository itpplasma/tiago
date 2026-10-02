#!/usr/bin/env python3
"""Repeat exactly two portable baseline/series CLI pairs on M16N08 inputs."""
import argparse
import importlib.util
import json
from pathlib import Path
import shlex
import shutil
from statistics import median
from types import SimpleNamespace

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--baseline-checkout", type=Path, required=True)
parser.add_argument("--candidate-checkout", type=Path, required=True)
parser.add_argument("--data-root", type=Path, required=True)
parser.add_argument("--work", type=Path, required=True)
parser.add_argument("--cpu-single", required=True)
parser.add_argument("--cpu-four", required=True)
args = parser.parse_args()
root = args.data_root.resolve()
work = args.work.resolve()
work.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("bench", root / "benchmarks/xdiagno/bench.py")
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)
bench.WORK = work
bench.ENV.pop("OMP_PROC_BIND", None)
bench.ENV.pop("OMP_PLACES", None)
bench.ENV.update(OMP_DYNAMIC="FALSE", OMP_DISPLAY_AFFINITY="TRUE")
runners = {}
for name, checkout in [("baseline", args.baseline_checkout),
                       ("candidate", args.candidate_checkout)]:
    wrapper = work / f"{name}-cli"
    wrapper.write_text(
        "#!/bin/bash\nset -e\n"
        f"cd {shlex.quote(str(checkout.resolve()))}\n"
        f"benchmark_cpus={shlex.quote(args.cpu_single)}\n"
        'if [[ "${OMP_NUM_THREADS:-1}" == 4 ]]; then '
        f"benchmark_cpus={shlex.quote(args.cpu_four)}; fi\n"
        'exec env -u OMP_PROC_BIND -u OMP_PLACES taskset -c "$benchmark_cpus" '
        'fo exec tiago_vacuum_cli "$@" '
        f"2>{shlex.quote(str(work / (name + '.affinity.log')))}\n")
    wrapper.chmod(0o755)
    runners[name] = bench.Runner(SimpleNamespace(
        tiago=str(wrapper), xdiagno=str(root / "benchmarks/xdiagno/_work/bin/xdiagno"),
        xdiagno_patched="", ncpu=4, xvmec=""))
source = root / "benchmarks/xdiagno/_work/runs/geom_m16n08"
target = work / "geom_m16n08"
shutil.copytree(source, target, dirs_exist_ok=True)
coil, flux, seg = [target / p for p in ["coils.in", "flux.diagno", "seg.diagno"]]
rows = []
for threads in [1, 4]:
    _, xf, xs = runners["baseline"].xdiagno_run(target, coil, flux, seg, 16, threads)
    reference = bench.merged(xf, xs)
    row = dict(case="geom_m16n08", samples=16, threads=threads,
               baseline=[], candidate=[], accuracy=[], affinity=[])
    for repeat in range(2):
        order = ["baseline", "candidate"] if repeat == 0 else ["candidate", "baseline"]
        for name in order:
            elapsed, fl, sg = runners[name].tiago_run(target, coil, flux, seg, 16, threads)
            row[name].append(elapsed)
            row["accuracy"].append(dict(implementation=name,
                                        **bench.compare(bench.merged(fl, sg), reference)))
            row["affinity"].append(dict(implementation=name, repeat=repeat,
                                       log=(work / f"{name}.affinity.log").read_text()))
    row["paired_speedups"] = [a / b for a, b in zip(row["baseline"], row["candidate"])]
    row["median"] = {name: median(row[name]) for name in ["baseline", "candidate"]}
    rows.append(row)
    (work / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps({k: row[k] for k in ["threads", "median", "paired_speedups"]}),
          flush=True)
