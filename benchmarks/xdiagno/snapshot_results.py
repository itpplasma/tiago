#!/usr/bin/env python3
"""Save a benchmark run with compiler settings and exact input/binary hashes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--reference-flags", required=True)
    parser.add_argument("--tiago-flags", required=True)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    repo = here.parents[1]
    work = here / "_work"
    args.destination.mkdir(parents=True, exist_ok=True)
    suites = json.loads((work / "results/results.json").read_text())
    names = ["results.json", "results.md"]
    if "equilibrium" in suites:
        names.append("equilibrium_jacobian.csv")
    for name in names:
        source = work / "results" / name
        if source.exists():
            shutil.copy2(source, args.destination / name)
    paths = [repo / "build/tiago_vacuum_cli", here / "bench.py"]
    paths += list((work / "bin").glob("x*"))
    paths += [p for p in (work / "data").iterdir() if p.is_file()]
    cpu = platform.processor()
    if Path("/proc/cpuinfo").exists():
        cpu = next((line.split(":", 1)[1].strip() for line in
                    Path("/proc/cpuinfo").read_text().splitlines()
                    if line.startswith("model name")), cpu)
    provenance = {
        "tiago_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
        "tiago_patch_sha256": hashlib.sha256(subprocess.check_output(
            ["git", "diff", "--binary", "HEAD"], cwd=repo)).hexdigest(),
        "reference_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=work / "STELLOPT",
            text=True).strip(),
        "compiler": subprocess.check_output(
            ["gfortran", "--version"], text=True).splitlines()[0],
        "platform": platform.platform(), "cpu": cpu,
        "threads": args.threads,
        "suites": list(suites),
        "environment": {key: os.environ.get(key) for key in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
        "timing": "one observation per case; whole process including MPI startup "
                  "and input loading; shared workstation",
        "tiago_flags": args.tiago_flags,
        "reference_flags": args.reference_flags,
        "hashes": {str(path.relative_to(repo)): hashlib.sha256(
            path.read_bytes()).hexdigest() for path in sorted(paths)},
    }
    (args.destination / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
