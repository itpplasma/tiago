#!/usr/bin/env python3
"""Repeated synthetic recorded-signal fits; no experimental source is accessed."""
import argparse
import concurrent.futures
import csv
import os
from pathlib import Path
import subprocess


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True,
                        help='directory containing the noiseless fit.csv')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--realisations', type=int, default=64)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--workers', type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    sensors = Path(__file__).with_name('sensors.awk')
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), MPLBACKEND='Agg')
    workers = args.workers or max(1, (os.cpu_count() or 1)//args.threads)
    baseline = (args.baseline/'fit.csv').resolve()
    if not baseline.exists():
        raise FileNotFoundError(baseline)

    def run(index):
        root = (args.output/f'realisation-{index:03d}').resolve()
        root.mkdir(exist_ok=True)
        subprocess.run(['awk', '-v', 'R0=1.42', '-v', 'rs=.85', '-v', 'nphi=4',
                        '-v', 'nth=4', '-v', 'nprobe=12', '-v', 'dir='+str(root),
                        '-f', str(sensors)], check=True, env=env)
        with (root/'fit.log').open('w') as log:
            subprocess.run([str(args.binary.resolve()), str(args.input.resolve()),
                            str(root), str(root/'work'), str(root),
                            str(20261028+index), str(baseline)],
                           check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
        return index, read(root/'work/profiles.csv'), read(root/'work/calibrations.csv')

    rows, calibrations = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run, i) for i in range(args.realisations)]
        for future in concurrent.futures.as_completed(futures):
            index, result, calibration = future.result()
            for row in result:
                truth, mean, sd = (float(row[k]) for k in ('truth', 'map', 'sd'))
                row.update(realisation=index, covered=int(abs(mean-truth)<=1.95996398454*sd))
                rows.append(row)
            for row in calibration:
                truth, mean, sd = (float(row[k]) for k in ('truth', 'map', 'sd'))
                row.update(realisation=index, covered=int(abs(mean-truth)<=1.95996398454*sd))
                calibrations.append(row)
            print(f'completed realisation {index+1}/{args.realisations}', flush=True)
    with (args.output/'coverage.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r:(int(r['realisation']),int(r['profile']),float(r['s']))))
    with (args.output/'calibration-coverage.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(calibrations[0]))
        writer.writeheader()
        writer.writerows(sorted(calibrations,key=lambda r:(int(r['realisation']),r['parameter'])))
    for profile in ('1','2'):
        selected = [r for r in rows if r['profile']==profile]
        rates = []
        for s in sorted({r['s'] for r in selected}):
            values = [int(r['covered']) for r in selected if r['s']==s]
            rates.append(sum(values)/len(values))
        print(f'profile {profile}: pointwise 95% coverage {min(rates):.3f}..{max(rates):.3f}')


if __name__ == '__main__':
    main()
