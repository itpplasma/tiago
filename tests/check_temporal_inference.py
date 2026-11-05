"""Independent dense Gaussian-conditioning oracle for native temporal inference.

Only synthetic observations are generated. The Fortran fixture assembles Markov
innovation residuals; this oracle instead conditions a dense stationary prior
in observation space. All numbers are dimensionless.
"""
import argparse
import csv
import math
from pathlib import Path
import random
import subprocess


def transpose(a):
    return list(map(list, zip(*a)))


def multiply(a, b):
    columns = transpose(b)
    return [[sum(x*y for x, y in zip(row, col)) for col in columns] for row in a]


def cholesky(a):
    n = len(a)
    lower = [[0.0]*n for _ in a]
    for i in range(n):
        for j in range(i+1):
            value = a[i][j]-sum(lower[i][k]*lower[j][k] for k in range(j))
            lower[i][j] = math.sqrt(value) if i == j else value/lower[j][j]
    return lower


def solve(lower, rhs):
    n = len(rhs)
    out = [0.0]*n
    for i in range(n):
        out[i] = (rhs[i]-sum(lower[i][j]*out[j] for j in range(i)))/lower[i][i]
    for i in reversed(range(n)):
        out[i] = (out[i]-sum(lower[j][i]*out[j] for j in range(i+1, n)))/lower[i][i]
    return out


def prior(nt, sd, rho):
    n = 4*nt+3
    cov = [[0.0]*n for _ in range(n)]
    for offset, scale, correlation in [(0, .28, .88)]+[
            (nt+3+s*nt, sd, rho) for s in range(3)]:
        for i in range(nt):
            for j in range(nt):
                cov[offset+i][offset+j] = scale**2*correlation**abs(i-j)
    for i in range(nt, nt+3):
        cov[i][i] = .35**2
    return cov


def condition(cov, design, y, noise):
    cross = multiply(cov, transpose(design))
    marginal = multiply(design, cross)
    for i, sigma in enumerate(noise):
        marginal[i][i] += sigma**2
    lower = cholesky(marginal)
    alpha = solve(lower, y)
    mean = [sum(a*b for a, b in zip(row, alpha)) for row in cross]
    solved = [solve(lower, row) for row in cross]
    reduction = multiply(cross, transpose(solved))
    posterior = [[value-reduction[i][j] for j, value in enumerate(row)]
                 for i, row in enumerate(cov)]
    logz = -.5*(len(y)*math.log(2*math.pi)
                  + 2*sum(math.log(lower[i][i]) for i in range(len(y)))
                  + sum(a*b for a, b in zip(y, alpha)))
    return mean, posterior, logz


def synthetic(nt, stress):
    rng = random.Random(817)

    def ar(sd, rho):
        state = [rng.gauss(0, sd)]
        for _ in range(1, nt):
            state.append(rho*state[-1]+rng.gauss(0, sd*math.sqrt(1-rho**2)))
        return state

    plasma = ar(.28, .88)
    offsets = [rng.gauss(0, .35) for _ in range(3)]
    drift = [ar(.60 if stress else .12, .99 if stress else .94) for _ in range(3)]
    truth = plasma+offsets+sum(drift, [])
    design, y, noise = [], [], []
    # Drop 28 of 120 sensor rows for the default forty-time case.
    for s in range(3):
        for t in range(nt):
            row_number = s*nt+t
            if row_number % 4 == 0 and row_number >= 8:
                continue
            row = [0.0]*len(truth)
            for column in (t, nt+s, nt+3+s*nt+t):
                row[column] = 1.0
            design.append(row)
            y.append(sum(a*b for a, b in zip(row, truth))+rng.gauss(0, .045))
            noise.append(.045)
    suspect_count = len(y)
    for t in sorted({round(i*(nt-1)/4) for i in range(5)}):
        row = [0.0]*len(truth)
        row[t] = 1.0
        design.append(row)
        y.append(plasma[t]+rng.gauss(0, .035))
        noise.append(.035)
    return truth, design, y, noise, suspect_count


def check(executable, directory, nt, stress, anchors):
    directory.mkdir(parents=True, exist_ok=True)
    truth, design, y, noise, suspect_count = synthetic(nt, stress)
    if not anchors:
        design, y, noise = design[:suspect_count], y[:suspect_count], noise[:suspect_count]
    with (directory/'input.txt').open('w') as out:
        out.write(f'{nt} {len(y)}\n')
        for observed, sigma, row in zip(y, noise, design):
            out.write(' '.join(map(str, [observed, sigma]+row))+'\n')
    subprocess.run([executable, str(directory)], check=True)
    with (directory/'modes.csv').open() as source:
        modes = list(csv.DictReader(source))
    assert len(modes) == 4
    assert {mode["family"] for mode in modes} == {f"covariance-{i}" for i in range(1, 5)}
    oracles = [condition(prior(nt, sd, rho), design, y, noise)
               for sd, rho in ((.12, .30), (.12, .94), (.30, .30), (.30, .94))]
    peak = max(result[2] for result in oracles)
    weights = [math.exp(result[2]-peak) for result in oracles]
    total = sum(weights)
    weights = [value/total for value in weights]
    errors = [0.0]*4
    for mode in modes:
        family = int(mode['family'].split('-')[-1])-1
        mean, cov, logz = oracles[family]
        prefix = directory/('mode_'+mode['mode'])
        with Path(str(prefix)+'_parameters.csv').open() as source:
            parameters = list(csv.DictReader(source))
        with Path(str(prefix)+'_covariance.csv').open() as source:
            native_cov = [list(map(float, row)) for row in csv.reader(source)]
        assert len(parameters) == len(mean)
        assert len(native_cov) == len(mean)
        assert all(len(row) == len(mean) for row in native_cov)
        differences = [max(abs(float(row['map'])-value)
                           for row, value in zip(parameters, mean)),
                       max(abs(a-b) for r, s in zip(native_cov, cov) for a, b in zip(r, s)),
                       abs(float(mode['log_evidence'])-logz),
                       abs(float(mode['weight'])-weights[family])]
        errors = [max(a, b) for a, b in zip(errors, differences)]
        assert float(mode['psis_ess']) > 19.99
        assert int(mode['psis_failures']) == 0
    for error, tolerance in zip(errors, (2e-7, 2e-10, 2e-8, 2e-8)):
        assert error < tolerance, (directory.name, errors)
    mean = [sum(w*result[0][i] for w, result in zip(weights, oracles))
            for i in range(len(truth))]
    sd = [math.sqrt(sum(w*(result[1][i][i]+(result[0][i]-mean[i])**2)
                        for w, result in zip(weights, oracles))) for i in range(nt)]
    rmse = math.sqrt(sum((a-b)**2 for a, b in zip(mean[:nt], truth[:nt]))/nt)
    coverage = sum(abs(a-b) <= 1.96*s for a, b, s in zip(mean[:nt], truth[:nt], sd))/nt
    print(f'{directory.name}: observations={len(y)}, suspect={suspect_count}, '
          f'max errors(mean,cov,logZ,weight)={errors}; '
          f'weights={weights}; plasma RMSE={rmse:.6g}, coverage={coverage:.3f}')
    return sd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('directory', type=Path)
    parser.add_argument('--times', type=int, default=40)
    args = parser.parse_args()
    anchored = check(args.executable, args.directory/'anchored', args.times, False, True)
    unanchored = check(args.executable, args.directory/'unanchored', args.times, False, False)
    assert sum(anchored) < sum(unanchored), 'anchors must reduce marginal uncertainty'
    check(args.executable, args.directory/'outside-support', args.times, True, True)


if __name__ == '__main__':
    main()
