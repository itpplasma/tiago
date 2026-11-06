"""Run TIAGO's conjugate joint-signal model and a small withheld-sensor study."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.stats import t


def predict(executable, samples, queries, *, mean0=None, psi0=None,
            kappa0=.05, nu0=None):
    """Samples are complete rows; query indices are zero-based observed columns."""
    samples = np.asarray(samples, dtype=float)
    n, d = samples.shape
    mean0 = np.zeros(d) if mean0 is None else np.asarray(mean0)
    nu0 = d+3 if nu0 is None else nu0
    psi0 = (nu0-d-1)*np.eye(d) if psi0 is None else np.asarray(psi0)
    with tempfile.TemporaryDirectory() as directory:
        source, target = Path(directory)/'input.txt', Path(directory)/'output.txt'
        with source.open('w') as stream:
            stream.write(f'{d} {n} {len(queries)} {kappa0:.17g} {nu0:.17g}\n')
            np.savetxt(stream, mean0[None], fmt='%.17g')
            np.savetxt(stream, psi0, fmt='%.17g')
            np.savetxt(stream, samples, fmt='%.17g')
            for indices, values in queries:
                stream.write(f'{len(indices)}\n')
                if len(indices):
                    np.savetxt(stream, (np.asarray(indices)+1)[None], fmt='%d')
                    np.savetxt(stream, np.asarray(values)[None], fmt='%.17g')
        subprocess.run([str(executable), str(source), str(target)], check=True)
        result = []
        with target.open() as stream:
            for index in range(len(queries)):
                number, count, df, log_observed = stream.readline().split()
                if int(number) != index+1:
                    raise ValueError('unexpected prediction index')
                count = int(count)
                missing = (np.fromstring(stream.readline(), sep=' ', dtype=int)-1
                           if count else np.empty(0, int))
                location = (np.fromstring(stream.readline(), sep=' ')
                            if count else np.empty(0))
                scale = np.array([np.fromstring(stream.readline(), sep=' ')
                                  for _ in range(count)]).reshape(count, count)
                result.append(dict(missing=missing, location=location, scale=scale,
                                   df=float(df), log_observed=float(log_observed)))
    return result


def scalar_summary(predictions, truth):
    centre = np.array([p['location'][0] for p in predictions])
    half = np.array([t.ppf(.975, p['df'])*np.sqrt(p['scale'][0, 0])
                     for p in predictions])
    return dict(rmse=float(np.sqrt(np.mean((truth-centre)**2))),
                coverage_95=float(np.mean(abs(truth-centre) <= half)),
                mean_interval_width=float(np.mean(2*half))), centre, half


def study(executable, output):
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(4182)
    def sensors(n, noise=.18):
        latent = rng.normal(size=(n, 2))
        return (latent@np.array([[1., .8, 1.1], [0., .3, -.2]])
                + rng.normal(size=(n, 3))*np.array([.15, .2, noise]))
    training = sensors(80)
    nominal = sensors(600)
    shifted = nominal.copy()
    shifted[:, 2] += .9
    changed = nominal.copy()
    changed[:, 2] = -nominal[:, 2]
    cases = dict(stable=nominal, target_offset=shifted, correlation_reversal=changed)
    outcomes, arrays = {}, {}
    for name, values in cases.items():
        queries = [(np.array([0, 1]), row[:2]) for row in values]
        result = predict(executable, training, queries)
        summary, centre, half = scalar_summary(result, values[:, 2])
        marginal = predict(executable, training, [(np.array([], int), np.array([]))])[0]
        baseline = marginal['location'][2]
        summary['marginal_rmse'] = float(np.sqrt(np.mean((values[:, 2]-baseline)**2)))
        outcomes[name] = summary
        arrays[name+'_truth'] = values[:, 2]
        arrays[name+'_prediction'] = centre
        arrays[name+'_half95'] = half
    # The classes differ only in the missing sensor; their observed likelihoods agree.
    shifted_training = training.copy()
    shifted_training[:, 2] += 2.
    shifted_prior = np.array([0., 0., 2.])
    queries = [(np.array([0, 1]), row[:2]) for row in nominal[:20]]
    first = predict(executable, training, queries)
    second = predict(executable, shifted_training, queries, mean0=shifted_prior)
    difference = max(abs(a['log_observed']-b['log_observed'])
                     for a, b in zip(first, second))
    if difference > 1e-12:
        raise AssertionError('unobserved class difference changed observed likelihood')
    # Two measurements x+noise with sigma=.5 and x~N(0,1).
    # Their average has noise variance .125 but is already determined by them.
    correct_variance = 1/(1+2/.5**2)
    naive_variance = 1/(1+2/.5**2+1/.125)
    outcomes['information_checks'] = dict(
        missing_only_class_log_likelihood_difference=difference,
        equal_class_prior_posterior=.5,
        correct_latent_variance=correct_variance,
        variance_after_double_counting_prediction=naive_variance)
    outcomes['provenance'] = dict(
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        training_rows=80, test_rows=600, seed=4182,
        interval='Pointwise 95% Student-t posterior predictive for a recorded signal',
        independence='Independent synthetic rows; target-only changes are deliberately unmodelled')
    (output/'results.json').write_text(json.dumps(outcomes, indent=2)+'\n')
    np.savez_compressed(output/'predictions.npz', **arrays)
    print(json.dumps(outcomes, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    study(args.executable.resolve(), args.output)
