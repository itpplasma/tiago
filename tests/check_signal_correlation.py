"""Sequential conjugate-update and joint-density conditional oracles."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.stats import multivariate_t


def predict(executable, root, samples, mean0, psi0, kappa0, nu0, queries):
    source, target = root/'input.txt', root/'output.txt'
    n, d = samples.shape
    with source.open('w') as f:
        f.write(f'{d} {n} {len(queries)} {kappa0:.17g} {nu0:.17g}\n')
        np.savetxt(f, mean0[None, :], fmt='%.17g')
        np.savetxt(f, psi0, fmt='%.17g')
        np.savetxt(f, samples, fmt='%.17g')
        for indices, values in queries:
            f.write(f'{len(indices)}\n')
            if len(indices):
                np.savetxt(f, (np.asarray(indices)+1)[None, :], fmt='%d')
                np.savetxt(f, np.asarray(values)[None, :], fmt='%.17g')
    subprocess.run([str(executable), str(source), str(target)], check=True)
    results = []
    with target.open() as f:
        for k in range(len(queries)):
            index, nm, df, log_observed = f.readline().split()
            assert int(index) == k+1
            nm = int(nm)
            missing = np.fromstring(f.readline(), sep=' ', dtype=int)-1 if nm else np.empty(0, int)
            location = np.fromstring(f.readline(), sep=' ') if nm else np.empty(0)
            scale = np.array([np.fromstring(f.readline(), sep=' ') for _ in range(nm)]).reshape(nm, nm)
            results.append((missing, location, scale, float(df), float(log_observed)))
    return results


def sequential_posterior(samples, mean, psi, kappa, nu):
    mean, psi = mean.copy(), psi.copy()
    for sample in samples:
        delta = sample-mean
        psi += kappa/(kappa+1)*np.outer(delta, delta)
        mean += delta/(kappa+1)
        kappa += 1
        nu += 1
    df = nu-len(mean)+1
    return mean, (kappa+1)/(kappa*df)*psi, df


def main(executable):
    rng = np.random.default_rng(612)
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for count, prior_df in [(0, .3), (0, 3.5), (2, 3.5), (41, 3.5)]:
            d = 3
            a = rng.normal(size=(d, d))
            psi = a@a.T+np.eye(d)
            mean = rng.normal(size=d)
            samples = rng.normal(size=(count, d))
            kappa, nu = .3, prior_df+d-1
            queries = [(np.array([], int), np.array([])), (np.array([2, 0]), np.array([.7, -.2])),
                       (np.array([0, 1, 2]), np.array([.2, .4, -.1])),
                       (np.array([1]), np.array([.3]))]
            expected_mean, expected_scale, expected_df = sequential_posterior(samples, mean, psi, kappa, nu)
            results = []
            for threads in [1, 2]:
                os.environ['OMP_NUM_THREADS'] = str(threads)
                results.append(predict(executable, root, samples, mean, psi, kappa, nu, queries))
            for result in results:
                np.testing.assert_allclose(result[0][1], expected_mean, atol=2e-13)
                np.testing.assert_allclose(result[0][2], expected_scale, atol=2e-13)
                np.testing.assert_allclose(result[0][3], expected_df, atol=2e-13)
                assert result[0][4] == 0
                for query, item in zip(queries[1:], result[1:]):
                    observed, values = query
                    marginal = multivariate_t(expected_mean[observed], expected_scale[np.ix_(observed, observed)], expected_df)
                    np.testing.assert_allclose(item[4], marginal.logpdf(values), atol=3e-12)
                    if len(item[0]):
                        def independent_density(x):
                            full = np.zeros(d)
                            full[observed] = values
                            full[item[0]] = x
                            return np.exp(multivariate_t.logpdf(full, expected_mean, expected_scale, expected_df)-marginal.logpdf(values))
                        for x in np.linspace(-2, 3, 11):
                            point = x+np.arange(len(item[0]))*.4
                            actual = multivariate_t.pdf(point, item[1], item[2], item[3])
                            np.testing.assert_allclose(actual, independent_density(point), atol=3e-12)
                        if len(item[0]) == 1:
                            np.testing.assert_allclose(quad(independent_density, -np.inf, np.inf)[0], 1., atol=2e-9)
                    else:
                        assert item[1].size == item[2].size == 0
        # A numerical change of coordinates must transform the prior as well.
        offset, units = np.array([-1., 4., 2.]), np.array([.3, 3., 20.])
        indices, values = np.array([2, 0]), np.array([.7, -.2])
        original = predict(executable, root, samples, mean, psi, kappa, nu,
                           [(indices, values)])[0]
        transformed = predict(executable, root, (samples-offset)/units,
            (mean-offset)/units, psi/np.outer(units, units), kappa, nu,
            [(indices, (values-offset[indices])/units[indices])])[0]
        missing = original[0]
        np.testing.assert_allclose(original[1], offset[missing]+units[missing]*transformed[1], atol=2e-13)
        np.testing.assert_allclose(original[2], transformed[2]*np.outer(units[missing], units[missing]), atol=2e-13)
        np.testing.assert_allclose(original[4], transformed[4]-np.log(units[indices]).sum(), atol=2e-13)
        # Invalid training prior, duplicate indices and nonfinite training rows.
        for bad_nu, bad_query, bad_samples in [
                (1., [(np.array([0]), np.array([0.]))], np.zeros((2, 3))),
                (6., [(np.array([0, 0]), np.array([0., 0.]))], np.zeros((2, 3))),
                (6., [(np.array([0]), np.array([0.]))], np.full((2, 3), np.nan))]:
            try:
                predict(executable, root, bad_samples, np.zeros(3), np.eye(3), 1., bad_nu, bad_query)
            except subprocess.CalledProcessError:
                pass
            else:
                raise AssertionError('invalid model accepted')
    print('Sequential training, marginal/conditional density, normalization, masks and threading: PASS')


if __name__ == '__main__':
    main(sys.argv[1])
