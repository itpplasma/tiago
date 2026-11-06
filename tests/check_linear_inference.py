"""Observation-space Gaussian oracle for the native latent-space inference."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


def main(executable):
    rng = np.random.default_rng(1975)
    cases = []
    for m, n in [(9, 3), (2, 5), (21, 6)]:
        a = rng.normal(size=(m, n))
        sd = np.exp(rng.normal(size=n))
        y = rng.normal(size=m)
        cases.append((a, y, sd))
    a = np.array([[1.], [0.]])
    cases.append((a, np.array([2., 1.e8]), np.ones(1)))
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source, target = root/'input.txt', root/'output.txt'
        with source.open('w') as f:
            f.write(f'{len(cases)}\n')
            for a, y, sd in cases:
                m, n = a.shape
                f.write(f'{m} {n} {-m*np.log(2*np.pi)/2:.17g}\n')
                np.savetxt(f, sd[None, :], fmt='%.17g')
                np.savetxt(f, np.column_stack([y, a]), fmt='%.17g')
        outputs = []
        for threads in [1, 2]:
            env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1')
            subprocess.run([executable, str(source), str(target)], env=env, check=True)
            with target.open() as f:
                results = []
                for k, (a, y, sd) in enumerate(cases):
                    header = f.readline().split()
                    assert int(header[0]) == k+1
                    n = len(sd)
                    mean = np.fromstring(f.readline(), sep=' ')
                    cov = np.array([np.fromstring(f.readline(), sep=' ') for _ in sd])
                    prior = np.diag(sd**2)
                    cross = prior @ a.T
                    observation = np.eye(len(y)) + a @ cross
                    expected_mean = cross @ np.linalg.solve(observation, y)
                    expected_cov = prior - cross @ np.linalg.solve(observation, cross.T)
                    np.testing.assert_allclose(mean, expected_mean, atol=2e-8, rtol=2e-8)
                    np.testing.assert_allclose(cov, expected_cov, atol=2e-12, rtol=2e-12)
                    evidence = -.5*(len(y)*np.log(2*np.pi) + np.linalg.slogdet(observation)[1]
                                     + y @ np.linalg.solve(observation, y))
                    np.testing.assert_allclose(float(header[2]), evidence, atol=2e-10, rtol=2e-15)
                    results.append(np.r_[float(header[2]), mean, cov.ravel()])
                outputs.append(np.concatenate(results))
        # Independent accuracy is checked above; scheduling may change roundoff.
        roundoff = 32*np.finfo(float).eps
        np.testing.assert_allclose(outputs[0], outputs[1], atol=roundoff, rtol=roundoff)
        source.write_text('1\n1 1 0\n-1\n0 1\n')
        result = subprocess.run([executable, str(source), str(target)], capture_output=True)
        assert result.returncode != 0
    print('Gaussian evidence/moments, underdetermined model, constant floor and threading: PASS')


if __name__ == '__main__':
    main(sys.argv[1])
