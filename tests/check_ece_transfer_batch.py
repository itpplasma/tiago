"""Transfer-equation oracles for the batch thermal ECE executable."""
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from ece_transfer_batch import transfer  # noqa: E402


def ray(n, te, ne=4e19, b0=2.42, b1=2.58, width=0.2, frequency=140e9, **extra):
    return dict(frequency=frequency, te=np.broadcast_to(te, (n,)).astype(float),
                ne=np.broadcast_to(ne, (n,)).astype(float), bfield=np.linspace(b0, b1, n),
                length=np.full(n, width/n), **extra)


def main(executable):
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    rng = np.random.default_rng(1983)
    # Isothermal medium: T = Te (1 - exp(-tau)) without incoming radiation
    # (Rybicki and Lightman 1979, Eq. 1.30), at thin, grey and thick depths.
    iso = [ray(200, 1500., ne=1e17), ray(200, 1500., ne=3e18), ray(200, 1500., ne=4e19, damping_limit=1e3)]
    for r, out in zip(iso, transfer(executable, iso)):
        assert out['status'] == 0
        expected = r['te'][0]*(-np.expm1(-out['optical_depth']))
        np.testing.assert_allclose(out['temperature'], expected, rtol=1e-12)
    depths = [o['optical_depth'] for o in transfer(executable, iso)]
    assert depths[0] < 0.1 < depths[1] < 3 < depths[2]
    # Splitting every cell into two identical halves leaves the solution unchanged.
    varied = ray(120, 1000. + 2000.*rng.random(120), ne=1e19*(1 + 3*rng.random(120)), damping_limit=1e3)
    halves = {k: (np.repeat(v, 2) if k != 'length' else np.repeat(v/2, 2))
              if isinstance(v, np.ndarray) else v for k, v in varied.items()}
    whole, split = transfer(executable, [varied, halves])
    np.testing.assert_allclose(split['temperature'], whole['temperature'], rtol=1e-11)
    np.testing.assert_allclose(split['optical_depth'], whole['optical_depth'], rtol=1e-11)
    # Reverse-mode gradients against central differences of the executable.
    base = transfer(executable, [varied])[0]
    for key, step, grad in [('te', 1e-3, 'te_grad'), ('ne', 1e13, 'ne_grad'), ('bfield', 1e-7, 'bfield_grad'),
                            ('length', 1e-7, 'length_grad')]:
        cells = rng.choice(120, 6, replace=False)
        rays = []
        for c in cells:
            for sign in (1, -1):
                r = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in varied.items()}
                r[key][c] += sign*step
                rays.append(r)
        out = transfer(executable, rays)
        fd = np.array([(out[2*i]['temperature'] - out[2*i+1]['temperature'])/(2*step) for i in range(len(cells))])
        scale = np.max(np.abs(base[grad]))
        np.testing.assert_allclose(base[grad][cells], fd, rtol=1e-5, atol=1e-6*scale)
    # Malformed and out-of-domain rays fail alone, without stopping the batch.
    bad = ray(10, -5.)
    hot = ray(10, 2e4, ne=1e17)
    good = ray(10, 1500., ne=1e17)
    out = transfer(executable, [bad, good, hot])
    assert out[0]['status'] != 0 and out[2]['status'] != 0 and out[1]['status'] == 0
    # The weak-damping limit is opt-in: the default rejects a dense, hot core
    # cell, a relaxed limit accepts it, and accepted rays do not depend on it.
    core = ray(200, 2000., ne=3e19)
    strict, relaxed = transfer(executable, [core, dict(core, damping_limit=1e3)])
    assert strict['status'] == 2 and relaxed['status'] == 0
    np.testing.assert_allclose(relaxed['temperature'], 2000.*(-np.expm1(-relaxed['optical_depth'])), rtol=1e-12)
    thin = ray(200, 1500., ne=1e17)
    a, b = transfer(executable, [thin, dict(thin, damping_limit=1e3)])
    assert a['status'] == 0 and b['status'] == 0
    assert a['temperature'] == b['temperature'] and np.array_equal(a['te_grad'], b['te_grad'])
    assert transfer(executable, [dict(thin, damping_limit=0.)])[0]['status'] != 0
    print('ECE transfer batch checks passed')


if __name__ == '__main__':
    main(sys.argv[1])
