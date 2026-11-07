"""Python transport for TIAGO's batch thermal ECE transfer executable."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np

O_MODE, X_MODE = 1, 2


def transfer(executable, rays):
    """Each ray supplies frequency [Hz], mode, harmonic, optional wall
    reflectivity and temperature [eV], an optional weak-damping limit of the
    opacity (default 0.1, see docs/ece-opacity.md), and per-cell te [eV], ne [m^-3],
    bfield [T] and length [m], ordered from the observer into the plasma.
    Returns per ray status, temperature [eV], optical depth and, on success,
    the gradient of temperature with respect to te, ne, bfield and length."""
    with tempfile.TemporaryDirectory() as directory:
        source, target = Path(directory)/'input.txt', Path(directory)/'output.txt'
        with source.open('w') as stream:
            stream.write(f'{len(rays)}\n')
            for ray in rays:
                cells = np.c_[ray['te'], ray['ne'], ray['bfield'], ray['length']]
                stream.write(f"{len(cells)} {ray['frequency']:.17g} {ray.get('mode', X_MODE)} "
                             f"{ray.get('harmonic', 2)} {ray.get('wall_rho', 0.):.17g} "
                             f"{ray.get('wall_te', 0.):.17g} {ray.get('damping_limit', 0.1):.17g}\n")
                np.savetxt(stream, cells, fmt='%.17g')
        subprocess.run([str(executable), str(source), str(target)], check=True)
        results = []
        with target.open() as stream:
            for index, ray in enumerate(rays):
                number, status, temperature, depth = stream.readline().split()
                if int(number) != index+1:
                    raise ValueError('unexpected ECE transfer ray index')
                out = dict(status=int(status), temperature=float(temperature),
                           optical_depth=float(depth))
                if out['status'] == 0:
                    grad = np.array([np.fromstring(stream.readline(), sep=' ')
                                     for _ in range(len(ray['te']))])
                    out.update(te_grad=grad[:, 0], ne_grad=grad[:, 1], bfield_grad=grad[:, 2],
                               length_grad=grad[:, 3])
                results.append(out)
    return results
