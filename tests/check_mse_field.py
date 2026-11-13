"""tiago_mse_field on a nested equilibrium given both as VMEC wout and as SPECTRE.

Six channels meet at points of the independent VMEC Fourier/SciPy field table
(qa_vmec_cartesian.dat). The expected angle is computed here from that table's
field with the pi-polarization definition of docs/mse.md; both field sources
must reproduce it within the 2e-3 relative field agreement of the two codes.
"""
import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROWS = [61, 121, 181, 241, 301, 361]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]


def unit(a):
    n = math.sqrt(dot(a, a))
    return [x/n for x in a]


def channel(point, angle):
    v = [math.cos(angle + 2), math.sin(angle + 2), 0.1]
    n = [math.cos(angle), math.sin(angle), 0.3]
    return v, n, [p - 0.7*x for p, x in zip(point, v)], [p + 0.9*x for p, x in zip(point, n)]


def expected(v, n, b):
    n = unit(n)
    z0 = unit(sub([0, 0, 1], [n[2]*x for x in n]))
    h = cross(n, z0)
    e = cross(unit(v), b)
    g = math.atan2(dot(h, e), dot(z0, e))
    return (g + math.pi/2) % math.pi - math.pi/2


def main(executable, wout, spectre, table):
    lines = [ln.split() for ln in Path(table).read_text().splitlines() if not ln.startswith('#')]
    rows = [[float(x) for x in lines[k - 1]] for k in ROWS]
    with tempfile.TemporaryDirectory() as tmp:
        channels = Path(tmp)/'channels.txt'
        text, reference = [str(len(rows))], []
        for k, row in enumerate(rows):
            v, n, beam, sight = channel(row[:3], 0.31*(k + 1))
            text.append(' '.join([f'QA{k}'] + [repr(x) for x in beam + v + sight + n + [0, 0, 1, 1e-9]]))
            reference.append((row[:3], expected(v, n, row[3:])))
        channels.write_text('\n'.join(text) + '\n')
        for option, path in [('--wout', wout), ('--spectre', spectre)]:
            out = Path(tmp)/'out.csv'
            subprocess.run([executable, '--mse', str(channels), option, path, '--output', str(out)],
                           check=True)
            result = list(csv.DictReader(out.open()))
            assert [r['label'] for r in result] == [f'QA{k}' for k in range(len(rows))]
            for r, (point, gamma) in zip(result, reference):
                midpoint = [float(r[c]) for c in ('x_m', 'y_m', 'z_m')]
                assert max(abs(a - b) for a, b in zip(midpoint, point)) < 1e-8, (option, midpoint, point)
                error = float(r['gamma_rad']) - gamma
                error = (error + math.pi/2) % math.pi - math.pi/2
                print(f'{option} {r["label"]} gamma {float(r["gamma_rad"]):+.6f} reference {gamma:+.6f}')
                assert abs(error) < 2e-3, (option, r['label'], error)
    print('check_mse_field passed')


if __name__ == '__main__':
    main(*sys.argv[1:])
