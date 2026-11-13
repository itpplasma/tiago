"""tiago_mse_field on a nested equilibrium given both as VMEC wout and as SPECTRE.

Six channels meet at points of the independent VMEC Fourier/SciPy field table
(qa_vmec_cartesian.dat). The expected angle is computed here from that table's
field with the pi-polarization definition of docs/mse.md; both field sources
must reproduce it within the 2e-3 relative field agreement of the two codes.
"""
import csv
import math
import os
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


def write_channels(path, points):
    text, geometries = [str(len(points))], []
    for k, point in enumerate(points):
        v, n, beam, sight = channel(point, 0.31*(k + 1))
        text.append(' '.join([f'QA{k}'] + [repr(x) for x in beam + v + sight + n + [0, 0, 1, 1e-9]]))
        geometries.append((v, n))
    path.write_text('\n'.join(text) + '\n')
    return geometries


def write_coils(path):
    lines = ['periods 1', 'begin filament', 'mirror NIL']
    for group, (a, current) in enumerate([(0.8, 2), (1.3, -3)], 1):
        for x, y in [(-a, -a), (a, -a), (a, a), (-a, a)]:
            lines.append(f'{x + 4} {y} 0 {current}')
        lines.append(f'{4 - a} {-a} 0 0 {group} square{group}')
    path.write_text('\n'.join(lines + ['end']) + '\n')


def run(executable, channels, out, *options):
    subprocess.run([executable, '--mse', str(channels), *map(str, options), '--output', str(out)],
                   check=True)
    rows = list(csv.DictReader(out.open()))
    fields = [[float(row[key]) for key in ('b_x_T', 'b_y_T', 'b_z_T')] for row in rows]
    return rows, fields


def check(rows, fields, points, geometries):
    assert [row['label'] for row in rows] == [f'QA{k}' for k in range(len(points))]
    for row, field, point, (v, n) in zip(rows, fields, points, geometries):
        midpoint = [float(row[key]) for key in ('x_m', 'y_m', 'z_m')]
        assert max(abs(a - b) for a, b in zip(midpoint, point)) < 1e-8
        error = float(row['gamma_rad']) - expected(v, n, field)
        error = (error + math.pi/2) % math.pi - math.pi/2
        assert abs(error) < 2e-12, (row['label'], error)


def vacuum_cases(executable, grid):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        coils = tmp/'squares.coils'
        write_coils(coils)
        extcur = tmp/'currents.input'
        extcur.write_text('&indata\n  EXTCUR = 1250, -400\n/\n')
        points = [[4, 0, z] for z in (0.2, 0.5, 1.0)]
        channels = tmp/'channels.txt'
        geometries = write_channels(channels, points)
        retained_coils = coils.read_text()
        for options in (['--output', coils],
                        ['--output', tmp/'angles.csv', '--field-jacobian', coils]):
            attempt = subprocess.run([executable, '--mse', str(channels), '--coils', str(coils),
                                      *map(str,options)], capture_output=True, text=True)
            assert attempt.returncode != 0 and 'output aliases input' in attempt.stderr
            assert coils.read_text() == retained_coils
        new_output = tmp/'new.csv'
        attempt = subprocess.run([executable, '--mse', str(channels), '--coils', str(coils),
                                  '--output', str(new_output), '--field-jacobian',
                                  str(tmp)+'/./new.csv'], capture_output=True, text=True)
        assert attempt.returncode != 0 and 'output files must differ' in attempt.stderr
        output = tmp/'protected.csv'
        for alias in (output, tmp/'protected-link.csv', tmp/'protected-hard.csv'):
            output.write_text('retained output\n')
            if alias != output:
                try:
                    if alias.name.endswith('hard.csv'):
                        os.link(output, alias)
                    else:
                        alias.symlink_to(output)
                except OSError:
                    continue
            attempt = subprocess.run([executable, '--mse', str(channels), '--coils', str(coils),
                                      '--output', str(output), '--field-jacobian', str(alias)],
                                     capture_output=True, text=True)
            assert attempt.returncode != 0 and 'output files must differ' in attempt.stderr
            assert output.read_text() == 'retained output\n'
        for bad in ('2,,3', ',2', '2,', '2 /', '2*1', '2 -0.5', 'nan,2',
                    ','.join(['1']*129), ','.join(['1.0']*400)):
            attempt = subprocess.run([executable, '--mse', str(channels), '--mgrid', str(grid),
                                      '--mgrid-currents', bad, '--output', str(tmp/'bad.csv')],
                                     capture_output=True, text=True)
            assert attempt.returncode != 0 and '--mgrid-currents' in attempt.stderr, bad
        for currents, options in [([2, -3], []), ([1250, -400], ['--extcur', extcur])]:
            jacobian = tmp/'jacobian.csv'
            rows, fields = run(executable, channels, tmp/'out.csv', '--coils', coils,
                               '--field-jacobian', jacobian, *options)
            check(rows, fields, points, geometries)
            derivatives = list(csv.DictReader(jacobian.open()))
            for field, (v, n), row in zip(fields, geometries, derivatives):
                strength = math.sqrt(dot(field, field))
                for c, key in enumerate(('dgamma_dBx_rad_per_T', 'dgamma_dBy_rad_per_T',
                                         'dgamma_dBz_rad_per_T')):
                    h = 1e-5
                    plus = [x/strength for x in field]
                    minus = plus[:]
                    plus[c] += h
                    minus[c] -= h
                    delta = expected(v, n, plus) - expected(v, n, minus)
                    delta = (delta + math.pi/2) % math.pi - math.pi/2
                    fd = delta/(2*h)
                    scaled = float(row[key])*strength
                    assert abs(scaled - fd) < 1e-7*max(1, abs(fd)), (key, scaled, fd)
            for point, field in zip(points, fields):
                z = point[2]
                # Closed-form field on the symmetry axis of each square loop.
                bz = sum(4e-7*2*a*a*current/((a*a + z*z)*math.sqrt(2*a*a + z*z))
                         for a, current in zip((0.8, 1.3), currents))
                assert max(abs(field[0]), abs(field[1]), abs(field[2] - bz)) < 2e-13

        nodes = [(7, 5, 2), (12, 10, 4), (16, 14, 7)]
        points, reference = [], []
        for i, j, k in nodes:
            radius, z = 2 + (i - 1)/10, -1 + (j - 1)/9
            phi = 2*math.pi/5*(k - 1)/12
            points.append([radius*math.cos(phi), radius*math.sin(phi), z])
            bcyl = [sum(scale*(group*(1 + 0.1*c) + 0.3*i - 0.2*j + 0.05*k*c - 0.7*group*j*k)
                        for group, scale in [(1, 2), (2, -0.5)]) for c in (1, 2, 3)]
            reference.append([bcyl[0]*math.cos(phi) - bcyl[1]*math.sin(phi),
                              bcyl[0]*math.sin(phi) + bcyl[1]*math.cos(phi), bcyl[2]])
        geometries = write_channels(channels, points)
        rows, fields = run(executable, channels, tmp/'out.csv', '--mgrid', grid,
                           '--mgrid-currents', '2,-0.5')
        check(rows, fields, points, geometries)
        assert max(abs(a - b) for got, want in zip(fields, reference)
                   for a, b in zip(got, want)) < 1e-10
        print('vacuum coil and mgrid MSE checks passed')


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
        coils = Path(tmp)/'squares.coils'
        write_coils(coils)
        _, coil_fields = run(executable, channels, Path(tmp)/'coils.csv', '--coils', coils)
        for option, path in [('--wout', wout), ('--spectre', spectre)]:
            out = Path(tmp)/'out.csv'
            subprocess.run([executable, '--mse', str(channels), option, path, '--output', str(out)],
                           check=True)
            result = list(csv.DictReader(out.open()))
            base_fields = [[float(row[key]) for key in ('b_x_T', 'b_y_T', 'b_z_T')]
                           for row in result]
            assert [r['label'] for r in result] == [f'QA{k}' for k in range(len(rows))]
            for r, (point, gamma) in zip(result, reference):
                midpoint = [float(r[c]) for c in ('x_m', 'y_m', 'z_m')]
                assert max(abs(a - b) for a, b in zip(midpoint, point)) < 1e-8, (option, midpoint, point)
                error = float(r['gamma_rad']) - gamma
                error = (error + math.pi/2) % math.pi - math.pi/2
                print(f'{option} {r["label"]} gamma {float(r["gamma_rad"]):+.6f} reference {gamma:+.6f}')
                assert abs(error) < 2e-3, (option, r['label'], error)
            combined, fields = run(executable, channels, out, option, path, '--coils', coils)
            geometries = [(channel(row[:3], 0.31*(k + 1))[:2]) for k, row in enumerate(rows)]
            check(combined, fields, [row[:3] for row in rows], geometries)
            assert max(abs(total - base - coil) for total_row, base_row, coil_row
                       in zip(fields, base_fields, coil_fields)
                       for total, base, coil in zip(total_row, base_row, coil_row)) < 1e-12
    print('check_mse_field passed')


if __name__ == '__main__':
    if sys.argv[2] == '--vacuum':
        vacuum_cases(sys.argv[1], sys.argv[3])
    else:
        main(*sys.argv[1:])
