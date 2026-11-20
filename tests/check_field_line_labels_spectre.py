"""tiago_field_line_labels on the interior field of a nested SPECTRE state.

The SPECTRE state qa_nested.h5 takes its interfaces from the VMEC equilibrium
of the wout file, so its field lines stay on the VMEC flux surfaces up to the
field difference of the two codes. Lines started at the midplane points of the
independent VMEC table (y = z = 0) must therefore all stay inside, and their
mean label in the VMEC chart must equal the label s of their start point.
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

TOLERANCE = 2e-3


def main():
    exe, wout, spectre, table = sys.argv[1:5]
    radii = []
    for line in Path(table).read_text().splitlines():
        if line.startswith('#'):
            continue
        x, y, z = (float(v) for v in line.split()[:3])
        if x > 0 and y == 0 and z == 0:
            radii.append(x)
    assert len(radii) >= 4, radii
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / 'radii.txt').write_text(''.join(f'{r!r}\n' for r in radii))
        (tmp / 'labels.nml').write_text(
            f"&field_line_labels\n spectre = '{spectre}'\n wout = '{wout}'\n"
            f" radii = '{tmp / 'radii.txt'}'\n turns = 20\n steps = 240\n"
            f" output = '{tmp / 'labels.json'}'\n/\n")
        subprocess.run([exe, str(tmp / 'labels.nml')], check=True)
        result = json.loads((tmp / 'labels.json').read_text())
    case = result['cases'][0]
    for r, s, label, lost in zip(radii, result['s_own'], case['label'], case['lost']):
        print(f'R = {r:.6f} m: s = {s:.6f}, label = {label}, lost = {lost}')
        assert not lost, r
        assert abs(label - s) < TOLERANCE, (r, s, label)


if __name__ == '__main__':
    main()
