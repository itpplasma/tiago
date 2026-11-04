"""Prepare the pinned public VMEC++ CTH-like fixture for native TIAGO tests."""
import hashlib
import json
import math
import os
from pathlib import Path
import sys

source, work = map(Path, sys.argv[1:3])
work.mkdir(parents=True, exist_ok=True)
hashes = {
    'mgrid_cth_like.nc': '2ccb8c566a9b3685681940ba94a0b20e9994073a936e41976d654d69090344ae',
    'cth_like_free_bdy.json': '97c3f78f581bdc86aec9a2bbb347ba0a27b07a5798b5d4488fbf181dc024c8ad',
    'coils.cth_like': 'bb9f2dc51d54d2c2388948f4346cd5c25574e79e083570b3751476e13598c2be',
}
for filename, expected in hashes.items():
    with (source / filename).open('rb') as stream:
        digest = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    if digest.hexdigest() != expected:
        raise ValueError(f'{filename} differs from the pinned public VMEC++ fixture')
indata = json.loads((source / 'cth_like_free_bdy.json').read_text())
indata['mgrid_file'] = os.path.relpath(source / 'mgrid_cth_like.nc', work)
(work / 'input.json').write_text(json.dumps(indata))
# VMEC++ also accepts isolated circular coils; TIAGO consumes filaments.
# Keep the published turns, expanding each circle to a closed polygon.
lines = (source / 'coils.cth_like').read_text().splitlines()
points, coils = [], []
for line in lines[lines.index('mirror NUL') + 1:]:
    if line.strip() == 'end':
        break
    fields = line.split()
    point = [float(v) for v in fields[:4]]
    points.append(point)
    if len(fields) > 4:
        group = int(fields[4])
        if len(points) == 1:
            radius, _, height, turns = points[0]
            points = [[radius * math.cos(t), radius * math.sin(t), height, turns]
                      for t in [2 * math.pi * k / 128 for k in range(129)]]
            points[-1][3] = 0
        coils.append((points, group))
        points = []
with (work / 'coils.coils').open('w') as stream:
    stream.write('periods 5\nbegin filament\nmirror NIL\n')
    for points, group in coils:
        for i, point in enumerate(points):
            stream.write(' '.join(f'{v:.16e}' for v in point))
            if i == len(points) - 1:
                stream.write(f' {group} public-cth')
            stream.write('\n')
    stream.write('end\n')
