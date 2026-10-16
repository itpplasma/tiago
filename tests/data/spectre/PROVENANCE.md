# SPECTRE nested-limit control

Verified source: `plasma/proj/spectre-orbits`, HEAD
`fc5ef9c57ac510ae336e59b3973908debd9c6d6d`,
`data/qa_nested/qa_nested.h5`. Its recorded generator is SPECTRE
`f51638520d2dbdee9b0b4e4174ace5b2349b7460`, eight fixed VMEC interfaces,
`mpol=ntor=6`, iota constraints, `angle_fix=True`. Its per-volume Beltrami
residuals were recorded near `1e-14`; its interfaces were **not** force-relaxed.
This control exercises a nested field representation, not a converged island
equilibrium. The QA VMEC has zero pressure and approximately zero total
toroidal current (`-2.43e-9 A`); it does not validate a finite-current conversion.

The VMEC input is the public simsopt fixture at commit
`9e027eac38028d57aa23777be52a781aa860e347`,
`tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc`.
CTest downloads that exact revision. Hashes:

| File | SHA-256 |
| --- | --- |
| Original SPECTRE HDF5 | `69351709f1ddc3e1409e9d791806cd41068d10a1a15cb7102a4afc462c384038` |
| Committed `qa_nested.h5` | `572fb80ee5997c4d9ece77876c9cc141206f932b16690f0a25feef9b2dac1663` |
| VMEC wout | `333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771` |
| `qa_vmec_cartesian.dat` | `391bb76b8c099f4f3a1566df02f00fff47d67b14fda829aeb7720327b17e7a8b` |

`qa_nested.h5` retains the original dataset values and numeric types for:

- `input/physics/`: `Igeometry`, `Nfp`, `Nvol`, `Mpol`, `Ntor`, `Lrad`,
  `pressure`, `tflux`, `pflux`, `phiedge`;
- `output/`: `Mvol`, `mn`, `im`, `in`, `mu`, `Rbc`, `Zbs`, `Rbs`, `Zbc`;
- `vector_potential/`: `Ate`, `Aze`, `Ato`, `Azo`.

Only unused datasets and attributes were removed; retained arrays were copied
with `h5py.create_dataset(name, data=original[name][()])`, with contiguous
storage and no compression. No equilibrium was recomputed for this packaging.

The independent Cartesian VMEC oracle uses `VmecField` from
`spectre-orbits/python/vmec_spectre_samefield.py`, SHA-256
`256150759268860e9c569b8a173344821165fc05510b2377bce5d6bd8f7a190b`.
It reconstructs VMEC B from the wout contravariant Fourier series and SciPy
cubic radial splines, without using libneo or SPECTRE field evaluation.
To reproduce the 512 rows, with that module on `PYTHONPATH`:

```python
from pathlib import Path
import numpy as np
from vmec_spectre_samefield import VmecField
v = VmecField.from_file(Path("wout_qa.nc"))
rows = []
for iv in range(8):
    for frac in [0.3, 0.7]:
        for j in range(8):
            for k in range(4):
                q = [(iv + frac)/8, 2*np.pi*j/8, 2*np.pi*k/(4*v.nfp)]
                x, b = v.evaluate(np.array(q))
                rows.append(np.r_[x, b])
np.savetxt("qa_vmec_cartesian.dat", rows, fmt="%.17e",
           header="Cartesian x y z [m], Bx By Bz [T]; independent VMEC Fourier + SciPy CubicSpline oracle")
```

`test_spectre_vmec` compares fields at these same physical points. The
acceptance limits are maximum relative B `2e-3`, RMS relative B `6e-4`;
exterior maximum component errors `2e-3 T` and `2e-3 T m`; probe errors
`2e-7 Wb` for area `1e-4 m^2`; linked-loop flux error `5e-2 Wb`; and
area-weighted Rogowski integral error `5e-6 T m^3`. Signed linked flux is
explicitly supplied from VMEC, since libneo does not yet expose it for SPECTRE.
It is metadata for this comparison, not an independently inferred quantity.
The 64- and 128-point casing grids test the same tolerances. Representation
error remains: the order-6 conversion truncates VMEC's order-8 geometry.

The other source test uses libneo's pinned `G3V3L3Fi` fixture and its 90 native
SPECTRE reference rows, whose provenance is recorded in libneo's
`test/spectre/data/PROVENANCE.md`. Its interior maximum relative B tolerance
is `1e-8`; its exterior Cartesian finite-difference curl tolerance is `1e-6`.
