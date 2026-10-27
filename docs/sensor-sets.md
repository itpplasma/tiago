# Magnetic sensor sets

A sensor set gathers every magnetic sensor family of a machine in one
manifest, so that a reconstruction reads one file instead of one option per
family. The manifest only names geometry files in the existing DIAGNO
formats and tags each group with a sensor kind. Machine geometry therefore
stays in its own files, and private geometry (for example W7-X) is read at run
time from a directory outside this repository.

## Format

```text
tiago-sensors 1
root ${MACHINE_MAGNETICS}      # optional, default: the manifest's directory
group diamagnetic_loop  dialoops.diagno  select=DIA1,DIA2
group compensation_loop dialoops.diagno  select=COMP1,COMP2
group saddle_coil       saddle.diagno    turns=1
group flux_loop         fluxloops.diagno
group segmented_rogowski rogowski.diagno area=3.4e-4
group bprobe            probes.txt       rphiz
combine DIA1_COMPENSATED 1.0 DIA1 -0.5 COMP1 -0.5 COMP2
```

The first line is the format tag. Blank lines and text after `#` are ignored.
`root` may contain `${NAME}` references to environment variables; an unset
variable is an error that names the variable. Relative file names are taken
relative to `root`.

`group <kind> <file> [options]` reads all sensors of a file, or the labels in
`select=` (comma separated, exact labels). Options are `turns=` (signal
scale), `area=` (total effective area for Rogowski files without a fourth
column) and `rphiz` (B-probe positions given as R, phi, z).

| kind | file format | signal |
|---|---|---|
| `flux_loop` | DIAGNO flux loop | loop integral of A |
| `saddle_coil` | DIAGNO flux loop | loop integral of A |
| `compensation_loop` | DIAGNO flux loop | loop integral of A |
| `diamagnetic_loop` | DIAGNO flux loop | loop integral of A plus the linked plasma toroidal flux (DIAGNO idia = 1) |
| `segmented_rogowski` | DIAGNO segmented Rogowski | sum of effective area times B along the path |
| `bprobe` | DIAGNO B-probe | effective area times B.n |

Every signal is multiplied by its `turns`. The kind fixes whether a loop links
the plasma poloidally: `diamagnetic_loop` sets idia = 1, and any other loop
kind rejects a file entry with idia = 1. DIAGNO differences (idia < 0) are
rejected because a selection would break their indices; write them as
`combine` lines instead.

`combine <label> <c1> <sensor1> [<c2> <sensor2> ...]` adds the output
`sum_k c_k raw_k`. This covers compensated diamagnetic loops, where the
compensation-loop signals are subtracted with a calibration factor, and any
other fixed linear combination. Calibration factors that change per
experimental programme belong in the private data next to the measurements,
not in this repository.

Label clashes, unknown sensors in a combination and unknown kinds or options
are reported with the manifest line. `tiago_cli diag lint <manifest> --kind
sensors` prints the counts per kind and runs the loop and Rogowski checks.

## Signals and derivatives

`tiago_sensor_signals` evaluates the raw signals with the existing vacuum and
virtual-casing evaluators, in the order loops, Rogowski coils, probes, and
returns `y = C raw`. The map is linear, so its JVP is `C` and its VJP is `C^T`
(`sensor_set_t%combine`, `%combine_vjp`, `%combine_matrix`). The Jacobians are
exact:

- `sensor_set_coil_jacobian`: per unit EXTCUR of every coil group;
- `sensor_set_plasma_jacobian`: per VMEC boundary-field coefficient and
  optionally per boundary-shape coefficient (the existing reverse mode);
- `sensor_set_linked_flux_row`: per unit linked toroidal flux (phiedge).

## Verification

`tiago_sensor_set` builds a synthetic set (geometry written by the test into a
directory named by an environment variable) and checks saddle, flux and
compensation loops against Maxwell's mutual inductance of coaxial circles, the
probe against the on-axis field of a circular coil, the diamagnetic loop
against zero vacuum flux plus its linked flux, the combination and its VJP,
the coil Jacobian against the linear signals, and the rejection messages.
`tiago_sensor_set_plasma` adds a VMEC plasma and checks that the plasma part
of every signal equals the mode Jacobian times the boundary coefficients plus
the linked-flux row. `tiago_w7x_sensor_set` lints a private manifest named by
`TIAGO_W7X_SENSORS` and is skipped when the variable is unset; the W7-X
geometry and its manifest live in the private plasma/data repository. The
evaluators themselves are those of the DIAGNO
comparison in `benchmarks/xdiagno`.
