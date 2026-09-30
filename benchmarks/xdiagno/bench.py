#!/usr/bin/env python3
"""Accuracy and performance comparison of tiago_vacuum_cli against STELLOPT xdiagno.

Suites
  repo       the flux-loop / segmented-Rogowski cases shipped in tests/
  geometry   realistic generated sensor sets on the NCSX and M16N08 coil sets
  semantics  targeted DIAGNO-format features (open loops, iflflg, idia, EXTCUR, ...)
  plasma     plasma response of the NCSX VMEC equilibrium (xdiagno -vmec), plus an
             Ampere check against the VMEC toroidal current
  features   magnetic probes and per-coil-group response matrices (xdiagno -mutual)
  equilibrium  derivatives of plasma signals with respect to VMEC input parameters
             (finite differences over xvmec2000 runs; needs build_xdiagno.sh --vmec)

Both codes always receive identical coil, EXTCUR and diagnostic files. Results
are printed as Markdown and written to _work/results/.
See README.md for setup.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "_work"
DATA = WORK / "data"
SEG_AREA = 3.4e-4
EQ_STEP = 1.0e-4   # relative finite-difference step of the equilibrium suite
AGREE = 1.0e-5          # relative tolerance for "codes agree"
ENV = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")


# --------------------------------------------------------------------- I/O
def write_diag(path: Path, entries, with_area: bool, flags=None) -> None:
    """DIAGNO format: (I6) count, (3I6,A48) header, then x y z [eff_area] rows."""
    with path.open("w") as f:
        f.write(f"{len(entries):6d}\n")
        for i, (label, pts) in enumerate(entries):
            ifl, idia = flags[i] if flags else (0, 0)
            f.write(f"{len(pts):6d}{ifl:6d}{idia:6d} {label:<48s}\n")
            for p in pts:
                area = f" {SEG_AREA / (len(pts) - 1): .12E}" if with_area else ""
                f.write(f" {p[0]: .12E} {p[1]: .12E} {p[2]: .12E}{area}\n")


def read_diag(path: Path):
    lines = path.read_text().splitlines()
    n, i, out = int(lines[0].split()[0]), 1, []
    for _ in range(n):
        head = lines[i].split(None, 3)
        npts = int(head[0])
        pts = [list(map(float, lines[i + 1 + j].split()[:3])) for j in range(npts)]
        out.append((head[3].strip(), pts, (int(head[1]), int(head[2]))))
        i += npts + 1
    return out


def extcur_from_coils(coil: Path, out: Path) -> None:
    """EXTCUR(g) = first non-zero current of group g, i.e. reproduce the file currents."""
    cur, buf = {}, []
    for line in coil.read_text().splitlines():
        s = line.strip().lower()
        if not s or s.startswith(("periods", "begin", "mirror")):
            continue
        if s.startswith("end"):
            break
        v = line.split()
        buf.append(float(v[3]))
        if len(v) >= 5:
            cur.setdefault(int(v[4]), next((x for x in buf if x != 0.0), 0.0))
            buf = []
    out.write_text("&INDATA\n" + "".join(
        f"  EXTCUR({g}) = {cur[g]:.15E}\n" for g in sorted(cur)) + "/\n")


def read_diagno_out(path: Path) -> dict:
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    n = int(lines[0])
    vals = [float(x) for x in lines[1:1 + n]]
    labels = []
    for l in lines[1 + n:1 + 2 * n]:
        parts = l.strip().split(None, 1)
        labels.append(parts[1].strip() if len(parts) == 2 and parts[0].isdigit() else l.strip())
    return dict(zip(labels, vals))


def read_tiago_csv(path: Path) -> dict:
    with path.open() as f:
        return {r["label"]: float(r["value"]) for r in csv.DictReader(f)}


# ----------------------------------------------------------------- runners
class Runner:
    def __init__(self, args):
        self.tiago = str(Path(args.tiago).resolve())
        self.xdiagno = str(Path(args.xdiagno).resolve())
        self.xdiagno_patched = args.xdiagno_patched if Path(args.xdiagno_patched).exists() else None
        self.ncpu = args.ncpu
        self.xvmec = str(Path(args.xvmec).resolve())

    def xdiagno_run(self, d: Path, coil: Path, flux, seg, ns: int, nproc: int = 1, binary=None,
                    turns=None):
        turn_lines = "".join(f"  {name} = {', '.join(f'{v:.12g}' for v in values)},\n"
                             for name, values in (turns or {}).items())
        (d / "diagno.control").write_text(
            "&diagno_in\n"
            f"  flux_diag_file = '{flux or ''}',\n"
            f"  seg_rog_file = '{seg or ''}',\n"
            "  nu = 64, nv = 64, int_type = 'midpoint',\n"
            f"  int_step = {ns},\n"
            "  lrphiz = .false., lvc_field = .false., luse_extcur = .true., units = 1.0,\n"
            + turn_lines + "/\n")
        for f in ("diagno_flux.", "diagno_seg."):
            (d / f).unlink(missing_ok=True)
        exe = [binary or self.xdiagno]
        cmd = (exe if nproc == 1 else ["mpirun", "--oversubscribe", "-np", str(nproc)] + exe) \
            + ["-vac", "-coil", str(coil), "-noverb"]
        t0 = time.perf_counter()
        subprocess.run(cmd, cwd=d, env=ENV, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        dt = time.perf_counter() - t0
        fl = read_diagno_out(d / "diagno_flux.") if flux else {}
        sg = read_diagno_out(d / "diagno_seg.") if seg else {}
        return dt, fl, sg

    def tiago_run(self, d: Path, coil: Path, flux, seg, ns: int, nthreads: int = 1, nfp: int = 1,
                  turn_files=None, gauss=False):
        out = d / "tiago"
        shutil.rmtree(out, ignore_errors=True)
        cmd = [self.tiago, "--coils", str(coil), "--coil-extcur", str(d / "input."),
               "--seg-area", str(SEG_AREA),
               "--samples", str(ns), "--nfp", str(nfp), "--output-dir", str(out)]
        if flux:
            cmd += ["--flux", str(flux)]
        if seg:
            cmd += ["--segrog", str(seg)]
        for option, path in (turn_files or {}).items():
            cmd += [option, str(path)]
        if gauss:
            cmd.append("--gauss")
        t0 = time.perf_counter()
        p = subprocess.run(cmd, env=dict(ENV, OMP_NUM_THREADS=str(nthreads)),
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        dt = time.perf_counter() - t0
        if p.returncode != 0:
            raise RuntimeError(f"tiago failed: {p.stderr.strip()}")
        fl = read_tiago_csv(out / "tiago_flux.csv") if flux else {}
        sg = read_tiago_csv(out / "tiago_segrog.csv") if seg else {}
        return dt, fl, sg


# ---------------------------------------------------------------- metrics
def compare(tiago: dict, ref: dict):
    """Relative error per signal, normalised by max(|ref|, 1e-6 * max|ref|)."""
    common = [k for k in ref if k in tiago]
    finite = [k for k in common if np.isfinite(ref[k]) and np.isfinite(tiago[k])]
    scale = max([abs(ref[k]) for k in finite] or [1.0])
    err = np.array([abs(tiago[k] - ref[k]) / max(abs(ref[k]), 1e-6 * scale) for k in finite])
    return dict(n=len(ref), missing=len(set(ref) ^ set(tiago)),
                nonfinite=len(common) - len(finite),
                med=float(np.median(err)) if err.size else float("nan"),
                max=float(err.max()) if err.size else float("nan"))


def merged(a: dict, b: dict) -> dict:
    return {**{f"flux:{k}": v for k, v in a.items()}, **{f"seg:{k}": v for k, v in b.items()}}


def prepare(name: str) -> Path:
    d = WORK / "runs" / name
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    return d


def read_turns(path: Path) -> dict:
    turns = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("#"):
            turns[parts[0]] = float(parts[1])
    return turns


def run_case(r: Runner, name, coil, flux_entries, seg_entries, nfp=1, flags=None,
             samples=(2, 6, 16), use_patched=False, turn_files=None):
    d = prepare(name)
    shutil.copy(coil, d / "coils.in")
    coil = d / "coils.in"
    extcur_from_coils(coil, d / "input.")
    flux = seg = None
    if flux_entries:
        flux = d / "flux.diagno"
        write_diag(flux, flux_entries, False, flags)
    if seg_entries:
        seg = d / "seg.diagno"
        write_diag(seg, seg_entries, True)
    binary = r.xdiagno_patched if use_patched and r.xdiagno_patched else None
    turns = None
    if turn_files:   # DIAGNO wants arrays in file order, Tiago label -> scale files
        ft = read_turns(turn_files["--flux-turns"])
        st = read_turns(turn_files["--segrog-turns"])
        turns = {"flux_turns": [ft.get(l, 1.0) for l, _ in flux_entries],
                 "segrog_turns": [st.get(l, 1.0) for l, _ in seg_entries]}
    xrun = lambda ns, n: r.xdiagno_run(d, coil, flux, seg, ns, n, binary, turns)
    trun = lambda ns, n: r.tiago_run(d, coil, flux, seg, ns, n, nfp, turn_files)
    _, xf, xs = xrun(64, r.ncpu)
    _, tf, ts = trun(64, r.ncpu)
    x64, t64 = merged(xf, xs), merged(tf, ts)
    rows = []
    for ns in samples:
        tx1, xf, xs = xrun(ns, 1)
        txn, _, _ = xrun(ns, r.ncpu)
        tt1, tf, ts = trun(ns, 1)
        ttn, _, _ = trun(ns, r.ncpu)
        x, t = merged(xf, xs), merged(tf, ts)
        tg, gf, gs = r.tiago_run(d, coil, flux, seg, ns, r.ncpu, nfp, turn_files, gauss=True)
        rows.append(dict(samples=ns, t_tiago_gauss_n=tg, quad_err_tiago_gauss=compare(merged(gf, gs), t64)["med"],
                         t_xdiagno_1=tx1, t_xdiagno_n=txn, t_tiago_1=tt1,
                         t_tiago_n=ttn, vs_xdiagno=compare(t, x),
                         quad_err_xdiagno=compare(x, x64)["med"],
                         quad_err_tiago=compare(t, t64)["med"]))
    return dict(case=name, nflux=len(flux_entries or []), nseg=len(seg_entries or []),
                coil_points=sum(1 for l in coil.read_text().splitlines() if len(l.split()) >= 4),
                xdiagno_binary="patched" if binary else "stock", rows=rows)


# ------------------------------------------------------------------ suites
def suite_repo(r: Runner, quick: bool):
    samples = (6,) if quick else (2, 6, 16)
    t = ROOT / "tests"
    cases = [(n, t / "data" / f"coils_{c}.coils", t / "data/fluxloop_sample.diagno",
              t / "data/segrog_sample.diagno", 1)
             for n, c in (("sample", "sample"), ("varying_current", "varying"),
                          ("negative_current", "negative"))]
    cases += [("ncsx_nfp1", DATA / "coils.NCSX_nfp1", t / "cases/ncsx_nfp1/fluxloop.diagno",
               t / "cases/ncsx_nfp1/segrog.diagno", 1),
              ("ncsx_nfp3", DATA / "coils.NCSX", t / "cases/ncsx_nfp3/fluxloop.diagno",
               t / "cases/ncsx_nfp3/segrog.diagno", 3)]
    out = []
    for name, coil, flux, seg, nfp in cases:
        fl, sg = read_diag(flux), read_diag(seg)
        turn_files = None
        if (flux.parent / "flux_turns.csv").exists():
            turn_files = {"--flux-turns": flux.parent / "flux_turns.csv",
                          "--segrog-turns": flux.parent / "segrog_turns.csv"}
        out.append(run_case(r, f"repo_{name}", coil, [(l, p) for l, p, _ in fl],
                            [(l, p) for l, p, _ in sg], nfp=nfp,
                            flags=[f for *_, f in fl], samples=samples, use_patched=nfp > 1,
                            turn_files=turn_files))
    return out


def torus(R0, r, phi, th):
    rr = R0 + r * np.cos(th)
    return [rr * np.cos(phi), rr * np.sin(phi), r * np.sin(th)]


def sensor_set(R0, rs, nphi, nth, side=6):
    """Poloidal (diamagnetic) loops, toroidal loops, saddle loops, segmented Rogowskis."""
    flux, seg = [], []
    for k in range(nphi):
        ph = 2 * np.pi * k / nphi
        flux.append((f"DIA_{k:02d}", [torus(R0, rs, ph, t) for t in np.linspace(0, 2 * np.pi, 49)]))
    for j, th in enumerate(np.linspace(-np.pi / 2, np.pi / 2, 5)):
        flux.append((f"TOR_{j:02d}", [torus(R0, rs, p, th) for p in np.linspace(0, 2 * np.pi, 97)]))
    dph, dth = 0.8 * 2 * np.pi / nphi, 0.8 * 2 * np.pi / nth
    s = np.linspace(0, 1, side)[:-1]
    for k in range(nphi):
        for j in range(nth):
            p0, t0 = 2 * np.pi * k / nphi, 2 * np.pi * j / nth
            pts = ([torus(R0, rs, p0 + dph * u, t0) for u in s]
                   + [torus(R0, rs, p0 + dph, t0 + dth * u) for u in s]
                   + [torus(R0, rs, p0 + dph * (1 - u), t0 + dth) for u in s]
                   + [torus(R0, rs, p0, t0 + dth * (1 - u)) for u in s])
            flux.append((f"SAD_{k:02d}_{j:02d}", pts + [pts[0]]))
            ph = 2 * np.pi * (k + 0.5) / nphi
            seg.append((f"SEG_{k:02d}_{j:02d}",
                        [torus(R0, rs, ph, t0 + dth * u) for u in np.linspace(0, 1, 8)]))
    return flux, seg


def suite_geometry(r: Runner, quick: bool):
    samples = (6,) if quick else (2, 6, 16)
    specs = [("ncsx", DATA / "coils.NCSX", 1.44, 0.50, 12, 8),
             ("m16n08", DATA / "coils.M16N08", 3.00, 0.65, 16, 8)]
    return [run_case(r, f"geom_{n}", c, *sensor_set(R0, rs, nphi, nth), samples=samples)
            for n, c, R0, rs, nphi, nth in specs]


SQUARE = [[0.4, 0.4, 0.3], [0.6, 0.4, 0.3], [0.6, 0.6, 0.3], [0.4, 0.6, 0.3]]
BIG = [[0.3, 0.3, 0.3], [0.7, 0.3, 0.3], [0.7, 0.7, 0.3], [0.3, 0.7, 0.3]]
TWO_GROUPS = """periods 3
begin filament
mirror NIL
 0 0 0 1
 1 0 0 1
 1 1 0 1
 0 1 0 1
 0 0 0 0 1 A
 0 0 0.1 1
 1 0 0.1 1
 1 1 0.1 1
 0 1 0.1 1
 0 0 0.1 0 2 B
end
"""


def semantic_check(r: Runner, name, issue, what, flux=None, seg=None, flags=None,
                   coil_text=TWO_GROUPS, extcur=None, nfp=1, seg_rows=None, use_patched=False):
    d = prepare(f"sem_{name}")
    coil = d / "coils.in"
    coil.write_text(coil_text)
    if extcur is None:
        extcur_from_coils(coil, d / "input.")
    else:
        (d / "input.").write_text(extcur)
    fpath = spath = None
    if flux:
        fpath = d / "flux.diagno"
        write_diag(fpath, flux, False, flags)
    if seg_rows is not None:
        spath = d / "seg.diagno"
        spath.write_text(seg_rows)
    binary = r.xdiagno_patched if use_patched else None
    if use_patched and not binary:
        return dict(check=name, issue=issue, what=what, verdict="skipped (no xdiagno_patched)")
    _, xf, xs = r.xdiagno_run(d, coil, fpath, spath, 6, 1, binary)
    try:
        _, tf, ts = r.tiago_run(d, coil, fpath, spath, 6, 1, nfp)
    except RuntimeError as err:
        return dict(check=name, issue=issue, what=what, xdiagno=merged(xf, xs),
                    tiago=str(err), verdict="DIFFER")
    x, t = merged(xf, xs), merged(tf, ts)
    c = compare(t, x)
    ok = c["missing"] == 0 and c["nonfinite"] == 0 and c["max"] < AGREE
    return dict(check=name, issue=issue, what=what, xdiagno=x, tiago=t,
                verdict="agree" if ok else "DIFFER")


def suite_semantics(r: Runner):
    closed = SQUARE + [SQUARE[0]]
    period = [[2.9 * np.cos(p), 2.9 * np.sin(p), 0.0]
              for p in np.linspace(0, 2 * np.pi / 3, 20, endpoint=False)]
    dup_coil = TWO_GROUPS.replace(" 1 0 0 1\n", " 1 0 0 1\n 1 0 0 1\n", 1)
    return [
        semantic_check(r, "closed_loop", "-", "reference: closed square loop",
                       flux=[("SQ", closed)]),
        semantic_check(r, "open_polygon", "#7", "loop without repeated first point",
                       flux=[("SQ_OPEN", SQUARE)]),
        semantic_check(r, "iflflg_period", "#8", "iflflg=1 loop over one field period (nfp=3)",
                       flux=[("TOR_PERIOD", period)], flags=[(1, 0)], nfp=3, use_patched=True),
        semantic_check(r, "idia_plus1", "#9", "idia=1 on a horizontal loop",
                       flux=[("HORIZ", closed)], flags=[(0, 1)]),
        semantic_check(r, "idia_minus", "#9", "idia=-1: subtract flux of loop 1",
                       flux=[("SMALL", closed), ("BIG_MINUS_SMALL", BIG + [BIG[0]])],
                       flags=[(0, 0), (0, -1)]),
        semantic_check(r, "zero_length_coil_segment", "#10", "duplicated coil point",
                       flux=[("SQ", closed)], coil_text=dup_coil),
        semantic_check(r, "extcur_zero", "#11", "EXTCUR(2)=0 switches group B off",
                       flux=[("SQ", closed)], extcur="&INDATA\n EXTCUR(1) = 1.0\n EXTCUR(2) = 0.0\n/\n"),
        semantic_check(r, "extcur_array", "#11", "EXTCUR = 1.0, 5.0 (namelist array)",
                       flux=[("SQ", closed)], extcur="&INDATA\n EXTCUR = 1.0, 5.0\n/\n"),
        semantic_check(r, "extcur_slice", "#11", "EXTCUR(1:2) = 2.0 5.0 (array slice)",
                       flux=[("SQ", closed)], extcur="&INDATA\n EXTCUR(1:2) = 2.0 5.0\n/\n"),
        semantic_check(r, "extcur_partial", "#11", "only EXTCUR(2:) = 5.0 given",
                       flux=[("SQ", closed)], extcur="&INDATA\n EXTCUR(2:) = 5.0\n/\n"),
        semantic_check(r, "label_with_space", "#12", "label 'Loop A/upper'",
                       flux=[("Loop A/upper", closed)]),
        semantic_check(r, "segrog_point_area", "#12", "per-point eff_area 1e-4 / 5e-4",
                       seg_rows="     1\n     3     0     0 SEG_NONUNIF\n"
                                " 0.5 0.5 -0.4 1.0e-4\n 0.5 0.5 0.0 5.0e-4\n 0.5 0.5 0.4 0.0\n"),
    ]


def probe_set(R0, rs, n, seed=1):
    """Probes on the sensor torus with pseudo-random orientations (degrees)."""
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(n):
        ph, th = 2 * np.pi * k / n, 2 * np.pi * rng.random()
        rows.append(torus(R0, rs, ph, th) + [360 * rng.random(), 180 * rng.random(), 1.0e-3])
    return rows


def write_probes(path: Path, rows) -> None:
    path.write_text(f"{len(rows)}\n" + "".join(
        " ".join(f"{v: .12E}" for v in row) + "\n" for row in rows))


def read_probe_out(path: Path) -> dict:
    """diagno_bth.<id>: i x y z |B| signal."""
    out = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 6 and parts[0].isdigit():
            out[f"PROBE_{int(parts[0]):04d}"] = float(parts[5])
    return out


def read_mut(path: Path, labels) -> dict:
    """DIAGNO *_mut_file: 'nfl ncg' then 'i ig value' rows."""
    out = {}
    for line in path.read_text().splitlines()[1:]:
        parts = line.split()
        if len(parts) == 3:
            i, g = int(parts[0]), int(parts[1])
            out[(labels[i - 1], g)] = float(parts[2])
    return out


def suite_features(r: Runner):
    """B-probes (vacuum) and response matrices on the NCSX coil set."""
    import csv
    d = prepare("features_ncsx")
    coil = d / "coils.in"
    shutil.copy(DATA / "coils.NCSX", coil)
    extcur_from_coils(coil, d / "input.")
    flux, seg = sensor_set(1.44, 0.50, 6, 2)
    write_diag(d / "flux.diagno", flux, False)
    write_diag(d / "seg.diagno", seg, True)
    write_probes(d / "probes.diagno", probe_set(1.44, 0.50, 40))
    # Separate controls: naming any *_mut_file makes stock DIAGNO try to read
    # every mutual file, including the unnamed probe one (see patches/).
    common = ("&diagno_in\n  flux_diag_file = 'flux.diagno',\n  seg_rog_file = 'seg.diagno',\n"
              "  int_type = 'midpoint', int_step = 6,\n  luse_extcur = .true., units = 1.0,\n")

    def run(control, extra):
        (d / "diagno.control").write_text(common + control + "/\n")
        subprocess.run([r.xdiagno, "-vac", "-coil", "coils.in", "-noverb"] + extra, cwd=d, env=ENV,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    run("  bprobes_file = 'probes.diagno',\n", [])
    x_probe = read_probe_out(d / "diagno_bth.")
    run("  flux_mut_file = 'flux.mut',\n  rog_mut_file = 'seg.mut',\n", ["-mutual"])
    x_resp = {**{("flux", *k): v for k, v in read_mut(d / "flux.mut", [l for l, _ in flux]).items()},
              **{("segrog", *k): v for k, v in read_mut(d / "seg.mut", [l for l, _ in seg]).items()}}

    out = d / "tiago"
    subprocess.run([r.tiago, "--coils", str(coil), "--flux", str(d / "flux.diagno"),
                    "--segrog", str(d / "seg.diagno"),
                    "--coil-extcur", str(d / "input."), "--samples", "6",
                    "--bprobes", str(d / "probes.diagno"), "--response-out", "response.csv",
                    "--output-dir", str(out)], env=ENV, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t_probe = read_tiago_csv(out / "tiago_bprobes.csv")
    with (out / "response.csv").open() as f:
        t_resp = {(row["kind"], row["label"], int(row["group"])): float(row["value"])
                  for row in csv.DictReader(f)}
    t_resp = {k: v for k, v in t_resp.items() if k[0] != "bprobe"}
    # response columns must add up to the signals: sum_g M_g EXTCUR_g
    extcur = [float(v) for v in re.findall(r"=\s*(\S+)", (d / "input.").read_text())]
    t_flux = read_tiago_csv(out / "tiago_flux.csv")
    recon = {l: sum(t_resp[("flux", l, g + 1)] * extcur[g] for g in range(len(extcur))) for l in t_flux}
    return dict(case="features_ncsx", probes=compare(t_probe, x_probe),
                response=compare({"|".join(map(str, k)): v for k, v in t_resp.items()},
                                 {"|".join(map(str, k)): v for k, v in x_resp.items()}),
                n_response=len(x_resp), recon=compare(recon, t_flux))


def suite_plasma(r: Runner):
    """Plasma-only signals of the NCSX equilibrium (STELLOPT DIAGNO_TEST).

    Sensors sit on a torus of minor radius 0.85 m around R = 1.40 m, clear of
    the plasma; poloidal loops are diamagnetic (idia=1). Stock xdiagno needs a
    coil file with as many groups as the wout has EXTCUR values (DIAGNO bugs,
    see patches/); zero-current coils keep it out of the result.
    """
    import re
    import netCDF4
    d = prepare("plasma_ncsx")
    shutil.copy(DATA / "wout_ncsx.nc", d / "wout_ncsx.nc")
    wout = netCDF4.Dataset(d / "wout_ncsx.nc")
    nextcur = len(wout["extcur"][:]) if "extcur" in wout.variables else 1
    ctor = float(wout["ctor"][:])
    wout.close()
    (d / "zero.coils").write_text("periods 1\nbegin filament\nmirror NIL\n" + "".join(
        f" 100 {g} 0 0.0\n 101 {g} 0 0.0\n 100 {g} 0 0.0 {g} ZERO{g}\n"
        for g in range(1, nextcur + 1)) + "end\n")

    flux, seg = sensor_set(1.40, 0.85, 6, 4)
    t = np.linspace(0.0, 2.0 * np.pi, 101)
    seg.append(("AMPERE", [[1.40 + 0.9 * np.cos(x), 0.0, 0.9 * np.sin(x)] for x in t]))
    flags = [(0, 1) if label.startswith("DIA") else (0, 0) for label, _ in flux]
    write_diag(d / "flux.diagno", flux, False, flags)
    write_diag(d / "seg.diagno", seg, True)
    write_probes(d / "probes.diagno", probe_set(1.40, 0.85, 20))
    inp = (DATA / "input.ncsx").read_text()
    inp = re.sub(r"&DIAGNO_IN.*?/", "&DIAGNO_IN\n NU = 128\n NV = 32\n units = 1.\n"
                 " int_type = 'midpoint'\n int_step = 4\n flux_diag_file = 'flux.diagno'\n"
                 " seg_rog_file = 'seg.diagno'\n bprobes_file = 'probes.diagno'\n"
                 " vc_adapt_tol = 1.0E-6\n vc_adapt_rel = 1.0E-5\n/",
                 inp, flags=re.S | re.I)
    (d / "input.ncsx").write_text(inp)

    cmd = ["mpirun", "--oversubscribe", "-np", str(r.ncpu), r.xdiagno,
           "-vmec", "ncsx", "-coil", "zero.coils", "-noverb"]
    t0 = time.perf_counter()
    subprocess.run(cmd, cwd=d, env=ENV, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t_xd = time.perf_counter() - t0
    x = merged(read_diagno_out(d / "diagno_flux.ncsx"), read_diagno_out(d / "diagno_seg.ncsx"))
    x_probe = read_probe_out(d / "diagno_bth.ncsx")

    rows = []
    for grid in (32, 64):
        out = d / f"tiago_{grid}"
        t0 = time.perf_counter()
        subprocess.run([r.tiago, "--flux", str(d / "flux.diagno"), "--segrog", str(d / "seg.diagno"),
                        "--plasma-wout", str(d / "wout_ncsx.nc"), "--plasma-nphi", str(grid),
                        "--plasma-ntheta", str(grid), "--samples", "4",
                        "--bprobes", str(d / "probes.diagno"), "--output-dir", str(out)],
                       env=dict(ENV, OMP_NUM_THREADS=str(r.ncpu)), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        dt = time.perf_counter() - t0
        tv = merged(read_tiago_csv(out / "tiago_flux.csv"), read_tiago_csv(out / "tiago_segrog.csv"))
        rows.append(dict(grid=grid, t_tiago_n=dt,
                         flux=compare({k: v for k, v in tv.items() if k.startswith("flux")},
                                      {k: v for k, v in x.items() if k.startswith("flux")}),
                         seg=compare({k: v for k, v in tv.items() if k.startswith("seg")},
                                     {k: v for k, v in x.items() if k.startswith("seg")}),
                         probes=compare(read_tiago_csv(out / "tiago_bprobes.csv"), x_probe),
                         ampere_tiago=tv["seg:AMPERE"]))
    jacobian = plasma_jacobian_check(r, d, x_probe)
    area = SEG_AREA / 100   # AMPERE has 101 points, eff_area per segment
    return dict(case="plasma_ncsx", nflux=len(flux), nseg=len(seg), nprobe=20, t_xdiagno_n=t_xd,
                ampere_exact=4e-7 * np.pi * abs(ctor) * area, ampere_xdiagno=x["seg:AMPERE"],
                rows=rows, jacobian=jacobian)


def plasma_jacobian_check(r: Runner, d: Path, x_probe) -> dict:
    """--plasma-response-out: finite difference in bsupvmnc(0,0) of a perturbed wout."""
    import csv
    import netCDF4

    def run(wout, out, extra=()):
        subprocess.run([r.tiago, "--flux", str(d / "flux.diagno"), "--segrog", str(d / "seg.diagno"),
                        "--plasma-wout", str(wout), "--plasma-nphi", "32", "--plasma-ntheta", "32",
                        "--samples", "4", "--bprobes", str(d / "probes.diagno"),
                        "--output-dir", str(d / out), *extra], env=ENV, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return {**merged(read_tiago_csv(d / out / "tiago_flux.csv"),
                         read_tiago_csv(d / out / "tiago_segrog.csv")),
                **{"probe:" + k: v for k, v in read_tiago_csv(d / out / "tiago_bprobes.csv").items()}}

    base = run(d / "wout_ncsx.nc", "jac_base", ("--plasma-response-out", "jacobian.csv"))
    shutil.copy(d / "wout_ncsx.nc", d / "wout_pert.nc")
    wout = netCDF4.Dataset(d / "wout_pert.nc", "r+")
    mode = int(np.where((wout["xm_nyq"][:] == 0) & (wout["xn_nyq"][:] == 0))[0][0])
    delta = 1.0e-3 * abs(float(wout["bsupvmnc"][-1, mode]))
    wout["bsupvmnc"][-1, mode] += delta / 1.5   # boundary value = 1.5 b(ns) - 0.5 b(ns-1)
    wout.close()
    pert = run(d / "wout_pert.nc", "jac_pert")
    prefix = {"flux": "flux:", "segrog": "seg:", "bprobe": "probe:"}
    with (d / "jac_base" / "jacobian.csv").open() as f:
        jac = {prefix[row["kind"]] + row["label"]: float(row["value"]) for row in csv.DictReader(f)
               if row["coefficient"] == "bsupvmnc" and row["m"] == "0" and row["n"] == "0"}
    fd = {k: (pert[k] - base[k]) / delta for k in jac}
    return dict(n=len(jac), fd=compare(fd, jac))


def set_indata(text: str, name: str, value) -> str:
    """Set a one-line &INDATA entry, adding it if absent."""
    line = f"  {name} = {value}"
    pattern = rf"^[ \t]*{name}[ \t]*=.*$"
    if re.search(pattern, text, flags=re.M | re.I):
        return re.sub(pattern, line, text, count=1, flags=re.M | re.I)
    return re.sub(r"&INDATA", "&INDATA\n" + line, text, count=1, flags=re.I)


def boundary_coefficients(wout: Path) -> dict:
    """(coefficient, m, n) -> boundary value 1.5 b(ns) - 0.5 b(ns-1), as Tiago's Jacobian columns."""
    import netCDF4
    w = netCDF4.Dataset(wout)
    xm, xn = w["xm_nyq"][:], w["xn_nyq"][:]
    out = {}
    for name in ("bsupumnc", "bsupvmnc", "bsupumns", "bsupvmns"):
        if name in w.variables:
            b = w[name][:]
            for k in range(len(xm)):
                out[(name, int(round(xm[k])), int(round(xn[k])))] = 1.5 * b[-1, k] - 0.5 * b[-2, k]
    ctor, signgs = float(w["ctor"][:]), int(w["signgs"][:])
    w.close()
    return out, ctor, signgs


def suite_equilibrium(r: Runner):
    """Derivatives of plasma signals with respect to VMEC input parameters (#24).

    Central finite differences over fixed-boundary VMEC runs of the LI383
    low-resolution case (0.8 s per run). Checks: step independence (h vs h/2),
    the chain rule through the boundary-field Jacobian (exact at fixed boundary
    shape, plus d phiedge for idia = 1 loops) and Ampere (d signal / d p =
    mu0 * eff_area * d I_tor / d p).
    """
    if not Path(r.xvmec).is_file():
        return dict(skipped=f"{r.xvmec} not found; build it with ./build_xdiagno.sh --vmec")
    d = prepare("equilibrium_li383")
    base_text = (DATA / "input.li383_low_res").read_text()
    base_text = re.sub(r"^\s*LWOUTTXT.*$\n?", "", base_text, flags=re.M | re.I)
    base_text = set_indata(base_text, "FTOL_ARRAY", "1.0E-14")
    base_text = set_indata(base_text, "NITER", "20000")

    R0, rs = 1.42, 0.85
    flux, seg = sensor_set(R0, rs, 6, 4)
    t = np.linspace(0.0, 2.0 * np.pi, 101)
    seg.append(("AMPERE", [[R0 + 0.9 * np.cos(x), 0.0, 0.9 * np.sin(x)] for x in t]))
    flags = [(0, 1) if label.startswith("DIA") else (0, 0) for label, _ in flux]
    write_diag(d / "flux.diagno", flux, False, flags)
    write_diag(d / "seg.diagno", seg, True)
    write_probes(d / "probes.diagno", probe_set(R0, rs, 20))
    dia = {f"flux:{label}" for label, _ in flux if label.startswith("DIA")}

    def vmec(tag: str, text: str) -> Path:
        (d / f"input.{tag}").write_text(text)
        subprocess.run(["mpirun", "-np", "1", r.xvmec, tag], cwd=d, env=ENV, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return d / f"wout_{tag}.nc"

    def signals(wout: Path, tag: str, extra=()) -> dict:
        out = d / f"tiago_{tag}"
        p = subprocess.run([r.tiago, "--flux", str(d / "flux.diagno"), "--segrog", str(d / "seg.diagno"),
                            "--bprobes", str(d / "probes.diagno"), "--plasma-wout", str(wout),
                            "--samples", "4", "--gauss", "--output-dir", str(out), *extra],
                           env=ENV, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                           text=True)
        if "WARNING" in p.stderr:
            raise RuntimeError(f"sensor too close to the plasma: {p.stderr.strip()}")
        return {**merged(read_tiago_csv(out / "tiago_flux.csv"), read_tiago_csv(out / "tiago_segrog.csv")),
                **{"probe:" + k: v for k, v in read_tiago_csv(out / "tiago_bprobes.csv").items()}}

    t0 = time.perf_counter()
    base_wout = vmec("base", base_text)
    t_vmec = time.perf_counter() - t0
    base = signals(base_wout, "base", ("--plasma-response-out", "jacobian.csv"))
    prefix = {"flux": "flux:", "segrog": "seg:", "bprobe": "probe:"}
    jac = {}
    with (d / "tiago_base" / "jacobian.csv").open() as f:
        for row in csv.DictReader(f):
            jac.setdefault(prefix[row["kind"]] + row["label"], {})[
                (row["coefficient"], int(row["m"]), int(row["n"]))] = float(row["value"])
    _, ctor0, signgs = boundary_coefficients(base_wout)
    area = SEG_AREA / 100   # AMPERE: 101 points, eff_area per segment
    orient = np.sign(base["seg:AMPERE"] / (4e-7 * np.pi * ctor0 * area))
    to_current = orient * 4e-7 * np.pi * area   # AMPERE signal per ampere enclosed

    values = {name: float(re.search(rf"^\s*{name}\s*=\s*([-+.\dEeDd]+)", base_text, re.M | re.I)
                          .group(1).replace("D", "E").replace("d", "e"))
              if re.search(rf"^\s*{name}\s*=", base_text, re.M | re.I) else 1.0
              for name in ("PHIEDGE", "CURTOR", "PRES_SCALE")}
    rows, derivative = [], {}
    for name, v0 in values.items():
        h = EQ_STEP * abs(v0)
        fd, coef_fd, ctor_fd = {}, {}, {}
        for step in (h, h / 2):
            runs = {}
            for sgn in (1, -1):
                tag = f"{name.lower()}_{'p' if sgn > 0 else 'm'}{step:.3e}"
                wout = vmec(tag, set_indata(base_text, name, f"{v0 + sgn * step:.16E}"))
                runs[sgn] = (signals(wout, tag), *boundary_coefficients(wout)[:2])
            fd[step] = {k: (runs[1][0][k] - runs[-1][0][k]) / (2 * step) for k in base}
            coef_fd[step] = {c: (runs[1][1][c] - runs[-1][1][c]) / (2 * step) for c in runs[1][1]}
            ctor_fd[step] = (runs[1][2] - runs[-1][2]) / (2 * step)
        dS = fd[h / 2]
        chain = {k: sum(jac[k][c] * coef_fd[h / 2][c] for c in jac[k])
                 + (signgs if (name == "PHIEDGE" and k in dia) else 0.0) for k in dS}
        rows.append(dict(parameter=name, value=v0, step=h, steps=compare(fd[h], dS),
                         chain=compare(chain, dS), dI_tiago=dS["seg:AMPERE"] / to_current,
                         dI_vmec=ctor_fd[h / 2]))
        derivative[name] = dS

    # Resolution study: the plasma current seen by the Ampere loop vs VMEC's ctor.
    curtor = values["CURTOR"]
    convergence = []
    for ns in (16, 32, 64, 128):
        wout = vmec(f"ns{ns}", set_indata(base_text, "NS_ARRAY", str(ns)))
        s = signals(wout, f"ns{ns}")
        convergence.append(dict(ns=ns, vmec=boundary_coefficients(wout)[1] - curtor,
                                tiago=s["seg:AMPERE"] / to_current - curtor))
    angular = []
    for mpol, ntor in ((4, 3), (6, 4), (8, 6)):
        text = set_indata(set_indata(set_indata(base_text, "NS_ARRAY", "64"), "MPOL", mpol), "NTOR", ntor)
        wout = vmec(f"mpol{mpol}", text)
        s = signals(wout, f"mpol{mpol}")
        angular.append(dict(mpol=mpol, ntor=ntor,
                            tiago=s["seg:AMPERE"] / to_current - boundary_coefficients(wout)[1]))

    res_dir = WORK / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    with (res_dir / "equilibrium_jacobian.csv").open("w") as f:
        f.write("signal,parameter,value\n")
        for name, dS in derivative.items():
            for k, v in dS.items():
                f.write(f"{k},{name},{v:.16e}\n")
    show = ["flux:DIA_00", "flux:TOR_02", "flux:SAD_00_00", "seg:SEG_00_00", "seg:AMPERE", "probe:PROBE_001"]
    show = [k for k in show if k in base] or list(base)[:6]
    sensitivity = {k: {n: values[n] * derivative[n][k] / base[k] for n in values} for k in show}
    return dict(case="li383_low_res", nsignals=len(base), t_vmec=t_vmec, nruns=1 + 4 * len(values),
                rows=rows, sensitivity=sensitivity, ctor=ctor0, curtor=curtor, convergence=convergence,
                angular=angular)


# ----------------------------------------------------------------- report
def fmt(x):
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.1e}"


def report(results, ncpu) -> str:
    out = []
    perf = [c for s in ("repo", "geometry") for c in results.get(s, [])]
    if perf:
        out += ["### Accuracy and performance", "",
                f"Wall time in seconds including coil loading; `n` = {ncpu} MPI ranks (xdiagno) "
                "or OpenMP threads (Tiago). *Tiago vs xdiagno* is the relative difference at "
                "equal `int_step`/`--samples`; *quadrature error* is each code's median "
                "deviation from its own run at 64 samples per segment.", "",
                "| case | coil pts | flux/seg | samples | xdiagno 1 | xdiagno n | Tiago 1 "
                "| Tiago n | Tiago vs xdiagno median / max | missing / non-finite "
                "| quad. err xdiagno / Tiago |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|"]
        for c in perf:
            for row in c["rows"]:
                v = row["vs_xdiagno"]
                out.append(
                    f"| {c['case']}{' (xdiagno_patched)' if c['xdiagno_binary'] == 'patched' else ''} "
                    f"| {c['coil_points']} | {c['nflux']}/{c['nseg']} | {row['samples']} "
                    f"| {row['t_xdiagno_1']:.2f} | {row['t_xdiagno_n']:.2f} "
                    f"| {row['t_tiago_1']:.2f} | {row['t_tiago_n']:.2f} "
                    f"| {fmt(v['med'])} / {fmt(v['max'])} | {v['missing']} / {v['nonfinite']} "
                    f"| {fmt(row['quad_err_xdiagno'])} / {fmt(row['quad_err_tiago'])} |")
        out.append("")
        geo = results.get("geometry", [])
        if geo:
            out += ["### Quadrature: midpoint vs Gauss-Legendre (Tiago, `--gauss`)", "",
                    "Median deviation from the converged 64-sample midpoint result, and wall "
                    f"time on {ncpu} threads, at equal points per segment.", "",
                    "| case | samples | midpoint error | midpoint time | Gauss error | Gauss time |",
                    "|---|---:|---:|---:|---:|---:|"]
            for c in geo:
                for row in c["rows"]:
                    out.append(f"| {c['case']} | {row['samples']} | {fmt(row['quad_err_tiago'])} | "
                               f"{row['t_tiago_n']:.2f} | {fmt(row['quad_err_tiago_gauss'])} | "
                               f"{row['t_tiago_gauss_n']:.2f} |")
            out.append("")
    if "features" in results:
        f = results["features"]
        out += ["### Magnetic probes and response matrices (NCSX coils, vacuum)", "",
                "| quantity | Tiago vs xdiagno median / max | missing / non-finite |", "|---|---|---|",
                f"| 40 B-probes | {fmt(f['probes']['med'])} / {fmt(f['probes']['max'])} "
                f"| {f['probes']['missing']} / {f['probes']['nonfinite']} |",
                f"| response matrix ({f['n_response']} entries, xdiagno -mutual) | "
                f"{fmt(f['response']['med'])} / {fmt(f['response']['max'])} "
                f"| {f['response']['missing']} / {f['response']['nonfinite']} |",
                f"| sum_g M_g EXTCUR_g vs Tiago signals | {fmt(f['recon']['med'])} / "
                f"{fmt(f['recon']['max'])} | {f['recon']['missing']} / {f['recon']['nonfinite']} |", ""]
    if "plasma" in results:
        p = results["plasma"]
        out += ["### Plasma response (NCSX, plasma only)", "",
                f"{p['nflux']} flux loops, {p['nseg']} segmented Rogowskis and {p['nprobe']} "
                "B-probes; xdiagno -vmec "
                f"(adaptive virtual casing, tol 1e-6) on {ncpu} ranks took {p['t_xdiagno_n']:.0f} s. "
                "Tiago uses the VMEC boundary sheet current on a grid of `grid` x `grid` points "
                "per field period.", "",
                "| grid | Tiago n threads | flux: median / max vs xdiagno | Rogowski: median / max "
                "| B-probe: median / max |", "|---:|---:|---|---|---|"]
        for row in p["rows"]:
            out.append(f"| {row['grid']} | {row['t_tiago_n']:.2f} | {fmt(row['flux']['med'])} / "
                       f"{fmt(row['flux']['max'])} | {fmt(row['seg']['med'])} / {fmt(row['seg']['max'])} "
                       f"| {fmt(row['probes']['med'])} / {fmt(row['probes']['max'])} |")
        out += ["", f"Ampere loop around the plasma (x eff_area): mu0 I_tor = {p['ampere_exact']:.6e}, "
                f"xdiagno {p['ampere_xdiagno']:.6e}, Tiago {p['rows'][-1]['ampere_tiago']:.6e}.", "",
                f"Boundary-field Jacobian (`--plasma-response-out`) vs a finite difference in "
                f"`bsupvmnc(0,0)` over all {p['jacobian']['n']} signals: median "
                f"{fmt(p['jacobian']['fd']['med'])}, max {fmt(p['jacobian']['fd']['max'])}.", ""]
    if "equilibrium" in results:
        e = results["equilibrium"]
        if "skipped" in e:
            out += ["### Equilibrium-parameter derivatives", "", f"Skipped: {e['skipped']}", ""]
        else:
            out += ["### Equilibrium-parameter derivatives (LI383 low resolution, plasma only)", "",
                    f"{e['nsignals']} signals; {e['nruns']} fixed-boundary VMEC runs of "
                    f"{e['t_vmec']:.1f} s each. Central differences with step h and h/2; *chain rule* "
                    "is the boundary-field Jacobian (`--plasma-response-out`) times the finite "
                    "difference of the boundary coefficients (plus d phiedge on idia = 1 loops); "
                    "*dI/dp* is the change of the current enclosed by the AMPERE loop (signal / "
                    "(mu0 eff_area)) next to that of VMEC's `ctor`.", "",
                    "| parameter | value | h | h vs h/2 median / max | chain rule median / max "
                    "| dI/dp Tiago / VMEC ctor |", "|---|---:|---:|---|---|---|"]
            for row in e["rows"]:
                out.append(f"| {row['parameter']} | {row['value']:.6g} | {row['step']:.2e} "
                           f"| {fmt(row['steps']['med'])} / {fmt(row['steps']['max'])} "
                           f"| {fmt(row['chain']['med'])} / {fmt(row['chain']['max'])} "
                           f"| {row['dI_tiago']:.5g} / {row['dI_vmec']:.5g} |")
            out += ["", f"Radial resolution: enclosed current minus CURTOR = {e['curtor']:.6g} A "
                    "(the converged value). VMEC's `ctor` extrapolates the covariant B_u to the "
                    "boundary; Tiago's sheet uses the extrapolated contravariant B^u, B^v.", "",
                    "| ns | VMEC ctor - CURTOR [A] | Tiago Ampere loop - CURTOR [A] |", "|---:|---:|---:|"]
            for c in e["convergence"]:
                out.append(f"| {c['ns']} | {c['vmec']:+.1f} | {c['tiago']:+.1f} |")
            out += ["", "The remaining difference between the two is angular truncation: the "
                    "extrapolated contravariant field does not conserve the sheet current exactly "
                    "(the current through a poloidal cross-section varies with phi), which "
                    "vanishes with MPOL/NTOR:", "",
                    "| MPOL / NTOR (ns 64) | Tiago Ampere loop - VMEC ctor [A] |", "|---|---:|"]
            for c in e["angular"]:
                out.append(f"| {c['mpol']} / {c['ntor']} | {c['tiago']:+.1f} |")
            params = [row["parameter"] for row in e["rows"]]
            out += ["", "Relative sensitivity (p / S) dS/dp of some signals:", "",
                    "| signal | " + " | ".join(params) + " |", "|---|" + "---:|" * len(params)]
            for k, s in e["sensitivity"].items():
                out.append(f"| {k} | " + " | ".join(f"{s[p]:+.3f}" for p in params) + " |")
            out.append("")
    if "semantics" in results:
        out += ["### DIAGNO-format semantics", "",
                "| check | issue | what | xdiagno | Tiago | verdict |", "|---|---|---|---|---|---|"]
        for c in results["semantics"]:
            def show(v):
                if isinstance(v, dict):
                    return "<br>".join(f"{k}={val:.4e}" for k, val in v.items())
                return str(v) if v is not None else ""
            out.append(f"| {c['check']} | {c['issue']} | {c['what']} | {show(c.get('xdiagno'))} "
                       f"| {show(c.get('tiago'))} | {c['verdict']} |")
        out.append("")
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("suites", nargs="*",
                    help="repo, geometry, semantics, features, plasma, equilibrium (default: all)")
    ap.add_argument("--tiago", default=str(ROOT / "build/tiago_vacuum_cli"))
    ap.add_argument("--xdiagno", default=str(WORK / "bin/xdiagno"))
    ap.add_argument("--xdiagno-patched", default=str(WORK / "bin/xdiagno_patched"))
    ap.add_argument("--xvmec", default=str(WORK / "bin/xvmec2000"))
    ap.add_argument("--ncpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--quick", action="store_true", help="one sample count (6) instead of 2/6/16")
    args = ap.parse_args()
    suites = ["repo", "geometry", "semantics", "features", "plasma", "equilibrium"]
    args.suites = args.suites or suites
    unknown = set(args.suites) - set(suites)
    if unknown:
        ap.error(f"unknown suite(s): {', '.join(sorted(unknown))}")

    r = Runner(args)
    results = {}
    if "repo" in args.suites:
        results["repo"] = suite_repo(r, args.quick)
    if "geometry" in args.suites:
        results["geometry"] = suite_geometry(r, args.quick)
    if "semantics" in args.suites:
        results["semantics"] = suite_semantics(r)
    if "plasma" in args.suites:
        results["plasma"] = suite_plasma(r)
    if "features" in args.suites:
        results["features"] = suite_features(r)
    if "equilibrium" in args.suites:
        results["equilibrium"] = suite_equilibrium(r)

    md = report(results, args.ncpu)
    res_dir = WORK / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    (res_dir / "results.json").write_text(json.dumps(results, indent=1, default=str))
    (res_dir / "results.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
