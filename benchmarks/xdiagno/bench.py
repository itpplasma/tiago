#!/usr/bin/env python3
"""Accuracy and performance comparison of tiago_vacuum_cli against STELLOPT xdiagno.

Suites
  repo       the flux-loop / segmented-Rogowski cases shipped in tests/
  geometry   realistic generated sensor sets on the NCSX and M16N08 coil sets
  semantics  targeted DIAGNO-format features (open loops, iflflg, idia, EXTCUR, ...)
  plasma     plasma response of the NCSX VMEC equilibrium (xdiagno -vmec), plus an
             Ampere check against the VMEC toroidal current

Both codes always receive identical coil, EXTCUR and diagnostic files. Results
are printed as Markdown and written to _work/results/.
See README.md for setup.
"""
from __future__ import annotations

import argparse
import csv
import json
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
                  turn_files=None):
        out = d / "tiago"
        shutil.rmtree(out, ignore_errors=True)
        cmd = [self.tiago, str(coil), str(flux or ""), str(seg or ""),
               "--coil-extcur", str(d / "input."), "--seg-area", str(SEG_AREA),
               "--samples", str(ns), "--nfp", str(nfp), "--output-dir", str(out)]
        for option, path in (turn_files or {}).items():
            cmd += [option, str(path)]
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
        rows.append(dict(samples=ns, t_xdiagno_1=tx1, t_xdiagno_n=txn, t_tiago_1=tt1,
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
        semantic_check(r, "label_with_space", "#12", "label 'Loop A/upper'",
                       flux=[("Loop A/upper", closed)]),
        semantic_check(r, "segrog_point_area", "#12", "per-point eff_area 1e-4 / 5e-4",
                       seg_rows="     1\n     3     0     0 SEG_NONUNIF\n"
                                " 0.5 0.5 -0.4 1.0e-4\n 0.5 0.5 0.0 5.0e-4\n 0.5 0.5 0.4 0.0\n"),
    ]


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
    inp = (DATA / "input.ncsx").read_text()
    inp = re.sub(r"&DIAGNO_IN.*?/", "&DIAGNO_IN\n NU = 128\n NV = 32\n units = 1.\n"
                 " int_type = 'midpoint'\n int_step = 4\n flux_diag_file = 'flux.diagno'\n"
                 " seg_rog_file = 'seg.diagno'\n vc_adapt_tol = 1.0E-6\n vc_adapt_rel = 1.0E-5\n/",
                 inp, flags=re.S | re.I)
    (d / "input.ncsx").write_text(inp)

    cmd = ["mpirun", "--oversubscribe", "-np", str(r.ncpu), r.xdiagno,
           "-vmec", "ncsx", "-coil", "zero.coils", "-noverb"]
    t0 = time.perf_counter()
    subprocess.run(cmd, cwd=d, env=ENV, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t_xd = time.perf_counter() - t0
    x = merged(read_diagno_out(d / "diagno_flux.ncsx"), read_diagno_out(d / "diagno_seg.ncsx"))

    rows = []
    for grid in (32, 64):
        out = d / f"tiago_{grid}"
        t0 = time.perf_counter()
        subprocess.run([r.tiago, "", str(d / "flux.diagno"), str(d / "seg.diagno"),
                        "--plasma-wout", str(d / "wout_ncsx.nc"), "--plasma-nphi", str(grid),
                        "--plasma-ntheta", str(grid), "--samples", "4", "--output-dir", str(out)],
                       env=dict(ENV, OMP_NUM_THREADS=str(r.ncpu)), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        dt = time.perf_counter() - t0
        tv = merged(read_tiago_csv(out / "tiago_flux.csv"), read_tiago_csv(out / "tiago_segrog.csv"))
        rows.append(dict(grid=grid, t_tiago_n=dt,
                         flux=compare({k: v for k, v in tv.items() if k.startswith("flux")},
                                      {k: v for k, v in x.items() if k.startswith("flux")}),
                         seg=compare({k: v for k, v in tv.items() if k.startswith("seg")},
                                     {k: v for k, v in x.items() if k.startswith("seg")}),
                         ampere_tiago=tv["seg:AMPERE"]))
    area = SEG_AREA / 100   # AMPERE has 101 points, eff_area per segment
    return dict(case="plasma_ncsx", nflux=len(flux), nseg=len(seg), t_xdiagno_n=t_xd,
                ampere_exact=4e-7 * np.pi * abs(ctor) * area, ampere_xdiagno=x["seg:AMPERE"],
                rows=rows)


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
    if "plasma" in results:
        p = results["plasma"]
        out += ["### Plasma response (NCSX, plasma only)", "",
                f"{p['nflux']} flux loops and {p['nseg']} segmented Rogowskis; xdiagno -vmec "
                f"(adaptive virtual casing, tol 1e-6) on {ncpu} ranks took {p['t_xdiagno_n']:.0f} s. "
                "Tiago uses the VMEC boundary sheet current on a grid of `grid` x `grid` points "
                "per field period.", "",
                "| grid | Tiago n threads | flux: median / max vs xdiagno | Rogowski: median / max vs xdiagno |",
                "|---:|---:|---|---|"]
        for row in p["rows"]:
            out.append(f"| {row['grid']} | {row['t_tiago_n']:.2f} | {fmt(row['flux']['med'])} / "
                       f"{fmt(row['flux']['max'])} | {fmt(row['seg']['med'])} / {fmt(row['seg']['max'])} |")
        out += ["", f"Ampere loop around the plasma (x eff_area): mu0 I_tor = {p['ampere_exact']:.6e}, "
                f"xdiagno {p['ampere_xdiagno']:.6e}, Tiago {p['rows'][-1]['ampere_tiago']:.6e}.", ""]
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
    ap.add_argument("suites", nargs="*", help="repo, geometry, semantics, plasma (default: all)")
    ap.add_argument("--tiago", default=str(ROOT / "build/tiago_vacuum_cli"))
    ap.add_argument("--xdiagno", default=str(WORK / "bin/xdiagno"))
    ap.add_argument("--xdiagno-patched", default=str(WORK / "bin/xdiagno_patched"))
    ap.add_argument("--ncpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--quick", action="store_true", help="one sample count (6) instead of 2/6/16")
    args = ap.parse_args()
    args.suites = args.suites or ["repo", "geometry", "semantics", "plasma"]
    unknown = set(args.suites) - {"repo", "geometry", "semantics", "plasma"}
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

    md = report(results, args.ncpu)
    res_dir = WORK / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    (res_dir / "results.json").write_text(json.dumps(results, indent=1, default=str))
    (res_dir / "results.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
