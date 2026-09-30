#!/usr/bin/env python3
"""Accuracy and performance comparison of tiago_vacuum_cli against STELLOPT xdiagno.

Suites
  repo       the flux-loop / segmented-Rogowski cases shipped in tests/
  geometry   realistic generated sensor sets on the NCSX and M16N08 coil sets
  semantics  targeted DIAGNO-format features (open loops, iflflg, idia, EXTCUR, ...)

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
        self.xdiagno_vacnfp = args.xdiagno_vacnfp if Path(args.xdiagno_vacnfp).exists() else None
        self.ncpu = args.ncpu

    def xdiagno_run(self, d: Path, coil: Path, flux, seg, ns: int, nproc: int = 1, binary=None):
        (d / "diagno.control").write_text(
            "&diagno_in\n"
            f"  flux_diag_file = '{flux or ''}',\n"
            f"  seg_rog_file = '{seg or ''}',\n"
            "  nu = 64, nv = 64, int_type = 'midpoint',\n"
            f"  int_step = {ns},\n"
            "  lrphiz = .false., lvc_field = .false., luse_extcur = .true., units = 1.0,\n/\n")
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

    def tiago_run(self, d: Path, coil: Path, flux, seg, ns: int, nthreads: int = 1, nfp: int = 1):
        out = d / "tiago"
        shutil.rmtree(out, ignore_errors=True)
        cmd = [self.tiago, str(coil), str(flux or ""), str(seg or ""),
               "--coil-extcur", str(d / "input."), "--seg-area", str(SEG_AREA),
               "--samples", str(ns), "--nfp", str(nfp), "--output-dir", str(out)]
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


def run_case(r: Runner, name, coil, flux_entries, seg_entries, nfp=1, flags=None,
             samples=(2, 6, 16), vacnfp=False):
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
    binary = r.xdiagno_vacnfp if vacnfp and r.xdiagno_vacnfp else None
    _, xf, xs = r.xdiagno_run(d, coil, flux, seg, 64, r.ncpu, binary)
    _, tf, ts = r.tiago_run(d, coil, flux, seg, 64, r.ncpu, nfp)
    x64, t64 = merged(xf, xs), merged(tf, ts)
    rows = []
    for ns in samples:
        tx1, xf, xs = r.xdiagno_run(d, coil, flux, seg, ns, 1, binary)
        txn, _, _ = r.xdiagno_run(d, coil, flux, seg, ns, r.ncpu, binary)
        tt1, tf, ts = r.tiago_run(d, coil, flux, seg, ns, 1, nfp)
        ttn, _, _ = r.tiago_run(d, coil, flux, seg, ns, r.ncpu, nfp)
        x, t = merged(xf, xs), merged(tf, ts)
        rows.append(dict(samples=ns, t_xdiagno_1=tx1, t_xdiagno_n=txn, t_tiago_1=tt1,
                         t_tiago_n=ttn, vs_xdiagno=compare(t, x),
                         quad_err_xdiagno=compare(x, x64)["med"],
                         quad_err_tiago=compare(t, t64)["med"]))
    return dict(case=name, nflux=len(flux_entries or []), nseg=len(seg_entries or []),
                coil_points=sum(1 for l in coil.read_text().splitlines() if len(l.split()) >= 4),
                xdiagno_binary="vacnfp" if binary else "stock", rows=rows)


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
        out.append(run_case(r, f"repo_{name}", coil, [(l, p) for l, p, _ in fl],
                            [(l, p) for l, p, _ in sg], nfp=nfp,
                            flags=[f for *_, f in fl], samples=samples, vacnfp=nfp > 1))
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
                   coil_text=TWO_GROUPS, extcur=None, nfp=1, seg_rows=None, vacnfp=False):
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
    binary = r.xdiagno_vacnfp if vacnfp else None
    if vacnfp and not binary:
        return dict(check=name, issue=issue, what=what, verdict="skipped (no xdiagno_vacnfp)")
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
                       flux=[("TOR_PERIOD", period)], flags=[(1, 0)], nfp=3, vacnfp=True),
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
                    f"| {c['case']}{' (xdiagno_vacnfp)' if c['xdiagno_binary'] == 'vacnfp' else ''} "
                    f"| {c['coil_points']} | {c['nflux']}/{c['nseg']} | {row['samples']} "
                    f"| {row['t_xdiagno_1']:.2f} | {row['t_xdiagno_n']:.2f} "
                    f"| {row['t_tiago_1']:.2f} | {row['t_tiago_n']:.2f} "
                    f"| {fmt(v['med'])} / {fmt(v['max'])} | {v['missing']} / {v['nonfinite']} "
                    f"| {fmt(row['quad_err_xdiagno'])} / {fmt(row['quad_err_tiago'])} |")
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
    ap.add_argument("suites", nargs="*", help="repo, geometry, semantics (default: all)")
    ap.add_argument("--tiago", default=str(ROOT / "build/tiago_vacuum_cli"))
    ap.add_argument("--xdiagno", default=str(WORK / "bin/xdiagno"))
    ap.add_argument("--xdiagno-vacnfp", default=str(WORK / "bin/xdiagno_vacnfp"))
    ap.add_argument("--ncpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--quick", action="store_true", help="one sample count (6) instead of 2/6/16")
    args = ap.parse_args()
    args.suites = args.suites or ["repo", "geometry", "semantics"]
    unknown = set(args.suites) - {"repo", "geometry", "semantics"}
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

    md = report(results, args.ncpu)
    res_dir = WORK / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    (res_dir / "results.json").write_text(json.dumps(results, indent=1, default=str))
    (res_dir / "results.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
