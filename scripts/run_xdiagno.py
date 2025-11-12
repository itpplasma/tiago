#!/usr/bin/env python3
"""Cross-code harness against STELLOPT's xdiagno.

This script does the following:
1. Runs Tiago's vacuum solver CLI to dump flux-loop and segmented Rogowski
   predictions as CSV files.
2. If the environment variable TIAGO_XDIAGNO points to an executable, invokes
   it with the provided diagnostics to generate reference CSV files.
3. Compares both sets, writing a diff report inside the chosen output folder.

When TIAGO_XDIAGNO is unset, the harness prints a skip message and exits with
code 125 so CTest can mark the run as skipped.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import urllib.request
from collections import OrderedDict
from pathlib import Path
from time import perf_counter
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]

SKIP_EXIT_CODE = 125


def ensure_coil_file(path: Path, url: str | None) -> Tuple[Path, float]:
    if path.exists():
        payload = path.read_text()
    else:
        if not url:
            raise FileNotFoundError(f"coil file {path} missing; provide --coil-url")
        with urllib.request.urlopen(url) as source:
            payload = source.read().decode("utf-8")
    if payload_is_stellopt(payload):
        sanitized = payload if payload.endswith("\n") else payload + "\n"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(sanitized)
        return path, 0.0
    neo_payload, coil_current = convert_coil_payload(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(neo_payload)
    return path, coil_current


def convert_coil_payload(payload: str) -> Tuple[str, float]:
    lines = payload.splitlines()
    first_token = ""
    for raw in lines:
        stripped = raw.strip()
        if stripped:
            first_token = stripped.split()[0]
            break
    try:
        int(first_token)
    except ValueError:
        coords: List[List[float]] = []
        currents: List[float] = []
        for raw in lines:
            parts = raw.split()
            if len(parts) < 4:
                continue
            try:
                values = list(map(float, parts[:4]))
            except ValueError:
                continue
            coords.append(values)
            currents.append(values[3])
        if not coords:
            raise ValueError("unable to parse coil payload")
        output = [str(len(coords))]
        for x, y, z, current in coords:
            output.append(f"{x:.15E} {y:.15E} {z:.15E} {current:.15E}")
        mean_current = sum(currents) / len(currents)
        return "\n".join(output) + "\n", mean_current
    sanitized = payload if payload.endswith("\n") else payload + "\n"
    # Extract current from first numeric line for downstream metadata
    representative = 0.0
    for raw in lines:
        parts = raw.split()
        if len(parts) >= 4:
            try:
                representative = abs(float(parts[3]))
            except ValueError:
                continue
            break
    return sanitized, representative


def load_turns(path: str | None) -> Dict[str, float]:
    mapping: Dict[str, float] = {}
    if not path:
        return mapping
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(file_path)
    with file_path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cleaned = line.replace(",", " ")
            parts = [tok for tok in cleaned.split() if tok]
            if len(parts) < 2:
                continue
            label = parts[0]
            value = float(parts[1])
            mapping[label] = value
    return mapping


def build_turn_array(labels: List[str], mapping: Dict[str, float]) -> List[float]:
    if not mapping:
        return []
    return [mapping.get(label, 1.0) for label in labels]


def format_turn_list(values: List[float]) -> str:
    return ", ".join(f"{val:.12g}" for val in values)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coil", required=True, help="Path to coil geometry")
    parser.add_argument("--flux", required=True, help="Flux-loop diagnostics")
    parser.add_argument(
        "--segrog", required=True, help="Segmented Rogowski diagnostics"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Directory where CSV outputs and diff reports are stored",
    )
    parser.add_argument(
        "--tiago-bin",
        required=True,
        help="Path to tiago_vacuum_cli or equivalent Tiago forward solver",
    )
    parser.add_argument(
        "--seg-area",
        default="3.40e-4",
        help="Effective area (m^2) used for segmented Rogowski diagnostics",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=6,
        help="Samples per segment for both Tiago and xdiagno integration",
    )
    parser.add_argument(
        "--nfp",
        type=int,
        default=1,
        help="Field periods assumed by Tiago",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-3,
        help="Absolute tolerance for flux/voltage comparisons",
    )
    parser.add_argument(
        "--coil-extcur",
        default=None,
        help="Optional VMEC input or list supplying EXTCUR values",
    )
    parser.add_argument(
        "--rel-tolerance",
        type=float,
        default=1e-2,
        help="Relative tolerance; applies to |tiago-legacy| / max(|legacy|, 1e-30)",
    )
    parser.add_argument(
        "--coil-url",
        default=None,
        help="Optional URL used to download the coil file if it is missing",
    )
    parser.add_argument(
        "--label",
        default="case",
        help="Short label to differentiate artifact filenames",
    )
    parser.add_argument(
        "--flux-turns",
        default=None,
        help="Optional text file mapping flux-loop labels to turn counts",
    )
    parser.add_argument(
        "--segrog-turns",
        default=None,
        help="Optional text file mapping segmented Rogowski labels to turn counts",
    )
    parser.add_argument(
        "--legacy-flux",
        default=None,
        help=(
            "Optional path to xdiagno flux CSV if it doesn't emit the default "
            "diagno_flux.csv"
        ),
    )
    parser.add_argument(
        "--legacy-segrog",
        default=None,
        help="Optional path to xdiagno segmented Rogowski CSV",
    )
    parser.add_argument(
        "--extra-xdiagno-args",
        default=[],
        nargs=argparse.REMAINDER,
        help="Any trailing arguments passed verbatim to xdiagno",
    )
    return parser.parse_args()


def run_tiago(args: argparse.Namespace, out_dir: Path) -> Tuple[Path, Path, float]:
    flux_out = out_dir / "tiago_flux.csv"
    segrog_out = out_dir / "tiago_segrog.csv"
    cmd = [
        args.tiago_bin,
        args.coil,
        args.flux,
        args.segrog,
        "--output-dir",
        str(out_dir),
        "--flux-out",
        str(flux_out),
        "--segrog-out",
        str(segrog_out),
        "--seg-area",
        args.seg_area,
        "--samples",
        str(args.samples),
        "--nfp",
        str(args.nfp),
    ]
    if args.coil_extcur:
        cmd.extend(["--coil-extcur", args.coil_extcur])
    if args.flux_turns:
        cmd.extend(["--flux-turns", args.flux_turns])
    if args.segrog_turns:
        cmd.extend(["--segrog-turns", args.segrog_turns])
    start = perf_counter()
    subprocess.run(cmd, check=True)
    duration = perf_counter() - start
    return flux_out, segrog_out, duration


def run_xdiagno(
    args: argparse.Namespace,
    out_dir: Path,
    flux_turns_map: Dict[str, float],
    segrog_turns_map: Dict[str, float],
) -> Tuple[Path, Path, float]:
    binary = resolve_xdiagno()
    if not binary:
        print(
            "TIAGO_XDIAGNO not set and xdiagno missing from PATH; "
            "skipping cross-code comparison",
            file=sys.stderr,
        )
        sys.exit(SKIP_EXIT_CODE)
    flux_path, flux_labels = write_diagno_flux(out_dir, args.flux)
    seg_path, seg_labels = write_diagno_segrog(out_dir, args.segrog, float(args.seg_area))
    write_diagno_control(
        out_dir,
        flux_path,
        seg_path,
        args.samples,
        build_turn_array(flux_labels, flux_turns_map),
        build_turn_array(seg_labels, segrog_turns_map),
        args.nfp,
    )
    coil_path_obj = Path(args.coil)
    if coil_is_stellopt(coil_path_obj):
        coil_path = str(coil_path_obj)
    else:
        coil_path = str(write_diagno_coils(out_dir, args.coil))
    prepare_extcur_file(args, out_dir, coil_path_obj)
    cmd = [binary, "-vac", "-coil", coil_path, "-noverb"]
    if args.extra_xdiagno_args:
        cmd.extend(args.extra_xdiagno_args)
    start = perf_counter()
    try:
        subprocess.run(cmd, check=True, cwd=str(out_dir))
    except FileNotFoundError as err:
        print(f"failed to execute {binary}: {err}", file=sys.stderr)
        sys.exit(SKIP_EXIT_CODE)
    except subprocess.CalledProcessError as err:
        print(f"xdiagno exited with {err.returncode}", file=sys.stderr)
        raise
    duration = perf_counter() - start
    flux_raw = out_dir / "diagno_flux."
    seg_raw = out_dir / "diagno_seg."
    flux_csv = Path(args.legacy_flux) if args.legacy_flux else out_dir / "diagno_flux.csv"
    seg_csv = Path(args.legacy_segrog) if args.legacy_segrog else out_dir / "diagno_segrog.csv"
    convert_diagno_output(flux_raw, flux_csv, expect_index=True)
    convert_diagno_output(seg_raw, seg_csv, expect_index=False)
    return flux_csv, seg_csv, duration


def resolve_xdiagno() -> str | None:
    env_value = os.environ.get("TIAGO_XDIAGNO")
    candidates = []
    if env_value:
        candidates.append(env_value)
    candidates.append("xdiagno")
    for candidate in candidates:
        if not candidate:
            continue
        path = shutil.which(candidate)
        if path:
            return path
        maybe_path = Path(candidate)
        if maybe_path.exists():
            return str(maybe_path)
    return None


def load_csv(path: Path) -> OrderedDict[str, float]:
    if not path.exists():
        raise FileNotFoundError(path)
    results: OrderedDict[str, float] = OrderedDict()
    with path.open() as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            results[row[0].strip()] = float(row[1])
    return results


def compare_arrays(
    tiago_path: Path,
    legacy_path: Path,
    abs_tol: float,
    rel_tol: float,
    report: Path,
) -> Tuple[int, List[Tuple[str, float, float]]]:
    tiago = load_csv(tiago_path)
    legacy = load_csv(legacy_path)
    failures = 0
    pairs: List[Tuple[str, float, float]] = []
    with report.open("w") as handle:
        handle.write("label,tiago,legacy,diff\n")
        for label, value in tiago.items():
            legacy_value = legacy.get(label)
            if label not in legacy:
                failures += 1
                handle.write(f"{label},{value},MISSING,NaN\n")
                continue
            if legacy_value is None or math.isnan(legacy_value):
                handle.write(f"{label},{value},NaN,IGNORED_LEGACY_NAN\n")
                continue
            if value is None or math.isnan(value):
                failures += 1
                handle.write(f"{label},NaN,{legacy_value},NaN\n")
                continue
            pairs.append((label, value, legacy_value))
            diff = abs(value - legacy_value)
            scale = max(abs(legacy_value), 1.0e-30)
            rel_err = diff / scale if scale > 0.0 else math.inf
            if diff > abs_tol or rel_err > rel_tol:
                failures += 1
                handle.write(
                    f"{label},{value},{legacy_value},{diff}"
                    f",abs={diff},rel={rel_err}\n"
                )
        extra = set(legacy.keys()) - set(tiago.keys())
        for label in extra:
            failures += 1
            handle.write(f"{label},MISSING,{legacy[label]},NaN\n")
            pairs.append((label, math.nan, legacy[label]))
    return failures, pairs


def main() -> int:
    args = parse_args()
    out_dir = Path(args.output).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    coil_path, coil_current = ensure_coil_file(Path(args.coil).resolve(), args.coil_url)
    args.coil = str(coil_path)
    args.coil_current = coil_current
    coil_path_obj = Path(args.coil)
    if not args.coil_extcur and not coil_is_stellopt(coil_path_obj):
        vmec_input = write_vmec_input(out_dir, args.nfp, coil_path_obj)
        args.coil_extcur = str(vmec_input)
    flux_turns = load_turns(args.flux_turns)
    segrog_turns = load_turns(args.segrog_turns)
    tiago_flux, tiago_seg, tiago_time = run_tiago(args, out_dir)
    legacy_flux, legacy_seg, xdiagno_time = run_xdiagno(
        args, out_dir, flux_turns, segrog_turns
    )
    flux_report = out_dir / "flux_diff.csv"
    seg_report = out_dir / "segrog_diff.csv"
    failures = 0
    flux_failures, flux_pairs = compare_arrays(
        tiago_flux,
        legacy_flux,
        args.tolerance,
        args.rel_tolerance,
        flux_report,
    )
    seg_failures, seg_pairs = compare_arrays(
        tiago_seg,
        legacy_seg,
        args.tolerance,
        args.rel_tolerance,
        seg_report,
    )
    failures += flux_failures + seg_failures
    combined_pairs = [(
        f"flux:{label}", tiago_val, legacy_val
    ) for label, tiago_val, legacy_val in flux_pairs]
    combined_pairs.extend(
        (
            f"seg:{label}", tiago_val, legacy_val
        )
        for label, tiago_val, legacy_val in seg_pairs
    )
    generate_plot(
        combined_pairs,
        out_dir / f"diagnostics_{args.label}.png",
        tiago_time,
        xdiagno_time,
    )
    plot_geometry(
        Path(args.coil),
        Path(args.flux),
        Path(args.segrog),
        out_dir / f"geometry_{args.label}.png",
    )
    if failures > 0:
        print(
            f"Detected {failures} mismatched diagnostics; see {out_dir}",
            file=sys.stderr,
        )
        return 2
    print("Tiago and xdiagno outputs agree within tolerance")
    return 0


def write_diagno_control(
    out_dir: Path,
    flux_path: Path,
    seg_path: Path,
    samples: int,
    flux_turns: List[float] | None = None,
    segrog_turns: List[float] | None = None,
    nfp: int = 1,
) -> Path:
    control_path = out_dir / "diagno.control"
    lines = [
        "&diagno_in",
        f"  flux_diag_file = '{flux_path}',",
        f"  seg_rog_file = '{seg_path}',",
        "  nu = 64,",
        "  nv = 64,",
        "  int_type = 'midpoint',",
        f"  int_step = {max(1, samples)},",
        "  lrphiz = .false.,",
        "  lvc_field = .false.,",
        "  luse_extcur = .true.,",
        "  units = 1.0,",
    ]
    if flux_turns:
        lines.append(f"  flux_turns = {format_turn_list(flux_turns)},")
    if segrog_turns:
        lines.append(f"  segrog_turns = {format_turn_list(segrog_turns)},")
    lines.append("/" )
    lines.append("")
    control_path.write_text("\n".join(lines))
    return control_path


def write_vmec_input(out_dir: Path, nfp: int, coil_path: Path) -> Path:
    source = PROJECT_ROOT / "tests" / "data" / "input.diagno.stub"
    target = out_dir / "input."
    lines = []

    def floats_close(a: float, b: float) -> bool:
        scale = max(abs(a), abs(b), 1.0)
        return abs(a - b) <= scale * 1.0e-9 + 1.0e-12

    unique_currents: List[float] = []
    with coil_path.open() as handle:
        header = handle.readline().strip()
        count = int(header.split()[0])
        for _ in range(count):
            line = handle.readline()
            if not line:
                break
            values = line.split()
            if len(values) < 4:
                continue
            curr = float(values[3])
            if abs(curr) < 1.0e-10:
                continue
            if not any(floats_close(curr, seen) for seen in unique_currents):
                unique_currents.append(curr)

    if not unique_currents:
        unique_currents.append(1.0)
    ncurr = len(unique_currents)

    for raw in source.read_text().splitlines():
        stripped = raw.strip()
        if stripped.upper().startswith("NFP"):
            lines.append(f"  NFP = {max(1, nfp)}")
        elif stripped.upper().startswith("NCURR"):
            lines.append(f"  NCURR = {ncurr}")
        elif stripped.upper().startswith("EXTCUR"):
            for i, curr in enumerate(unique_currents, start=1):
                lines.append(f"  EXTCUR({i:2d}) = {curr:.12E}")
        else:
            lines.append(raw)
    target.write_text("\n".join(lines) + "\n")
    return target


def write_diagno_coils(out_dir: Path, source_path: str) -> Path:
    dest = out_dir / "coils.tiago"
    src = Path(source_path).resolve()

    def floats_close(a: float, b: float) -> bool:
        scale = max(abs(a), abs(b), 1.0)
        return abs(a - b) <= scale * 1.0e-9 + 1.0e-12

    def coords_close(p: Tuple[float, float, float], q: Tuple[float, float, float]) -> bool:
        return all(floats_close(a, b) for a, b in zip(p, q))

    with src.open() as handle:
        header = handle.readline().strip()
        try:
            count = int(header.split()[0])
        except ValueError as exc:
            raise ValueError(f"invalid Tiago coil header '{header}'") from exc
        points: List[Tuple[float, float, float, float]] = []
        for _ in range(count):
            line = handle.readline()
            if not line:
                break
            values = line.split()
            if len(values) < 4:
                continue
            x, y, z, current = map(float, values[:4])
            points.append((x, y, z, current))

    if not points:
        raise ValueError("No coil points found")

    loops: List[List[Tuple[float, float, float, float]]] = []
    loop_buffer: List[Tuple[float, float, float, float]] = []
    for point in points:
        if not loop_buffer:
            loop_buffer.append(point)
            continue
        current_value = loop_buffer[0][3]
        if not floats_close(point[3], current_value):
            loops.append(loop_buffer[:])
            loop_buffer = [point]
            continue
        loop_buffer.append(point)
        if coords_close(point[:3], loop_buffer[0][:3]):
            loops.append(loop_buffer[:])
            loop_buffer = []

    if loop_buffer:
        loop_buffer.append((loop_buffer[0][0], loop_buffer[0][1], loop_buffer[0][2], 0.0))
        loops.append(loop_buffer[:])

    if not loops:
        raise ValueError("Unable to derive any closed coils")

    unique_currents: List[float] = []
    for loop in loops:
        current_value = loop[0][3]
        if abs(current_value) < 1.0e-10:
            continue
        if not any(floats_close(current_value, seen) for seen in unique_currents):
            unique_currents.append(current_value)

    def find_group_id(value: float) -> int:
        for idx, seen in enumerate(unique_currents, start=1):
            if floats_close(value, seen):
                return idx
        # if everything cancelled out (e.g., zero current), keep group 1
        return 1

    with dest.open("w") as handle:
        handle.write("periods 1\n")
        handle.write("begin filament\n")
        handle.write("mirror NIL\n")

        for loop in loops:
            if len(loop) < 2:
                continue
            current_value = loop[0][3]
            group_id = find_group_id(current_value)
            label = f"TIAGO_G{group_id:02d}"
            for x, y, z, curr in loop[:-1]:
                handle.write(
                    f" {x: .15E} {y: .15E} {z: .15E} {curr: .15E}\n"
                )
            close_x, close_y, close_z, _ = loop[0]
            handle.write(
                f" {close_x: .15E} {close_y: .15E} {close_z: .15E}"
                f" {0.0: .15E} {group_id:4d} {label}\n"
            )

        handle.write("end\n")

    return dest


def prepare_extcur_file(args: argparse.Namespace, out_dir: Path, coil_path: Path) -> None:
    target = out_dir / "input."
    if args.coil_extcur:
        source = Path(args.coil_extcur).expanduser().resolve()
        if source != target.resolve():
            shutil.copyfile(source, target)
        args.coil_extcur = str(source)
        return
    if coil_is_stellopt(coil_path):
        raise RuntimeError(
            "STELLOPT coil files require --coil-extcur to supply EXTCUR values"
        )
    vmec_input = write_vmec_input(out_dir, args.nfp, coil_path)
    args.coil_extcur = str(vmec_input)


def coil_is_stellopt(path: Path) -> bool:
    try:
        with path.open() as handle:
            text = handle.read(512)
    except FileNotFoundError:
        return False
    return payload_is_stellopt(text)


def payload_is_stellopt(payload: str) -> bool:
    for raw in payload.splitlines():
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith("!"):
            continue
        return stripped.lower().startswith("periods")
    return False


def write_diagno_segrog(
    out_dir: Path, source_path: str, area_value: float
) -> Tuple[Path, List[str]]:
    dest = out_dir / "segrog.diagno"
    src = Path(source_path).resolve()
    with src.open() as handle:
        lines = [line.rstrip("\n") for line in handle]
    if not lines:
        raise ValueError("segrog file is empty")
    labels: List[str] = []
    with dest.open("w") as out:
        total = int(lines[0].split()[0])
        out.write(f"{total:6d}\n")
        idx = 1
        while idx < len(lines):
            header = lines[idx].strip()
            idx += 1
            if not header:
                continue
            parts = header.split(maxsplit=3)
            if len(parts) < 4:
                raise ValueError("invalid segmented Rogowski header")
            nseg = int(parts[0])
            ifl = int(parts[1])
            idia = int(parts[2])
            label = parts[3][:48]
            out.write(f"{nseg:6d}{ifl:6d}{idia:6d} {label:<48}\n")
            labels.append(label.strip())
            effective = area_value / max(1, nseg - 1)
            for _ in range(nseg):
                if idx >= len(lines):
                    raise ValueError("unexpected end of segmented Rogowski data")
                coords = lines[idx].split()
                idx += 1
                if len(coords) < 3:
                    continue
                x, y, z = map(float, coords[:3])
                out.write(
                    f" {x: .10E} {y: .10E} {z: .10E} {effective: .10E}\n"
                )
    return dest, labels


def write_diagno_flux(out_dir: Path, source_path: str) -> Tuple[Path, List[str]]:
    dest = out_dir / "fluxloop.diagno"
    src = Path(source_path).resolve()
    with src.open() as handle:
        lines = [line.rstrip("\n") for line in handle]
    if not lines:
        raise ValueError("flux file is empty")
    labels: List[str] = []
    with dest.open("w") as out:
        total = int(lines[0].split()[0])
        out.write(f"{total:6d}\n")
        idx = 1
        while idx < len(lines):
            header = lines[idx].strip()
            idx += 1
            if not header:
                continue
            parts = header.split(maxsplit=3)
            if len(parts) < 4:
                raise ValueError("invalid flux header")
            nseg = int(parts[0])
            ifl = int(parts[1])
            idia = int(parts[2])
            label = parts[3][:48]
            out.write(f"{nseg:6d}{ifl:6d}{idia:6d} {label:<48}\n")
            labels.append(label.strip())
            for _ in range(nseg):
                if idx >= len(lines):
                    raise ValueError("unexpected end of flux file")
                coords = lines[idx].split()
                idx += 1
                if len(coords) < 3:
                    continue
                x, y, z = map(float, coords[:3])
                out.write(f" {x: .10E} {y: .10E} {z: .10E}\n")
    return dest, labels


def convert_diagno_output(raw_path: Path, csv_path: Path, expect_index: bool) -> None:
    entries = parse_diagno_file(raw_path, expect_index)
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "value"])
        for label, value in entries:
            writer.writerow([label, value])


def parse_diagno_file(raw_path: Path, expect_index: bool) -> List[Tuple[str, float]]:
    if not raw_path.exists():
        raise FileNotFoundError(raw_path)
    with raw_path.open() as handle:
        lines = [line.rstrip() for line in handle if line.strip()]
    if not lines:
        raise ValueError(f"diagno output {raw_path} is empty")
    count = int(lines[0].split()[0])
    values: List[float] = []
    idx = 1
    for _ in range(count):
        if idx >= len(lines):
            raise ValueError(f"diagno output {raw_path} truncated")
        values.append(float(lines[idx].split()[0]))
        idx += 1
    labels: List[str] = []
    while idx < len(lines) and len(labels) < count:
        line = lines[idx].strip()
        idx += 1
        if expect_index:
            parts = line.split(None, 1)
            if len(parts) == 2 and parts[0].lstrip("+-").isdigit():
                labels.append(parts[1].strip())
                continue
        labels.append(line)
    if len(labels) < count:
        labels.extend(f"diag_{i+1}" for i in range(len(labels), count))
    return list(zip(labels[:count], values))


def generate_plot(
    pairs: List[Tuple[str, float, float]],
    path: Path,
    tiago_time: float,
    xdiagno_time: float,
) -> None:
    if not pairs:
        return

    labels = [label for label, _, _ in pairs]
    tiago_vals = [abs(val) if val is not None and not math.isnan(val) else math.nan for _, val, _ in pairs]
    legacy_vals = [abs(val) if val is not None and not math.isnan(val) else math.nan for _, _, val in pairs]

    rel_errors = []
    for _, tiago_val, legacy_val in pairs:
        if legacy_val is None or math.isnan(legacy_val) or tiago_val is None or math.isnan(tiago_val):
            rel_errors.append(math.nan)
        else:
            denom = max(abs(legacy_val), 1.0e-30)
            rel_errors.append(abs(tiago_val - legacy_val) / denom)

    fig, (ax_abs, ax_rel) = plt.subplots(
        2, 1, figsize=(max(6, len(pairs) * 0.7), 6), sharex=True
    )

    ax_abs.plot(labels, legacy_vals, color="#f68d40", marker="o", label="xdiagno")
    ax_abs.plot(labels, tiago_vals, color="#18b5aa", marker="x", label="tiago")
    ax_abs.set_ylabel("|Signal| [SI]")
    ax_abs.set_ylim(bottom=0.0)
    ax_abs.grid(True, linestyle="--", alpha=0.3)
    ax_abs.legend()

    ax_rel.bar(range(len(labels)), rel_errors, color="#7f8c8d", alpha=0.8)
    ax_rel.set_xticks(range(len(labels)))
    ax_rel.set_xticklabels(labels, rotation=45, ha="right")
    ax_rel.set_ylabel("Relative error")
    ax_rel.set_ylim(bottom=0.0)
    ax_rel.grid(True, linestyle="--", alpha=0.3)
    runtime_text = (
        f"tiago: {tiago_time*1e3:.1f} ms\n"
        f"xdiagno: {xdiagno_time*1e3:.1f} ms"
    )
    ax_rel.text(0.02, 0.95, runtime_text, transform=ax_rel.transAxes, va="top")

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_geometry(coil_path: Path, flux_path: Path, seg_path: Path, out_path: Path) -> None:
    try:
        coil_points = read_simple_coils(coil_path)
        flux_loops = read_diag_flux(flux_path)
        seg_loops = read_diag_segrog(seg_path)
    except Exception as exc:  # pragma: no cover
        print(f"geometry plot skipped: {exc}", file=sys.stderr)
        return

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    if coil_points:
        xs, ys, zs = zip(*coil_points)
        ax.plot(xs, ys, zs, color="#7f8c8d", label="Coil", linewidth=2)

    for label, pts in flux_loops:
        xs, ys, zs = zip(*pts)
        ax.plot(xs, ys, zs, label=f"Flux:{label}")

    for label, pts in seg_loops:
        xs, ys, zs = zip(*pts)
        ax.plot(xs, ys, zs, linestyle="--", label=f"Seg:{label}")

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    set_equal_3d(ax)
    ax.view_init(elev=20, azim=35)
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def read_simple_coils(path: Path) -> List[Tuple[float, float, float]]:
    with path.open() as handle:
        header = handle.readline().strip()
        if header.lower().startswith("periods"):
            points = []
            for line in handle:
                stripped = line.strip().lower()
                if not stripped:
                    continue
                if stripped.startswith("begin") or stripped.startswith("mirror"):
                    continue
                if stripped.startswith("end"):
                    break
                vals = line.split()
                if len(vals) < 3:
                    continue
                x, y, z = map(float, vals[:3])
                points.append((x, y, z))
            return points
        count = int(header.split()[0])
        points = []
        for _ in range(count):
            line = handle.readline()
            if not line:
                break
            vals = line.split()
            if len(vals) < 3:
                continue
            x, y, z = map(float, vals[:3])
            points.append((x, y, z))
    return points


def read_diag_flux(path: Path) -> List[Tuple[str, List[Tuple[float, float, float]]]]:
    loops = []
    with path.open() as handle:
        total = int(handle.readline().split()[0])
        for _ in range(total):
            header = handle.readline()
            if not header:
                break
            parts = header.split(maxsplit=3)
            nseg = int(parts[0])
            label = parts[3].strip()
            pts = []
            for _ in range(nseg):
                vals = handle.readline().split()
                if len(vals) < 3:
                    continue
                pts.append(tuple(map(float, vals[:3])))
            loops.append((label, pts))
    return loops


def read_diag_segrog(path: Path) -> List[Tuple[str, List[Tuple[float, float, float]]]]:
    loops = []
    with path.open() as handle:
        total = int(handle.readline().split()[0])
        for _ in range(total):
            header = handle.readline()
            if not header:
                break
            parts = header.split(maxsplit=3)
            nseg = int(parts[0])
            label = parts[3].strip()
            pts = []
            for _ in range(nseg):
                vals = handle.readline().split()
                if len(vals) < 3:
                    continue
                pts.append(tuple(map(float, vals[:3])))
            loops.append((label, pts))
    return loops


def set_equal_3d(ax):
    xs = ax.get_xlim3d()
    ys = ax.get_ylim3d()
    zs = ax.get_zlim3d()
    ranges = [xs, ys, zs]
    centers = [0.5 * (r[0] + r[1]) for r in ranges]
    radius = 0.5 * max(r[1] - r[0] for r in ranges)
    ax.set_xlim3d(centers[0] - radius, centers[0] + radius)
    ax.set_ylim3d(centers[1] - radius, centers[1] + radius)
    ax.set_zlim3d(centers[2] - radius, centers[2] + radius)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as err:  # pragma: no cover
        print(f"run_xdiagno failed: {err}", file=sys.stderr)
        sys.exit(1)
