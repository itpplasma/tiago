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
import struct
import subprocess
import sys
import zlib
from collections import OrderedDict
from pathlib import Path
from typing import List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]

SKIP_EXIT_CODE = 125


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
        "--tolerance",
        type=float,
        default=1e-3,
        help="Absolute tolerance for flux/voltage comparisons",
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


def run_tiago(args: argparse.Namespace, out_dir: Path) -> Tuple[Path, Path]:
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
    ]
    subprocess.run(cmd, check=True)
    return flux_out, segrog_out


def run_xdiagno(args: argparse.Namespace, out_dir: Path) -> Tuple[Path, Path]:
    binary = resolve_xdiagno()
    if not binary:
        print(
            "TIAGO_XDIAGNO not set and xdiagno missing from PATH; "
            "skipping cross-code comparison",
            file=sys.stderr,
        )
        sys.exit(SKIP_EXIT_CODE)
    flux_path = write_diagno_flux(out_dir, args.flux)
    seg_path = write_diagno_segrog(out_dir, args.segrog, float(args.seg_area))
    write_diagno_control(out_dir, flux_path, seg_path)
    write_vmec_input(out_dir)
    coil_path = str(write_diagno_coils(out_dir, args.coil))
    cmd = [binary, "-vac", "-coil", coil_path, "-noverb"]
    if args.extra_xdiagno_args:
        cmd.extend(args.extra_xdiagno_args)
    try:
        subprocess.run(cmd, check=True, cwd=str(out_dir))
    except FileNotFoundError as err:
        print(f"failed to execute {binary}: {err}", file=sys.stderr)
        sys.exit(SKIP_EXIT_CODE)
    except subprocess.CalledProcessError as err:
        print(f"xdiagno exited with {err.returncode}", file=sys.stderr)
        raise
    flux_raw = out_dir / "diagno_flux."
    seg_raw = out_dir / "diagno_seg."
    flux_csv = Path(args.legacy_flux) if args.legacy_flux else out_dir / "diagno_flux.csv"
    seg_csv = Path(args.legacy_segrog) if args.legacy_segrog else out_dir / "diagno_segrog.csv"
    convert_diagno_output(flux_raw, flux_csv, expect_index=True)
    convert_diagno_output(seg_raw, seg_csv, expect_index=False)
    return flux_csv, seg_csv


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
    tiago_path: Path, legacy_path: Path, tolerance: float, report: Path
) -> Tuple[int, List[Tuple[str, float, float]]]:
    tiago = load_csv(tiago_path)
    legacy = load_csv(legacy_path)
    failures = 0
    pairs: List[Tuple[str, float, float]] = []
    with report.open("w") as handle:
        handle.write("label,tiago,legacy,diff\n")
        for label, value in tiago.items():
            legacy_value = legacy.get(label)
            pairs.append((label, value, legacy_value))
            if label not in legacy:
                failures += 1
                handle.write(f"{label},{value},MISSING,NaN\n")
                continue
            diff = abs(value - legacy_value)
            if diff > tolerance:
                failures += 1
                handle.write(f"{label},{value},{legacy_value},{diff}\n")
        extra = set(legacy.keys()) - set(tiago.keys())
        for label in extra:
            failures += 1
            handle.write(f"{label},MISSING,{legacy[label]},NaN\n")
            pairs.append((label, math.nan, legacy[label]))
    return failures, pairs


def main() -> int:
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    tiago_flux, tiago_seg = run_tiago(args, out_dir)
    legacy_flux, legacy_seg = run_xdiagno(args, out_dir)
    flux_report = out_dir / "flux_diff.csv"
    seg_report = out_dir / "segrog_diff.csv"
    failures = 0
    flux_failures, flux_pairs = compare_arrays(
        tiago_flux, legacy_flux, args.tolerance, flux_report
    )
    seg_failures, seg_pairs = compare_arrays(
        tiago_seg, legacy_seg, args.tolerance, seg_report
    )
    failures += flux_failures + seg_failures
    generate_plot(flux_pairs, out_dir / "flux_plot.png")
    generate_plot(seg_pairs, out_dir / "segrog_plot.png")
    if failures > 0:
        print(
            f"Detected {failures} mismatched diagnostics; see {out_dir}",
            file=sys.stderr,
        )
        return 2
    print("Tiago and xdiagno outputs agree within tolerance")
    return 0


def write_diagno_control(out_dir: Path, flux_path: Path, seg_path: Path) -> Path:
    control_path = out_dir / "diagno.control"
    payload = (
        "&diagno_in\n"
        f"  flux_diag_file = '{flux_path}',\n"
        f"  seg_rog_file = '{seg_path}',\n"
        "  nu = 64,\n"
        "  nv = 64,\n"
        "  int_type = 'midpoint',\n"
        "  int_step = 1,\n"
        "  lrphiz = .false.,\n"
        "  lvc_field = .false.,\n"
        "  luse_extcur = .true.,\n"
        "  units = 1.0,\n"
        "/\n"
    )
    control_path.write_text(payload)
    return control_path


def write_vmec_input(out_dir: Path) -> None:
    source = PROJECT_ROOT / "tests" / "data" / "input.diagno.stub"
    target = out_dir / "input."
    shutil.copyfile(source, target)


def write_diagno_coils(out_dir: Path, source_path: str) -> Path:
    dest = out_dir / "coils.tiago"
    src = Path(source_path).resolve()
    with src.open() as handle:
        header = handle.readline().strip()
        try:
            count = int(header.split()[0])
        except ValueError as exc:
            raise ValueError(f"invalid Tiago coil header '{header}'") from exc
        points = []
        for _ in range(count):
            line = handle.readline()
            if not line:
                break
            values = line.split()
            if len(values) < 4:
                continue
            x, y, z, current = map(float, values[:4])
            points.append((x, y, z, current))
    with dest.open("w") as handle:
        handle.write("periods 1\n")
        handle.write("begin filament\n")
        handle.write("mirror NIL\n")
        for x, y, z, current in points:
            handle.write(
                f" {x: .15E} {y: .15E} {z: .15E} {current: .15E}\n"
            )
        handle.write("end\n")
    return dest


def write_diagno_segrog(out_dir: Path, source_path: str, area_value: float) -> Path:
    dest = out_dir / "segrog.diagno"
    src = Path(source_path).resolve()
    with src.open() as handle:
        lines = [line.rstrip("\n") for line in handle]
    if not lines:
        raise ValueError("segrog file is empty")
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
    return dest


def write_diagno_flux(out_dir: Path, source_path: str) -> Path:
    dest = out_dir / "fluxloop.diagno"
    src = Path(source_path).resolve()
    with src.open() as handle:
        lines = [line.rstrip("\n") for line in handle]
    if not lines:
        raise ValueError("flux file is empty")
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
            for _ in range(nseg):
                if idx >= len(lines):
                    raise ValueError("unexpected end of flux file")
                coords = lines[idx].split()
                idx += 1
                if len(coords) < 3:
                    continue
                x, y, z = map(float, coords[:3])
                out.write(f" {x: .10E} {y: .10E} {z: .10E}\n")
    return dest


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


def generate_plot(pairs: List[Tuple[str, float, float]], path: Path) -> None:
    if not pairs:
        return
    width = max(480, len(pairs) * 40)
    height = 320
    data = bytearray([255] * width * height * 3)
    values = [v for (_, v, lv) in pairs for v in (v, lv) if v is not None and not math.isnan(v)]
    if not values:
        return
    min_val = min(values)
    max_val = max(values)
    if math.isclose(max_val, min_val):
        max_val += 1.0
        min_val -= 1.0
    margin = 40
    def to_point(idx: int, value: float) -> Tuple[int, int]:
        if len(pairs) == 1:
            pos = 0
        else:
            pos = idx / (len(pairs) - 1)
        x = margin + int(pos * (width - 2 * margin))
        span = height - 2 * margin
        y = height - margin - int(((value - min_val) / (max_val - min_val)) * span)
        return x, max(0, min(height - 1, y))

    def draw_line(points: List[Tuple[int, int]], color: Tuple[int, int, int]) -> None:
        if not points:
            return
        for px, py in points:
            stamp(px, py, color)
        for (x0, y0), (x1, y1) in zip(points, points[1:]):
            steps = max(abs(x1 - x0), abs(y1 - y0))
            if steps == 0:
                continue
            for s in range(steps + 1):
                t = s / steps
                xs = int(round(x0 + (x1 - x0) * t))
                ys = int(round(y0 + (y1 - y0) * t))
                stamp(xs, ys, color)

    def stamp(x: int, y: int, color: Tuple[int, int, int]) -> None:
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                xx = x + dx
                yy = y + dy
                if 0 <= xx < width and 0 <= yy < height:
                    offset = (yy * width + xx) * 3
                    data[offset : offset + 3] = bytes(color)

    tiago_points = []
    legacy_points = []
    for idx, (_, tiago_val, legacy_val) in enumerate(pairs):
        if tiago_val is not None and not math.isnan(tiago_val):
            tiago_points.append(to_point(idx, tiago_val))
        if legacy_val is not None and not math.isnan(legacy_val):
            legacy_points.append(to_point(idx, legacy_val))

    draw_axes(data, width, height, margin)
    draw_line(legacy_points, (246, 141, 64))
    draw_line(tiago_points, (24, 181, 170))
    write_png(path, width, height, data)


def draw_axes(buffer: bytearray, width: int, height: int, margin: int) -> None:
    color = (200, 200, 200)
    for x in range(margin, width - margin):
        idx = ((height - margin) * width + x) * 3
        buffer[idx : idx + 3] = bytes(color)
    for y in range(margin, height - margin):
        idx = (y * width + margin) * 3
        buffer[idx : idx + 3] = bytes(color)


def write_png(path: Path, width: int, height: int, data: bytearray) -> None:
    def chunk(tag: bytes, payload: bytes) -> bytes:
        crc = zlib.crc32(tag + payload) & 0xFFFFFFFF
        return struct.pack(
            ">I", len(payload)
        ) + tag + payload + struct.pack(
            ">I", crc
        )

    raw_rows = []
    row_bytes = width * 3
    for y in range(height):
        start = y * row_bytes
        raw_rows.append(b"\x00" + bytes(data[start : start + row_bytes]))
    raw = b"".join(raw_rows)
    png = b"\x89PNG\r\n\x1a\n"
    png += chunk(
        b"IHDR",
        struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0),
    )
    png += chunk(b"IDAT", zlib.compress(raw, 9))
    png += chunk(b"IEND", b"")
    with path.open("wb") as handle:
        handle.write(png)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as err:  # pragma: no cover
        print(f"run_xdiagno failed: {err}", file=sys.stderr)
        sys.exit(1)
