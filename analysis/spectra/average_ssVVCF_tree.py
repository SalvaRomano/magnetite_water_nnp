#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import numpy as np


def read_dat(path: Path):
    """
    Reads a .dat file with:
      - first two lines are header/comment lines (kept verbatim)
      - remaining lines are whitespace-separated numeric columns
    Returns: header_lines(list[str]), data(np.ndarray shape [N, M])
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header1 = f.readline()
        header2 = f.readline()
        if not header1 or not header2:
            raise ValueError(f"{path}: file is too short (needs 2 header lines).")
        try:
            data = np.loadtxt(f, dtype=float)
        except Exception as e:
            raise ValueError(f"{path}: failed to parse numeric data: {e}") from e

    if data.ndim == 1:
        data = data[None, :]

    return [header1.rstrip("\n"), header2.rstrip("\n")], data


def find_inputs(root: Path, rel_pattern: str):
    """
    Find files using a NON-recursive pattern relative to root.
    IMPORTANT: rel_pattern like 'trj_*/chunk_*/ssVVCF_wcf.dat' will not match deeper paths.
    """
    return sorted(root.glob(rel_pattern))


def main():
    p = argparse.ArgumentParser(
        description="Average columns of ssVVCF_wcf.dat across trj_*/chunk_* while excluding deeper subfolders."
    )
    p.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Root directory containing trj_* folders (default: current directory).",
    )
    p.add_argument(
        "--pattern",
        default="trj_*/chunk_*/ssVVCF_wcf.dat",
        help="Relative glob pattern from root. Default matches trj_*/chunk_*/ssVVCF_wcf.dat (no deeper paths).",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("ssVVCF_wcf_AVG.dat"),
        help="Output file path (default: ssVVCF_wcf_AVG.dat in current directory).",
    )
    p.add_argument(
        "--tol",
        type=float,
        default=1e-10,
        help="Tolerance for time-column equality check (default: 1e-10).",
    )
    args = p.parse_args()

    root = args.root.resolve()
    out = args.out.resolve()

    files = [f for f in find_inputs(root, args.pattern) if f.is_file()]

    if not files:
        print(
            f"ERROR: no files found under {root} matching pattern '{args.pattern}'",
            file=sys.stderr,
        )
        sys.exit(2)

    # Read first file as reference
    ref_header, ref_data = read_dat(files[0])
    ref_shape = ref_data.shape
    ref_time = ref_data[:, 0].copy()

    sum_data = ref_data.copy()
    n = 1

    for f in files[1:]:
        header, data = read_dat(f)

        if data.shape != ref_shape:
            raise ValueError(
                f"Shape mismatch:\n  {files[0]}: {ref_shape}\n  {f}: {data.shape}"
            )

        if not np.allclose(data[:, 0], ref_time, rtol=0.0, atol=args.tol):
            idx = np.argmax(np.abs(data[:, 0] - ref_time) > args.tol)
            raise ValueError(
                f"Time column mismatch beyond tol={args.tol} at row {idx}:\n"
                f"  {files[0]} time={ref_time[idx]}\n"
                f"  {f} time={data[idx,0]}"
            )

        sum_data += data
        n += 1

    avg_data = sum_data / float(n)
    avg_data[:, 0] = ref_time  # keep exact reference time values

    with out.open("w", encoding="utf-8") as w:
        w.write(ref_header[0] + "\n")
        w.write(ref_header[1] + "\n")
        np.savetxt(w, avg_data, fmt="%.16e")

    print(f"Wrote averaged file from {n} inputs:\n  {out}")
    print("Matched inputs (excluded deeper subfolders automatically):")
    for f in files:
        print(f"  {f.relative_to(root)}")


if __name__ == "__main__":
    main()

