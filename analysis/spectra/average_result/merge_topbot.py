#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import numpy as np


def parse_second_header_line(line: str):
    """
    Parse a header like:
     #time, full , oct-H2O-top, oct-H2O-bot, ... , rest
    Returns list of column names, stripped.
    """
    s = line.strip()
    if s.startswith("#"):
        s = s[1:].strip()
    # split by commas, strip whitespace
    cols = [c.strip() for c in s.split(",")]
    # drop empty tokens if any
    cols = [c for c in cols if c]
    return cols


def format_header_line(cols):
    # keep a simple, readable comma-separated header
    return " #"+", ".join(cols)


def main():
    ap = argparse.ArgumentParser(description="Average *-top and *-bot columns into a single column.")
    ap.add_argument("inp", type=Path, help="Input .dat file (already averaged across runs).")
    ap.add_argument("-o", "--out", type=Path, default=None, help="Output file (default: <inp_stem>_topbot.dat)")
    args = ap.parse_args()

    inp = args.inp.resolve()
    out = args.out.resolve() if args.out else inp.with_name(inp.stem + "_topbot.dat")

    # Read two header lines, then numeric table
    with inp.open("r", encoding="utf-8", errors="replace") as f:
        h1 = f.readline()
        h2 = f.readline()
        if not h1 or not h2:
            raise SystemExit(f"ERROR: {inp} is too short (needs 2 header lines).")
        data = np.loadtxt(f, dtype=float)

    if data.ndim == 1:
        data = data[None, :]

    colnames = parse_second_header_line(h2)
    ncols = data.shape[1]
    if len(colnames) != ncols:
        raise SystemExit(
            f"ERROR: header columns ({len(colnames)}) != data columns ({ncols}).\n"
            f"Header: {colnames}"
        )

    # Build mapping: base -> (top_idx, bot_idx)
    top_idx = {}
    bot_idx = {}
    for i, name in enumerate(colnames):
        if name.endswith("-top"):
            top_idx[name[:-4]] = i
        elif name.endswith("-bot"):
            bot_idx[name[:-4]] = i

    bases = sorted(set(top_idx) & set(bot_idx), key=lambda b: min(top_idx[b], bot_idx[b]))

    # Create new columns in original order:
    # - When we encounter a *-top whose base has a matching bot, output the averaged base once
    # - Skip the corresponding *-bot
    new_names = []
    new_cols = []

    i = 0
    while i < len(colnames):
        name = colnames[i]
        if name.endswith("-top"):
            base = name[:-4]
            if base in top_idx and base in bot_idx:
                ti = top_idx[base]
                bi = bot_idx[base]
                # only emit when we reach the top column index (to preserve order)
                if i == ti:
                    new_names.append(base)
                    new_cols.append(0.5 * (data[:, ti] + data[:, bi]))
                i += 1
                continue

        if name.endswith("-bot"):
            base = name[:-4]
            if base in top_idx and base in bot_idx:
                # skip the bot column entirely
                i += 1
                continue

        # default: keep column
        new_names.append(name)
        new_cols.append(data[:, i])
        i += 1

    new_data = np.column_stack(new_cols)

    # Write output (keep first header line verbatim; regenerate the second)
    

    # Update first header line to reflect the new number of columns (excluding the leading "#")
    new_ncols = len(new_names)

    # If the first line is like "#    ,  1   ,  2 , ...", regenerate it safely:
    h1_new = "#    , " + ", ".join(f"{i:>3d}" for i in range(1, new_ncols))  # 1..(N-1)

    with out.open("w", encoding="utf-8") as w:
        w.write(h1_new + "\n")
        w.write(format_header_line(new_names) + "\n")
        np.savetxt(w, new_data, fmt="%.16e")



    print(f"Wrote: {out}")
    if bases:
        print("Merged pairs:")
        for b in bases:
            print(f"  {b}-top + {b}-bot -> {b}")
    else:
        print("No *-top/*-bot pairs found to merge.")


if __name__ == "__main__":
    main()

