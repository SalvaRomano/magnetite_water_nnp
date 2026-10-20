#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./split_xyz_chunks.sh big.xyz outdir 50000
# Defaults:
#   outdir=chunks, N=50000

infile="${1:?Need input .xyz file}"
outdir="${2:-chunks}"
N="${3:-50000}"

mkdir -p "$outdir"

awk -v OUTDIR="$outdir" -v N="$N" '
function open_new_chunk() {
  chunk++
  frame_in_chunk = 0
  subdir = sprintf("%s/chunk_%01d", OUTDIR, chunk)
  cmd = "mkdir -p " quote(subdir)
  system(cmd)
  outfile = sprintf("%s/pos.xyz", subdir, chunk)
}

function quote(s,    t) { # minimal shell quoting for mkdir path
  gsub(/'\''/, "'\''\"'\''\"'\''", s)
  return "'"'"'" s "'"'"'"
}

BEGIN {
  chunk = 0
  frame_in_chunk = 0
  open_new_chunk()
}

# Each frame starts with natoms line
{
  natoms = $1 + 0
  if (natoms <= 0) {
    print "ERROR: Expected natoms as positive integer at line " NR > "/dev/stderr"
    exit 1
  }

  # Start a new chunk if we already wrote N frames into the current one
  if (frame_in_chunk >= N) {
    open_new_chunk()
  }

  # Write natoms line + comment line
  print $0 >> outfile
  if ((getline line2) <= 0) { exit 0 }
  print line2 >> outfile

  # Write natoms coordinate lines
  for (i = 0; i < natoms; i++) {
    if ((getline linec) <= 0) { exit 0 }
    print linec >> outfile
  }

  frame_in_chunk++
}
' "$infile"

