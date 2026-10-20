#!/usr/bin/env bash

input="$1"
output="$2"

if [[ -z "$input" || -z "$output" ]]; then
    echo "Usage: $0 input_log output_file"
    exit 1
fi

awk '
BEGIN {
    in_table = 0
    header_printed = 0
}

/^Step[[:space:]]+Time[[:space:]]+Temp[[:space:]]+PotEng[[:space:]]+TotEng[[:space:]]+Press[[:space:]]+Volume[[:space:]]+Density[[:space:]]+Lx[[:space:]]+Ly[[:space:]]+Lz[[:space:]]+Enthalpy/ {
    in_table = 1
    if (header_printed == 0) {
        print
        header_printed = 1
    }
    next
}

/^Loop/ {
    in_table = 0
    next
}

/^### NNP EW SUMMARY ###/ {
    next
}

in_table == 1 {
    if ($1 ~ /^[0-9]+$/) {
        print
    }
}
' "$input" > "$output"
