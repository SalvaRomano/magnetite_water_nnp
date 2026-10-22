awk -v start=300000 -v end=1530000 '

# When we hit a TIMESTEP header line
/^ITEM: TIMESTEP/ {
    header = $0               # save "ITEM: TIMESTEP"
    getline                    # read the timestep number
    ts_line = $0
    ts = ts_line + 0           # convert to number

    # Check bounds: include start, exclude end
    if (ts >= start && ts < end) {
        keep = 1
        print header
        print ts_line
    } else {
        keep = 0
    }
    next
}

# For all other lines: print only if the block is to be kept
keep { print }

' trjmagnetite.dump> trjcut.dump
cat trjcut.dump  continue/trjmagnetite.dump  >  trjtot.dump
