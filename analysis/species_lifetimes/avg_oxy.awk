BEGIN {
    FS=","
}

NR == 1 {
    # store header names
    for (i=1; i<=NF; i++)
        name[i] = $i
    next
}

{
    # accumulate sums
    for (i=1; i<=NF; i++)
        sum[i] += $i
    row_count++
}

END {
    for (i=1; i<=NF; i++)
        printf "%s\t%f\n", name[i], sum[i]/row_count
}

