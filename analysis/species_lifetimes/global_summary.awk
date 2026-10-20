awk '
{
    key=$1
    val=$2

    # Record the order only the first time we encounter the key
    if (!(key in seen)) {
        seen[key] = ++order
        keys[order] = key
    }

    count[key]++
    sum[key] += val
    sumsq[key] += val * val
}
END {
    for (i = 1; i <= order; i++) {
        key = keys[i]
        mean = sum[key] / count[key]
        sigma = sqrt(sumsq[key] / count[key] - mean*mean)
        printf "%s\t%f\t%f\n", key, mean, sigma
    }
}
' */avg*.txt > global_stats_2.txt

