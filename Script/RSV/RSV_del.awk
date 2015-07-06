BEGIN {
    OFS="\t"
}
{
    if (NR > 1) {
        print $2, $3-1, $4+1, $5, "DEL"
    }
}
