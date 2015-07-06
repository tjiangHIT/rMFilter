BEGIN {
    OFS="\t"
}
{
    if (NR>1) {
        print $2,$4,$4+1, $5*($6-1), "INS"
    }
}
