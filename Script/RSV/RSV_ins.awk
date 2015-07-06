#!/bin/awk
BEGIN {
    OFS="\t"
}
{
    if (NR>1) {
        print $5,$6-1,$6,$8,"INS"
    }
}
