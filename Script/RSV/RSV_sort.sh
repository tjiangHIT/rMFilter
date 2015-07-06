#!/bin/bash
# sort -n -k 2 chr1 > chr1_sort.RSV
chr=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
N=24

for ((i=0; i<$N; i++))
do
    sort -n -k 2 chr${chr[$i]} > chr${chr[$i]}_sort.RSV
done

