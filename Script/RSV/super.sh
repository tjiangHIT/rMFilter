#!/bin/bash

READ_DIR=/home/tjiang/2015_exam/results
RSV_DIR=/home/tjiang/2015_exam/results
PRO_DIR=/home/tjiang/2015_exam/results/RSV

# for RSV_msg
# 1. awk ins/del/dup/inv.csv -> RSV.ins/del/inv
# 2. cat RSV.ins/del/inv     -> RSV.msg
# 3. split RSV.msg           -> chr1.RSV/chr2.RSV ...
# 4. sort  chr1.RSV/chr2.RSV -> chr1_sort.RSV/chr2_sort.RSV

del=deletions.csv  
ins=insertions.csv  
inv=inversions.csv  
dup=tandemDuplications.csv

cd $RSV_DIR
# 1
awk -f $PRO_DIR/RSV_ins.awk $ins > RSV.ins
awk -f $PRO_DIR/RSV_del.awk $del > RSV.del
awk -f $PRO_DIR/RSV_dup.awk $dup > RSV.dup
awk -f $PRO_DIR/RSV_inv.awk $inv > RSV.inv
# 2
cat RSV.ins > RSV.msg
cat RSV.del >> RSV.msg
cat RSV.dup >> RSV.msg
cat RSV.inv >> RSV.msg
# 3
awk -f $PRO_DIR/split.awk RSV.msg
# 4
bash $PRO_DIR/RSV_sort.sh

