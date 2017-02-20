#!/bin/bash
script="compile_counts_matrix_reorg.py"
cdir=../results/tss_counts_t2000_mrc3_all/
g_in=$cdir/HD109n1genes.txt

Cfi=$cdir/HD109n1counts_rpmlog2.txt
for x in HD109n4  HD21n1  HD28n6  HD33n1  HD60n5  HD60n8
do
    Cfi="$Cfi,$cdir/${x}counts_rpmlog2.txt"
done

Cmat="$cdir/matrix/allreps.txt"
mkdir $cdir/matrix/

CMD="python $script --outfile $Cmat --genes_in $g_in --genes_out $g_in $Cfi"
echo $CMD
$CMD

script="create_pycluster_matrix.py"
a=HD109n1,HD109n4,HD21n1,HD28n6,HD33n1,HD60n5,HD60n8
Cmat="$cdir/matrix/allreps.txt"
b="$cdir/HD109n1genes.txt"
out="$cdir/matrix/allreps_pycluster.txt"

CMD="python $script --outfile $out $Cmat $a $b"
echo $CMD
$CMD
