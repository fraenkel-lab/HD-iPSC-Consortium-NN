#!/bin/bash
script="compile_counts_matrix_reorg.py"
cdir=../results/tss_counts_t2000_mrc3_all
g_in=$cdir/HD109n1genes.txt

Cfi=$cdir/HD109n1counts_rpmlog2.txt
for x in HD109n4 HD21n1  HD28n6  HD33n1  HD60n5  HD60n8
do
    Cfi="$Cfi,$cdir/${x}counts_rpmlog2.txt"
done

Cmat="$cdir/matrix/allreps_DEgenes.txt"
mkdir $cdir/matrix/
g_out=../data/DE_all_genes.txt

CMD="python $script --outfile $Cmat --genes_in $g_in --genes_out $g_out $Cfi"
echo $CMD
$CMD

script="create_pycluster_matrix.py"
a=HD109n1,HD109n4,HD21n1,HD28n6,HD33n1,HD60n5,HD60n8
Cmat="$cdir/matrix/allreps_DEgenes.txt"
out="$cdir/matrix/allreps_DEgenes_pycluster.txt"

CMD="python $script --outfile $out $Cmat $a $g_out"
echo $CMD
$CMD
