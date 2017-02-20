#!/bin/bash

for d in c s
do
    rdist=e
    cdist=${d}
    outdir=../results/cluster_t2000_mrc3_all/
    mkdir $outdir
    mat="../results/tss_counts_t2000_mrc3_all/matrix/allreps_pycluster.txt"
    clus_out=$outdir/all_reps_counts_log2_row_${rdist}_col_${cdist}
    #clus_out=$outdir/all_reps_counts_log2_col_${cdist}

    CMD="python pycluster_matrix.py --rdist $rdist --cdist $cdist --clusrow --cluscol --scale --outfile $clus_out $mat"
    #CMD="python pycluster_matrix.py --rdist $rdist --cdist $cdist --cluscol --scale --outfile $clus_out $mat"
    echo $CMD
    $CMD
done
