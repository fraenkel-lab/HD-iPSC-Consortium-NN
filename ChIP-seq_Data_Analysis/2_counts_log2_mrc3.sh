#!/bin/bash
script=plt_array_dist.py
cdir=../results/tss_counts_t2000_mrc3_all/

for x in HD109n1  HD109n4 HD21n1  HD28n6  HD33n1  HD60n5  HD60n8
do
    COUNTS=$cdir/${x}counts.txt
    readcounts=../data/mapped_read_counts.txt
    OUT=$cdir/${x}counts_rpmlog2.pdf
    CMD="wqsub.py python $script --outfile $OUT --log --sample $x --counts ${readcounts} $COUNTS"
    echo $CMD
    $CMD
done

