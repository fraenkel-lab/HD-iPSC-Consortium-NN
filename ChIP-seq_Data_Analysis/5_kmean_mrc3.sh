#!/bin/bash
n="HD109n1,HD109n4,HD21n1,HD28n6,HD33n1,HD60n5,HD60n8"

tss=2000
tags=40

outdir=../results/kmeans_cluster
mkdir $outdir
k=5
ma=3000
CMD="python kmeans_multi_inputs.py --ma $ma --k $k --outfile $outdir --TSS $tss --tags $tags $n"
echo $CMD
$CMD
