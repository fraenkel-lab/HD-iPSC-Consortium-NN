#!/bin/bash

n="HD109n1,HD109n4,HD21n1,HD28n6,HD33n1,HD60n5,HD60n8"
clusterdir='../results/kmeans_cluster'
outdir='../results/kmeans_cluster'
k=5

CMD="python count_genes_per_cluster.py --names $n --clusterdir $clusterdir --outdir $outdir -k $k"
echo $CMD
$CMD
