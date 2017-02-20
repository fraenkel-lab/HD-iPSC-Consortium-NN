#!/bin/bash
script=overlap_profiles.py
#in="../results/H3K27ac/HD28n6/k5_"
#in="../results/kmeans_cluster/HD28n6/k5_"
in="../results/kmeans_cluster/control_overlap_k5_"
nc=5

CMD="python overlap_control_clusters.py"
echo $CMD
$CMD

yl="#_genes_w/_lower_expression_in_HD"
downHD=../data/DE_down_genes.txt

all=../data/msns_all_gene_sym.txt
#out="../results/kmeans_cluster/HD28n6/k5_downHD_all_exp_iPSC_fisher.pdf"
out="../results/kmeans_cluster/k5_downHD_all_exp_iPSC_fisher.pdf"

CMD="python $script --ylabel $yl --outfile $out --k $nc $downHD $in $all"
echo $CMD
$CMD

yl="#_genes_w/_higher_expression_in_HD"
upHD=../data/DE_up_genes.txt
all=../data/msns_all_gene_sym.txt
#out="../results/kmeans_cluster/HD28n6/k5_upHD_all_exp_iPSC_fisher.pdf"
out="../results/kmeans_cluster/k5_upHD_all_exp_iPSC_fisher.pdf"

CMD="python $script --ylabel $yl --outfile $out --k $nc $upHD $in $all"
echo $CMD
$CMD

