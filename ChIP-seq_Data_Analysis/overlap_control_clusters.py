#! /usr/bin/python
# coding: utf-8

for sample in ('HD21n1','HD28n6','HD33n1'):
	clusters = open('../results/kmeans_cluster/%s/k5_1.txt'%sample, 'r')
	dn = open('../data/DE_down_genes.txt','r')
	up = open('../data/DE_up_genes.txt', 'r')
	cluster = []
	for line in clusters:
		cluster.append(line)
	dboth = []
	for line in dn:
		if line in cluster:
			dboth.append(line)
	upboth = []
	for line in up:
		if line in cluster:
			upboth.append(line)
	with open('../results/kmeans_cluster/%s/k5_1_downregulated.txt'%sample, 'w') as doverlap:
		for gene in dboth:
			doverlap.write(gene)
	with open('../results/kmeans_cluster/%s/k5_1_upregulated.txt'%sample, 'w') as uoverlap:
		for gene in upboth:
			uoverlap.write(gene)

clusters1 = '../results/kmeans_cluster/HD21n1'
clusters2 = '../results/kmeans_cluster/HD28n6'
clusters3 = '../results/kmeans_cluster/HD33n1'

for i in range(1,6):
	cluster1 = set(line.strip() for line in open('%s/k5_%i.txt'%(clusters1,i), 'r'))
	cluster2 = set(line.strip() for line in open('%s/k5_%i.txt'%(clusters2,i), 'r'))
	cluster3 = set(line.strip() for line in open('%s/k5_%i.txt'%(clusters3,i), 'r'))
	overlap = set.intersection(cluster1, cluster2,cluster3)
	
	with open('../results/kmeans_cluster/control_overlap_k5_%i.txt'%i, 'w') as overlap_file:
		for gene in overlap:
			overlap_file.write(gene+'\n')

