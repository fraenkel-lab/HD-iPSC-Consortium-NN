for l in ('HD109n1','HD28n6','HD60n5'):
	for k in range(1,6):
		full = open('../results/kmeans_cluster/%s/k5_%i.txt'%(l,k),'r')
		ds  = open('../results/kmeans_cluster_downsample/%s/k5_%i.txt'%(l,k),'r')

		out = open('../results/kmeans_cluster_downsample/%s/diff.txt'%l,'a')
		out.write('k%i\n'%k)

		fullgenes = [line.strip() for line in full]
		out.write('Out of %i genes in full sample\n'%len(fullgenes))
		overlap = set()
		for line in ds:
			if line.strip() in fullgenes:
				overlap.add(line.strip())
		out.write('%i genes overlap\n'%len(overlap))
