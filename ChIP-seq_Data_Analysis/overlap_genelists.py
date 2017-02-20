deg = '../data/DE_all_genes.txt'
degUp = '../data/DE_up_genes.txt'
degDown = '../data/DE_down_genes.txt'
chipHD = '../results/MAnorm/HDvsControl/output_tables/HD_genes.txt'
chipCtrl = '../results/MAnorm/HDvsControl/output_tables/Ctrl_genes.txt' 

out = open('../results/MAnorm/HDvsControl/output_tables/DE_overlap.txt', 'w')
uout = open('../results/MAnorm/HDvsControl/output_tables/up_DE_overlap.txt', 'w')
dout = open('../results/MAnorm/HDvsControl/output_tables/down_DE_overlap.txt', 'w')

with open(deg) as degfile:
	deggenes = [line.strip() for line in degfile]
	with open(chipHD) as chipHDfile:
		for line in chipHDfile:
			if line.strip() in deggenes:
				out.write(line)
	with open(chipCtrl) as chipCtrlfile:
		for line in chipCtrlfile:
			if line.strip() in deggenes:
				out.write(line)

with open(degUp) as degfile:
	deggenes = [line.strip() for line in degfile]
	with open(chipHD) as chipHDfile:
		for line in chipHDfile:
			if line.strip() in deggenes:
				uout.write(line)

with open(degDown) as degfile:
	deggenes = [line.strip() for line in degfile]
	with open(chipCtrl) as chipCtrlfile:
		for line in chipCtrlfile:
			if line.strip() in deggenes:
				dout.write(line)
