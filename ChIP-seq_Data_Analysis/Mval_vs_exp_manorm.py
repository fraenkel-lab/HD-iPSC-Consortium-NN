#Amanda Daigle 
#7/31/14

import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype']=42
import math
from scipy.stats.stats import pearsonr

mark = 'H3K4me3'
outfile = '../results/MAnorm/HDvsControl/output_figures/Mval_vs_exp.pdf'

treat_pks = open('../results/MAnorm/HDvsControl/output_tables/biased_peaks_of_HD_all_peaks_genes.txt', 'r')
cont_pks = open('../results/MAnorm/HDvsControl/output_tables/biased_peaks_of_Ctrl_all_peaks_genes.txt', 'r')
unbiased_pks = open('../results/MAnorm/HDvsControl/output_tables/unbiased_peaks_all_peaks_genes.txt', 'r')
exp_file = open('../data/msnsRPKM_FC.txt', 'r')

#Store m-value and gene symbol for all peaks. pks = {geneSymbol: {'exp':0.0, 'm':mValue, 'exp_given':0}}
pks={}
ms_for_this_peak=[]
dists_for_this_peak=[]
chrLast = ''
startLast = None
for f in (unbiased_pks,treat_pks,cont_pks):
    for line in f:
        if line.startswith('k'): continue
        words = line.strip().split('\t')
        #account for duplicate entries - if there are multiple dists for one peak, use the closest
        chrom = words[2]
        start = int(words[3])
        dist = int(words[15])
        m = float(words[6])
        gene = words[1]
        #first peak
        if startLast == None:
            chrLast = chrom
            startLast = start
            dists_for_this_peak=[dist]
            ms_for_this_peak=[m]
            continue
        #new chrom
        if chrom != chrLast:
            if chrom == 'chrY' or chrom == 'chrX': 
                continue
            chrLast = chrom
            startLast = start
            m = ms_for_this_peak[dists_for_this_peak.index(min(dists_for_this_peak))]
            pks[gene] = {'exp':0.0, 'm':m, 'exp_given':None}
            dists_for_this_peak=[dist]
            ms_for_this_peak=[m]
        else:
            #same peak
            if start == startLast:
                dists_for_this_peak.append(dist)
                ms_for_this_peak.append(m)
            #new peak
            else:
                startLast = start
                m = ms_for_this_peak[dists_for_this_peak.index(min(dists_for_this_peak))]
                pks[gene] = {'exp':0.0, 'm':m, 'exp_given':None}
                dists_for_this_peak=[dist]
                ms_for_this_peak=[m]
    #last peak
    m = ms_for_this_peak[dists_for_this_peak.index(min(dists_for_this_peak))]
    pks[gene] = {'exp':0.0, 'm':m, 'exp_given':None}

#Update exp value for all peaks
exp_file.readline()
for line in exp_file:
    words = line.strip().split('\t')
    gene = words[0]
    if gene in pks:
        exp = float(words[12])
        try:
            logexp = math.log(exp,2)
        except ValueError:
            #when fold change is 0
            pks[gene]['exp_given'] = 0
        pks[gene]['exp'] = logexp
        if logexp >= 0:
            pks[gene]['exp_given'] = 1
        else:
            pks[gene]['exp_given'] = -1


upreg_m = [pks[gene]['m'] for gene in pks.keys() if pks[gene]['exp_given'] !=None]
upreg_exp = [pks[gene]['exp'] for gene in pks.keys() if pks[gene]['exp_given'] !=None]
#downreg_m = [pks[gene]['m'] for gene in pks.keys() if pks[gene]['exp_given'] == -1]
#downreg_exp = [pks[gene]['exp'] for gene in pks.keys() if pks[gene]['exp_given'] == -1]
#noreg_m = [pks[gene]['m'] for gene in pks.keys() if pks[gene]['exp_given'] ==0]
#noreg_exp = [pks[gene]['exp'] for gene in pks.keys() if pks[gene]['exp_given'] ==0]

(r,p) = pearsonr(upreg_m, upreg_exp)
print '%i genes total' %len(upreg_m)
print 'r=%f, p=%f'%(r,p)

#from MAnorm.py
fig=plt.figure()
#plt.scatter(noreg_m, noreg_exp, c='0.5',label='No expression change')
plt.scatter(upreg_m, upreg_exp, c='b', edgecolor='none', alpha=0.3) #, label='Upregulated')
#plt.scatter(downreg_m, downreg_exp, c='b', label='Downregulated')
plt.figtext(0.4,0.8,'r=%.4f, p<1e-6'%(r))

#plt.grid(axis='y')
plt.xlabel('Normalized logFC in # of reads under peak (M value)')
plt.ylabel('Log fold change in expression (RPKM)')
plt.title('Expression Change vs. Differential %s Peaks'%(mark))
#plt.legend(loc=2)

fig.savefig(outfile)
