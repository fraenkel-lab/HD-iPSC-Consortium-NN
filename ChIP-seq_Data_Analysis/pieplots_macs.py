#Amanda Daigle 
#7/31/14

import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=14
import sys

mark = 'H3K4me3'
for n in ('HD', 'Ctrl'):    
    genefile = open('../results/MACS/%s/%s_no180s_mfold10,30_pval1e-5_genes.txt'%(n,n), 'r')
    peakfile = open('../results/MACS/%s/%s_no180s_mfold10,30_pval1e-5_peaks.bed'%(n,n), 'r')
    outfile = '../results/MACS/%s/%s_peaks_piechart.pdf'%(n,n)

    types = {'promoter':0, 'after':0, 'intron':0, 'exon': 0}
        
    #read mapped gene file, store closest map for each peak
    peaks={} #{chrom:{peakStart:[dist, type]}}
    for line in genefile:
        words = line.strip().split('\t')
        #ignore first line
        if words[0] == 'knownGeneID': continue
        chrom = words[2]
        start = int(words[3])-1
        dist = abs(int(words[12]))
        maptype = words[14]
        if maptype == 'gene':
            maptype = words[15]
        if chrom not in peaks:
            #new chrom
            peaks[chrom] = {start:[dist,maptype]}
        else:
            if start in peaks[chrom].keys():
                #account for duplicate entries - choose closest gene and store type
                if dist < peaks[chrom][start][0]:
                    #closer gene
                    peaks[chrom][start] = [dist, maptype]
            else: peaks[chrom][start] = [dist, maptype]

    #count types - 1 per peak in peak file
    types = {'promoter':0, 'after':0, 'intron':0, 'exon': 0, 'inter': 0}
    totalpks = 0
    for line in peakfile:
        words = line.strip().split('\t')
        chrom = words[0]
        start = int(words[1])
        if chrom in peaks:
            if start in peaks[chrom]:
                types[peaks[chrom][start][1]] += 1
            else:
                types['inter'] += 1
        else:
            types['inter'] += 1
        totalpks += 1
            
    # make a square figure and axes
    fig = plt.figure(figsize=(6,6))
    pie_ax = fig.add_axes((0.3,0.2,0.4,0.4))

    # The slices will be ordered and plotted counter-clockwise.
    labels = ['exon: %i'%types['exon'],'intron: %i'%types['intron'],'promoter: %i'%types['promoter'],'intergenic: %i'%types['inter'], 'after: %i'%types['after']]
    fracs = [types['exon'], types['intron'], types['promoter'], types['inter'], types['after']]
    plt.pie(fracs, labels=labels) #, autopct='%1.1f%%')
        
    plt.title('Total %s peaks in %s'%(mark,n))
    plt.figtext(.5,.1, 'Total: %i'%totalpks, ha='center')
    fig.savefig(outfile)
