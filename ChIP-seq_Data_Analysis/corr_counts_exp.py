#Amanda Daigle
#4/25/14

import math
from scipy.stats.stats import pearsonr
from optparse import OptionParser
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype']=42

def correlate(names, countsdir, expfilename, XYfilename):
    #load expression data
    expfile = open(expfilename, 'rb')
    e = []
    for line in expfile:
        e.append(line.split())
    for c in range(1,len(e[0])):
        e[0][c] = 'HD'+ e[0][c][:-3]
    expfile.close()	
    
    #load genes to ignore due to XY chrom
    XYfile = open(XYfilename, 'r')
    XY = set(line.strip() for line in XYfile.readlines())
    XYfile.close()

    for n in names.split(','):
        #create dictionary of {gene:counts}
        countsfile = open(countsdir+'%scounts.txt'%n, 'rb')
        countslist=countsfile.readlines()
        genesfile = open(countsdir+'%sgenes.txt'%n, 'rb')
        genes=genesfile.readlines()
        counts = {}
        for g in range(0,len(genes)):
            count = float(countslist[g].strip())
            counts[genes[g].strip()] = count
        
        #create dictionary of {gene:RPKM}
        exps={}
        namecol = e[0].index(n)
        for line in e[1:]:
            gene = line[0]
            exps[gene] = float(line[namecol])
 
        #correlate
        counts_overlap = []
        exp_overlap = []
        maxc =0.0
        maxc_exp = 0.0
        max_gene = ''
        for gene in exps.keys():
            if gene=='' or gene=='n/a': continue
            if gene in XY: continue
            #log2 exp with pseudocount +1
            exp = math.log(exps[gene]+1,2)
            #exp = exps[gene]
            try:
                count = math.log((counts[gene]+1),2)
                #count = counts[gene]
            except KeyError:
                #count = 0.0
                continue
            counts_overlap.append(count)
            exp_overlap.append(exp)
            if count > maxc:
                maxc = count
                maxc_exp = exp
                max_gene = gene
        genesfile.close()
        countsfile.close()
        print "The gene with the maximum counts for %s is %s" %(n,max_gene), maxc, maxc_exp 
        (r,p) = pearsonr(counts_overlap, exp_overlap)
        print '%i genes total.\n'%len(counts_overlap)
        print 'p = %f'%p
        plt.figure()
        plt.scatter(counts_overlap, exp_overlap, edgecolor='none', alpha=0.3)
        #plt.figtext(0.6, 0.8, 'r=%.4f, p=%.4f'%(r,p))
        #from observation I know all p values are <1e-6
        plt.figtext(0.6,0.8, 'r=%.4f, p<1e-6'%r)
        plt.xlabel('Log_2 (H3K4me3 counts)',fontsize='x-large')
        plt.ylabel('Log_2 (Expression) (RPKM)',fontsize='x-large')
        plt.savefig('%s/%s_counts_vs_exp.pdf'%(countsdir,n))

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--counts", dest="countsdir", default='../results/tss_counts_t2000_mrc3_all/')
    parser.add_option("--exp", dest="expfile", default='../data/msnsRPKM.txt')
    parser.add_option("--names", dest="names", default='HD109n1,HD109n4,HD21n1,HD28n6,HD33n1,HD60n5,HD60n8')
    parser.add_option("--XYgenes", dest="XYfile", default='../data/genes_on_chrXY.txt')
    
    (opts, args) = parser.parse_args()
    
    correlate(opts.names, opts.countsdir, opts.expfile, opts.XYfile)
    
if __name__=='__main__': main()
