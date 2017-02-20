#Amanda Kedaigle
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
        countsoverlaps = []
        expoverlaps = []
        for k in range(1,6):
            classfile = open('../results/kmeans_cluster/%s/k5_%i.txt'%(n,k), 'r')
            classgenes = set(line.strip() for line in classfile.readlines())
            classfile.close()
            counts_overlap = []
            exp_overlap = []
            maxc =0.0
            maxc_exp = 0.0
            max_gene = ''
            for gene in exps.keys():
                if gene=='' or gene=='n/a': continue
                if gene in XY: continue
                if gene not in classgenes: continue
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
            countsoverlaps.append(counts_overlap)
            expoverlaps.append(exp_overlap)
        genesfile.close()
        countsfile.close()

        #this only does for class 5 for now
        print "The gene with the maximum counts for %s is %s" %(n,max_gene), maxc, maxc_exp 
        rs = []
        ps = []
        for k in range(0,5):
            (r,p) = pearsonr(countsoverlaps[k], expoverlaps[k])
            rs.append(r)
            ps.append(p)
        print 'r = ', rs
        print 'p = ', ps
        plt.figure()
        colors = ['b','g','c','y','r']
        for k in range(4,-1,-1):
        #k = 0
            plt.scatter(countsoverlaps[k], expoverlaps[k], edgecolor='none', alpha=0.3, c=colors[k])
            plt.figtext(0.15,0.85-(.03*k), 'Class %i: r=%.4f, p=%.3g'%(k+1,rs[k],ps[k]),color=colors[k])
            plt.xlabel('Log_2 (H3K4me3 counts)',fontsize='x-large')
            plt.ylabel('Log_2 (Expression) (RPKM)',fontsize='x-large')
            plt.savefig('%s/%s_counts_vs_exp_classes_H3K4me3.pdf'%(countsdir,n))

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
