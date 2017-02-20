#Amanda Daigle
#4/25/14

from optparse import OptionParser
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype']=42

def count(clusterdir,outdir,names,k):
    k = int(k)
    #list of lists for each cell line
    countsdata = []
    for n in names.split(','):
        #List for each cellline with number of genes in order of clusterid
        counts = []
        for i in range(1,k+1):
            with open(clusterdir+'/%s/k%s_%s.txt'%(n,k,i), 'rb') as cf:
                for g,l in enumerate(cf):
                    pass
            count = g+1
            counts.append(count)
        countsdata.append(counts)
    
    #make table
    fig=plt.figure()
    ax = fig.add_subplot(111)
    #ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)
    ax.axis('off')
    colLabels=['Class %i' %x for x in range(1,k+1)]
    rows=names.split(',')
    the_table = ax.table(cellText=countsdata, 
                         colLabels=colLabels, 
                         rowLabels=rows,
                         loc='center')
    plt.title('Number of genes assigned to each cluster')
    plt.savefig("%s/counts_table.pdf"%outdir)


def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--clusterdir", dest="clusterdir", default='../results/kmeans_cluster')
    parser.add_option("--outdir", dest="outdir", default="../results/kmeans_cluster")
    parser.add_option("--names", dest="names", default='HD109n1,HD109n4,HD109n5,HD180n1,HD180n3,HD21n1,HD28n6,HD33n1,HD60n5,HD60n8')
    parser.add_option("-k", dest="k", default=5)
    
    (opts, args) = parser.parse_args()
    
    count(opts.clusterdir,opts.outdir,opts.names,opts.k)
    
if __name__=='__main__': main()
