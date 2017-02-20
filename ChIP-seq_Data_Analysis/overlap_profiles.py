from __future__ import division
import sys, math, cPickle,os
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype']=42
from optparse import OptionParser
from collections import defaultdict
from TAMO.util.Arith import hypgeomsummore
sys.path+=['/nfs/vendata/cwng/apps/fishers_exact_test','/nfs/vendata/cwng/apps/fishers_exact_test/src','/nfs/vendata/cwng/apps/fisher-0.1.4/fisher-0.1.4-py2.6-linux-x86_64.egg/fisher/']
from fisher import pvalue

def bobby_fisher(a,b,c,d):
    P = pvalue(a,b,c,d)
    #return P.right_tail
    return P.left_tail,P.right_tail

def contingency_table(C,G,ALL_CLUS,G_bg):
    a=len(C.intersection(G))
    b=len(C.intersection(G_bg))-a
    c=len(ALL_CLUS.intersection(G))-a
    d=len(ALL_CLUS.intersection(G_bg))-a-b-c
    return a,b,c,d

def f(a,b,c,k,outfile,title=None,ylabel=None):
    genes = np.array(open(a).readlines(),dtype='string')
    G=set(genes)
    genes_bg = np.array(open(c).readlines(),dtype='string')
    G_bg=set(genes_bg)
    over=[]
    n_clus=[]
    ALL_CLUS=set()
    CLUS=[]
    for i in xrange(k):
        clus_file=b+str(i+1)+'.txt'
        clus_genes=np.array(open(clus_file).readlines(),dtype='string')
        C=set(clus_genes)
        ALL_CLUS.update(C)
        CLUS.append(C)
        intersect=C.intersection(G)
        over.append(len(intersect))
        n_clus.append(len(C))
    tot_clus=np.sum(n_clus)
    over=np.array(over)
    #n_g=len(G)
    n_g=np.sum(over)
    E,O=[],[]
    P=[]
    U=[]
    for i,o,n_c,C in zip(xrange(k),over,n_clus,CLUS):
        e=n_g*(n_c/tot_clus)
        print 'In Class %d, '%(i+1)
        print 'Expected: %d, '%e,
        print 'Observed: %d,'%o
        #M=tot_clus
        #N=n_g
        #n=n_c
        #x=o
        #pval=hypgeomsummore(n,M,N,x)
        a,b,c,d=contingency_table(C,G,ALL_CLUS,G_bg)
        left,right=bobby_fisher(a,b,c,d)
        pval=right
        print a,b,c,d
        print left,right
        print 'Hypergeometric: %.3e'%pval
        E.append(e)
        O.append(o)
        P.append(right)
        U.append(left)
    
    ind = np.arange(k)  # the x locations for the groups
    width = 0.35       # the width of the bars
    
    colors=['blue','magenta','cyan','yellow','red','0.75','yellow','orange','magenta','pink']
    fig = plt.figure()
    ax = fig.add_subplot(111)
#    hide1= ax.bar(ind[0],E[0],width,color='white',hatch='/')
#    hide2= ax.bar(ind[0]+width,O[0],width,color='black')
    hide1= ax.bar(ind,E,width,color='w',hatch='/',alpha=0.7)
    hide2= ax.bar(ind[0]+width,O[0],width,color='w')
    #rects1 = ax.bar(ind, E, width, color='white',hatch='/')
    rects1 = ax.bar(ind, E, width, color=colors,hatch='/',alpha=0.3)
    rects2 = ax.bar(ind+width, O, width, color=colors)

    # add some
    ax.set_xlabel('Class',fontsize='xx-large')
    if ylabel: 
        yl=ylabel.replace('_',' ')
        ax.set_ylabel(yl,fontsize='xx-large')
    else: ax.set_ylabel('# genes w/ lower expression in HD',fontsize='xx-large')
    ax.set_xticks(ind+width)
    #ax.set_xticklabels( ('0', '1', '2', '3', '4') )
    ax.set_xticklabels( [str(x+1) for x in xrange(k)] )

    #extend ylim by 2 ticks
    Ylim=ax.get_ylim()
    Ytick=ax.get_yticklabels()
    n=len(Ytick)-1
    Ydist=Ylim[1]-Ylim[0]
    print Ylim,Ydist,n
    ax.set_ylim((Ylim[0],Ylim[1]+Ydist/n))

    def autolabel(rects,P):
        # attach some text labels
        i=0
        for rect,p,u in zip(rects,P,U):
            
            height = rect.get_height()
            if p<0.001: 
                ax.text(rect.get_x()+rect.get_width()/2., 1.15*height, 'p=%.1e'%p,
                            ha='center', va='bottom',fontsize='xx-large',color=colors[i])
                ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '**',
                            ha='center', va='bottom',fontsize='xx-large')
            elif u<0.001:
                ax.text(rect.get_x()+rect.get_width()/2., 1.15*height, 'p=%.1e'%u,
                            ha='center', va='bottom',fontsize='xx-large',color=colors[i])
                ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '*',
                            ha='center', va='bottom',fontsize='xx-large')
            i+=1

    autolabel(rects2,P)

    ax.legend( (hide1[0], hide2), ('Expected', 'Observed') , prop={'size':20})
#    ax.legend( (rects1[0], rects2[0]), ('Expected', 'Observed') )
    if title: ax.set_title(title.replace('_',' '),fontsize='large')
    fig.savefig(outfile)
    
def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--k", dest="k", type=int,default=5)
    parser.add_option("--ylabel", dest="ylabel", default=None)
    parser.add_option("--title", dest="title", default=None)
    parser.add_option("--outfile", dest="outfile", default=None)
    (opts, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")

    a=args[0]
    b=args[1]
    c=args[2]

    f(a,b,c,opts.k,opts.outfile,title=opts.title,ylabel=opts.ylabel)

if __name__=='__main__': main()
