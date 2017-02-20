#!/usr/bin/env python

#Chris Ng

from __future__ import division
import sys, math, cPickle,os
import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['pdf.fonttype']=42
from optparse import OptionParser
from collections import defaultdict
from matplotlib.pyplot import *
from numpy import array, zeros, loadtxt, save,savetxt,arange
from Pycluster import kcluster, distancematrix
import numpy as np

def kmeans(names,tsstags,tags,outdir,norm=True,mi=-2000,ma=6000,bin=100,k=5):
    Xall=[]
    N=names.split(',')
    if len(tsstags.split(','))==1: TSS=[int(tsstags)]*len(N)
    else: TSS=[int(x) for x in tsstags.split(',')]
    if len(tags.split(','))==1: TAGS=[int(tags)]*len(N)
    else: TAGS=[int(x) for x in tags.split(',')]
    INDEX=[]
    for n,tss,t in zip(N,TSS,TAGS):
        fi="../results/tss_profiles/%s/TSStags%i_tags%i_bin%i_avg200_pdist%i_profiles.txt"%(n,tss,t,bin,ma)
        X = loadtxt(fi)
        Xall.append(X)
        gfi="../results/tss_profiles/%s/TSStags%i_tags%i_bin%i_avg200_pdist%i_genes.txt"%(n,tss,t,bin,ma)

        for line in open(gfi):
            sym=line.strip('\n')
            INDEX.append((n,sym))
    mi=mi+1.5*bin
    ma=ma-bin
    INDEX=array(INDEX)
    Xcat=np.vstack(Xall)
    if norm: X_norm = (Xcat.T/Xcat.sum(axis=1).T).T
    print X_norm
    #cond_name = fn[0:fn.find('_')]
    #figure()
    figure()
    #Specify the corresponding gene list for the profile array
    #if genes:
    #    #G=[]
    #    Gfi=genes.split(',')
    #    for name,fi in zip(N,Gfi):
    #            #G.append(line.strip('\n'))
    #        #G=array(G)
    #    #INDEX.append(G)
    #    INDEX=array(INDEX)
    k_clusts = [int(x) for x in k.split(',')]
    for i, nclust in enumerate(k_clusts):
        clusterid, error, nfound = kcluster(X_norm,nclusters=nclust,npass=10)
        print 'nclusters=%d'%nclust
        print clusterid
        mean_clusts = []
        num_in_clusts = []
        
        #CIDS=set(clusterid)
        CIDS=xrange(nclust)

        #This was a quick/dirty solution to rank the clusters by 5' density
        RANKS,SCORES={},[]
        for cid in CIDS:
            print cid,
            mean_clust = X_norm[clusterid==cid].mean(axis=0)
            score=mean_clust[:30].sum()
            SCORES.append((score,cid))
        SCORES.sort()
        RANKS=[x[1] for x in SCORES]
        R={}
        for ind,cid in enumerate(RANKS): R[cid]=ind+1
        #Reorganize clusterid array by cluster ranks
        clusranks=array([R[cid] for cid in clusterid])
        
        new_cids=xrange(1,nclust+1)
        for name in N:
            name_dir='%s/%s'%(outdir,name)
            try: os.mkdir(name_dir)
            except:
                exception=1
        for cid in new_cids:
            mean_clust = X_norm[clusranks==cid].mean(axis=0)
            mean_clusts.append(mean_clust)
            num_in_clusts.append((clusranks==cid).sum())
            for name in N:
                IN_CLUS = INDEX[clusranks==cid]
                G_out=[]
                for n,sym in IN_CLUS:
                    if name==n: G_out.append(sym)
                name_dir='%s/%s'%(outdir,name)
                savetxt('%s/k%d_%d.txt'%(name_dir,nclust,cid),G_out,fmt='%s')
            
        save('%s/k%d.npy'%(outdir,nclust),array(mean_clusts))
        save('%s/k%d_ids.npy'%(outdir,nclust),array(clusranks))
        ax=subplot(len(k_clusts),1,i+1)
        xax=arange(mi,ma,bin)
        colors=['blue','magenta','cyan','yellow','red']
        for mean_clust, num_in_clust,cid,color in zip(mean_clusts,num_in_clusts,new_cids,colors):
            #print xax,mean_clust
            print xax.shape,mean_clust.shape
            ax.plot(xax,mean_clust,label='%d'%(cid),linewidth=2,color=color)
        ax.set_xlabel('Distance to TSS',fontsize='xx-large')
        ax.set_ylabel('Read density',fontsize='xx-large')
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize('large')
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize('large')
        ax.grid(True)
        legend(loc=1, prop={'size':18},markerscale=0.5)

    savefig('%s/kmean_clusts.pdf'%(outdir))

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--ma", type=int, default=6000, dest="ma")
    parser.add_option("--mi", type=int, default=-2000, dest="mi")
    parser.add_option("--TSS",  dest="TSS")
    #parser.add_option("--mark", dest="mark")
    parser.add_option("--tags", dest="tags")
    parser.add_option("--bin", type=int, default=100, dest="bin")
    parser.add_option("--k", type=str, default="5", dest="k")
    parser.add_option("--outfile", dest="outfile", default=None,help='Output prefix')
    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    names=args[0]

    kmeans(names,opts.TSS,opts.tags,opts.outfile,mi=opts.mi,ma=opts.ma,bin=opts.bin,k=opts.k)

if __name__=='__main__': main()
