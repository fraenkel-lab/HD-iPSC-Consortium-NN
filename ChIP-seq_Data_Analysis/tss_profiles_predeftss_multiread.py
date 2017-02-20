#Chris Ng

from __future__ import division
import sys, math, cPickle,os
import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from optparse import OptionParser
from collections import defaultdict

def bed5p2dict(bed,shift=85,pos=False,BED={},mrc=3):
    bedo=open(bed)
    #BED={}
    print 'Assembling 5prime reads from bedfile'
    for line in bedo:
        L=line.strip('\n').split('\t')
        chr=L[0]
        #try: start,stop,seq,strand=int(L[1]),int(L[2]),L[3],L[5]
        try: start,stop,seq,strand=int(L[1]),int(L[2]),L[3],L[5]
        except:
            try: start,stop,seq,strand=int(L[1]),int(L[2]),L[3],L[4]
            except:
                try:
                    L=line.split(' ')
                    start,stop,seq,strand=int(L[1]),int(L[2]),L[3],L[4]
                    #start,stop,seq,strand=int(L[1]),int(L[2]),L[3],L[5]
                except:
                    print 'line %s skipped'%(line)
                    continue
        if strand=='+':
            start=start+shift
        elif strand=='-':
            if pos: continue
            start=stop-1-shift
        else:
            print 'strand not found'
            continue
        if chr not in BED.keys():   BED[chr]={}
        if start in BED[chr]: 
            if BED[chr][start]==mrc: continue
            else: BED[chr][start]+=1
        else: BED[chr][start]=1
        #BED[chr][start]=1
    uniq_reads=0
    for chr in BED.keys(): uniq_reads+=len(BED[chr].keys())
    bedo.close()
    return BED, uniq_reads

def sam5p2dict(sam,shift=90,pos=False,swap=False,SAM=defaultdict(dict),mrc=3):
    samo=open(sam)
    #SAM=defaultdict(dict)
    print 'Assembling 5prime reads from samfile'
    for line in samo:
        if line[0]=='@': continue
        L=line.strip('\n').split('\t')
        flag=L[1]
        if flag=='0': strand='+'
        elif flag=='16': strand='-'
        elif flag=='4': continue
        else: 
            print 'strand %s not found'%flag
            continue
        chr,start=L[2],int(L[3])
        if strand=='+':
            start=start+shift
        elif strand=='-':
            if pos: continue
            start=start+36-shift
        else: continue
        if start in SAM[chr]:
            if SAM[chr][start]==mrc: continue
            else: SAM[chr][start]+=1
        else: SAM[chr][start]=1
        #SAM[chr][start]=1
    uniq_reads=0
    for chr in SAM.keys(): uniq_reads+=len(SAM[chr].keys())
    samo.close()
    return SAM, uniq_reads

def TSS_dist(sam_files,TSSfile,genome,outdir,mdist=-2000,pdist=5000,pos=False,genes=False,shift=90,bin=20,min_tags=40,mv_avg=False,norm=False,format='SAM'):

    if format=='SAM': 
        SAM=defaultdict(dict)
        sf=sam_files.split(',')
        for sam_file in sf:
            SAM,n_reads=sam5p2dict(sam_file,pos=pos,shift=shift,SAM=SAM)
    else:
        SAM={}
        sf=sam_files.split(',')
        for sam_file in sf:
            SAM,n_reads=bed5p2dict(sam_file,pos=pos,shift=shift,BED=SAM)

    if TSSfile: TSS=cPickle.load(open(TSSfile))
    else:
        if genome=='mm9':   TSS=cPickle.load(open('/nfs/vendata/cwng/pickles/mm9_tss_strand.pickle'))
        else:   TSS=cPickle.load(open('/nfs/vendata/cwng/pickles/hg18_tss_strand.pickle'))
    R=range(mdist,pdist+1)
    #bins=np.arange(mdist,pdist,bin)
    #if not os.path.isdir(outdir): os.mkdir(outdir)
    get_dist(SAM,TSS,R,outdir,n_reads,genes=genes,bins=bin,min_tags=min_tags,mv_avg=mv_avg,norm=norm)
        
def get_dist(SAM,TSS,R,outdir,n_reads,genes=False,bins=None,min_tags=40,mv_avg=False,norm='density'):
    print 'Assembling read dist adj to matches'
    if genes:
        G=[]
        FH=open(genes)
        for line in FH:
            gene=line.strip('\n').split('\t')[0]
            G.append(gene)
    n_tss=0
    P={}
    TAGS={}
    GENES=open(outdir+'genes.txt','w')
    for chr in TSS.keys():
        if chr not in SAM.keys(): continue
        for tuple in TSS[chr]:
            #tse is the end site
	    tss,tse,gene,s=tuple
            if genes:
                if gene not in G: continue
            r=calc_read_dist(SAM[chr],R,tss,s)
#            tags_adj=sum([r[x] for x in xrange(-1000,1000)])
#            if tags_adj<min_tags: continue
            tags=sum(r.values())
#            if gene in TAGS:
#                if TAGS[gene]>tags_adj: continue
#                else: TAGS[gene]=tags_adj
#            else: TAGS[gene]=tags_adj
            if tags<min_tags: continue
            n_tss+=1
            if norm=='density': r_vect=np.array([r[x]/tags for x in R])
            else: r_vect=np.array([r[x] for x in R],dtype=float)
            #elif norm=='max': r_vect/=np.sum(r_vect)
            #if mv_avg: r=moving_average_dict(r,mv_avg)
            if mv_avg: r_vect=moving_average(r_vect,mv_avg)
            #bincounts=BIN(r,bins,mv_avg=mv_avg)
            bincounts=BIN(r_vect,bins)
            #PROF.append(bincounts)
            P[gene]=bincounts
            #GENES.write(gene+'\n')
    print n_tss
    PROF=np.array(P.values())
    outfile=outdir+'profiles.txt'
    np.savetxt(outfile,PROF,delimiter='\t')
    for g in P.keys(): GENES.write(g+'\n')
#    plt.figure()
#    plt.hist(TAGS.values(),bins=40)
#    plt.xlabel('Tags per TSS window')
#    plt.ylabel('#')
#    plt.savefig(outdir+'tags_per_tss_hist.pdf')

def BIN(X,bins,mv_avg=False):
#    binned=np.zeros(len(bins))
    #bin_index=0
    k=len(X)//bins
#    for i in xrange(bins):
    binned=[np.sum(X[(i*bins):((i+1)*bins)]) for i in xrange(k)]
    #keys=X.keys()
    #keys.sort()
    #for k in keys:
        #v=X[k]
        #if bin_index==len(bins)-1: binned[bin_index]+=v
        #elif k<bins[bin_index+1]: binned[bin_index]+=v
        #elif k==bins[bin_index+1]:
        #    binned[bin_index]+=v
        #    bin_index+=1
    if mv_avg:
        binned_avg=moving_average(binned,mv_avg)
        return binned_avg
    else: return binned

def moving_average(X,k):
    mv_avg=np.zeros(len(X)-k+1)
    s= sum(X[0:k])
    mv_avg[0]=float(s/k)
    #yield s / k
    for i,x in enumerate(range(k, len(X))):
        s += X[x] - X[x-k]
        #yield s / k
        mv_avg[i+1]=float(s/k)
    #for i in np.arange(len(mv_avg)): mv_avg[i]=np.mean(X[i:i+k])
    return mv_avg

def moving_average_dict(X,k):
    mv_avg={}
    keys=X.keys()
    keys.sort()
    for i in keys[0:-k+1]:
        mv_avg[i]=np.mean([X[x] for x in np.arange(i,i+k)])
    return mv_avg

def calc_read_dist(BED,R,event,s,rank=False):
    r={}
    if s=='+': C=[event+x for x in R]
    elif s=='-': C=[event-x for x in R]
    for x,c in zip(R,C):
        if c in BED: r[x]=BED[c]
        else: r[x]=0
    #if rank:    
    #    rt=rank_trans(r)
    #    return rt
    #else:
    return r

def rank_trans(READS):
    d=defaultdict(list)
    d_rank={}
    for key,value in READS.items(): d[value].append(key)
    tags=d.keys()
    tags.sort()
    for rank,key in enumerate(tags):
        for coord in d[key]:
            d_rank[coord]=rank
    return d_rank

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--genome", dest="genome", default='mm9')
    parser.add_option("--mdist", dest="mdist", type=int, default=-2000)
    parser.add_option("--pdist", dest="pdist", type=int, default=3000)
    parser.add_option("--pos", dest="pos", action="store_true")
    parser.add_option("--norm", dest="norm")
    parser.add_option("--shift", dest="shift",type=int, default=100)
    parser.add_option("--genes", dest="genes")
    parser.add_option("--bin", dest="bin",type=int,default=20)
    parser.add_option("--mv_avg", dest="mv_avg",type=int)
    parser.add_option("--min_tags", dest="min_tags",type=int,default=40)
    parser.add_option("--format", dest="format",default='SAM')
    parser.add_option("--TSS", dest="TSSfile")
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    sam_file=args[0]
    outdir=args[1]
    
    genome=options.genome
    TSS_dist(sam_file,options.TSSfile,genome,outdir,mdist=options.mdist,pdist=options.pdist,pos=options.pos,genes=options.genes,shift=options.shift,bin=options.bin,min_tags=options.min_tags,mv_avg=options.mv_avg,norm=options.norm,format=options.format)

if __name__=='__main__':
  main()
