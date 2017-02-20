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

def TSS_dist(sam_files,TSSfile,genome,outdir,pos,genes=False,shift=90,format='SAM'):
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
    #R=range(mdist,pdist+1)
    get_dist(SAM,TSS,outdir,n_reads,genes=genes)
        
def get_dist(SAM,TSS,outdir,n_reads,genes=False):
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
    T=[]
    GENES=open(outdir+'genes_norm.txt','w')
    for chr in TSS.keys():
        if chr not in SAM.keys(): continue
        for tuple in TSS[chr]:
            tss,tse,gene,s=tuple
            if genes:
                if gene not in G: continue
            if s=='+': genelength=tse-tss
            elif s=='-': genelength=tss-tse
            R=np.arange(0,genelength+1)
            #if genelength<bins: continue
            r=calc_read_dist(SAM[chr],R,tss,s)
            tags=sum(r.values())
            T.append(tags/float(genelength))
            n_tss+=1
            GENES.write(gene+'\n')
    print n_tss
    outfile=outdir+'counts_norm.txt'
    np.savetxt(outfile,T,fmt='%f')

def calc_read_dist(BED,R,event,s):
    r={}
    if s=='+': C=[event+x for x in R]
    elif s=='-': C=[event-x for x in R]
    for x,c in zip(R,C):
        if c in BED: r[x]=BED[c]
        else: r[x]=0
    return r

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--genome", dest="genome", default='mm9')
    parser.add_option("--shift", dest="shift",type=int, default=100)
    parser.add_option("--genes", dest="genes", default=0)
    parser.add_option("--pos", dest="pos", action="store_true")
    parser.add_option("--format", dest="format",default='SAM')
    parser.add_option("--TSS", dest="TSSfile")
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    sam_file=args[0]
    outdir=args[1]
    
    genome=options.genome
    TSS_dist(sam_file,options.TSSfile,genome,outdir,options.pos,genes=options.genes,shift=options.shift,format=options.format)

if __name__=='__main__':
  main()
