#Chris Ng

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

def f(a,outfile,countsfile,sample,log=False,bins=40,xlabel=None,mask=None,mi=None,ma=None):
    A=np.loadtxt(a,dtype=float)
    counts = open(countsfile, 'rb')
    
    countsdict = {}
    for line in counts:
        words = line.strip().split()
        countsdict[words[0]] = words[1]
    count = float(countsdict[sample])
    
    if log:
        A=np.divide(A,(count/1000000.0))
        A=np.log2(np.add(A,1))
        np.savetxt(outfile.replace('.pdf','.txt'),A,fmt='%.4f')

    if mask:
        MASK=np.loadtxt(mask,dtype=int)
        #mdat = np.ma.masked_array(dat,np.isnan(dat))
        #A=np.ma.array(A,dtype=float,mask=MASK)
        A=np.ma.masked_array(A,dtype=float,mask=np.isnan(A))
        #A=np.array(A.compressed())
    plt.figure()
    if ma: plt.hist(A,bins=bins,range=(mi,ma))
    else: plt.hist(A,bins=bins)
    if log: plt.xlabel('log2 read counts')
    elif xlabel: plt.xlabel(xlabel.replace('_',' '))
    else: plt.xlabel('# of read counts')
    plt.ylabel('#')
    plt.savefig(outfile)

    print A.shape
    print A
    mean,std=A.mean(),A.std()
    print mean,std
    mm = np.mean(A)
    #print mm.filled(np.nan)
    ms = np.std(A)
    #print ms.filled(np.nan)
    print mm,ms
    mean,std=mm,ms

    out=open(outfile.replace('.pdf','summary.txt'),'w')
    out.write('Mean\tStdDev\n')
    out.write('%.2f\t%.2f\n'%(mean,std))

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--mask", dest="mask")
    parser.add_option("--log", action="store_true", dest="log")
    parser.add_option("--counts", dest="countsfile")
    parser.add_option("--sample", dest="sample")
    #parser.add_option("--thres", dest="thres", type=float, default=0.7)
    parser.add_option("--mi", dest="mi", type=float, default=None)
    parser.add_option("--ma", dest="ma", type=float, default=None)
    parser.add_option("--bins", dest="bins", type=int, default=40)
    parser.add_option("--xlabel", dest="xlabel", default=None)
    parser.add_option("--outfile", dest="outfile", default=None)
    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    a=args[0]

    f(a,opts.outfile,opts.countsfile,opts.sample,log=opts.log,bins=opts.bins,xlabel=opts.xlabel,mask=opts.mask,mi=opts.mi,ma=opts.ma)

if __name__=='__main__': main()
