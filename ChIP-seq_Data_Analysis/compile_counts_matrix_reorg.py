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

def f(a,outfile,genes_in,genes_out):
    #load each counts array and append to matrix
    Alist=a.split(',')
    for i,file in enumerate(Alist):
        A=np.loadtxt(file,dtype=float)
        if i==0: M=np.transpose([A])
        else:
            print M.shape,A.shape
            M=np.hstack((M, np.transpose([A])))
    
    Fin=np.loadtxt(genes_in,dtype=str)
    Fout=np.loadtxt(genes_out,dtype=str)
    D={}
    for i,g in enumerate(Fin):
        D[g]=M[i]
    Mnew=[]
    skipped_genes=0
    for g in Fout:
	try:
            Mnew.append(D[g])
        except KeyError:
            skipped_genes += 1
    print "Skipped %i genes" %skipped_genes
    np.savetxt(outfile,Mnew,fmt='%.3f')
    
def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    #parser.add_option("-n", action="store_true", dest="n")
    #parser.add_option("--thres", dest="thres", type=float, default=0.7)
    #parser.add_option("--dist", dest="dist", type=int, default=1)
    parser.add_option("--genes_in", dest="genes_in", default=None)
    parser.add_option("--genes_out", dest="genes_out", default=None)
    parser.add_option("--outfile", dest="outfile", default=None)
    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    a=args[0]

    f(a,opts.outfile,opts.genes_in,opts.genes_out)

if __name__=='__main__': main()
