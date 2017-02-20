from __future__ import division
import sys, math, cPickle,os
import numpy as np
from optparse import OptionParser
from collections import defaultdict
#from Pycluster import *
import Pycluster

def f(a,outfile,log2,log10,clusrow,cluscol,scale,m='m',rdist='c',cdist='c',mask=None):
    #make record object
    handle = open(a)
    record = Pycluster.read(handle)
    #log transform
    if log2:
        record.data=np.log2(record.data+np.ones(record.data.shape))
    if log10:
        record.data=np.log10(record.data+np.ones(record.data.shape))
    #hierarchical cluster record
    if clusrow:
        print "clusterNg rows..."
        rowtree = record.treecluster(transpose=0, method=m, dist=rdist)
    if cluscol:
        print "clusterNg columns..."
        coltree = record.treecluster(transpose=1, method=m, dist=cdist)
    #save record to jtreeview file
    print "scaleNg and saveNg..."
    if clusrow and cluscol:
        if scale:
            rowtree.scale()
            coltree.scale()
        record.save(outfile, rowtree, coltree)
    elif clusrow:
        if scale: rowtree.scale()
        record.save(outfile, geneclusters=rowtree)
    elif cluscol: 
        if scale: coltree.scale()
        record.save(outfile, expclusters=coltree)

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--scale", action="store_true", dest="scale")
    parser.add_option("--clusrow", action="store_true", dest="clusrow")
    parser.add_option("--cluscol", action="store_true", dest="cluscol")
    parser.add_option("--log2", action="store_true", dest="log2")
    parser.add_option("--log10", action="store_true", dest="log10")
    #parser.add_option("--log", type=int, dest="log")
    parser.add_option("--mask", dest="mask", default=None)
    parser.add_option("--rdist", dest="rdist", default='c')
    parser.add_option("--cdist", dest="cdist", default='c')
    parser.add_option("--m", dest="m", default='m')
    parser.add_option("--outfile", dest="outfile", default=None)
    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    a=args[0]

    f(a,opts.outfile,opts.log2,opts.log10,opts.clusrow,opts.cluscol,opts.scale,m=opts.m,rdist=opts.rdist,cdist=opts.cdist,mask=opts.mask)

if __name__=='__main__': main()
