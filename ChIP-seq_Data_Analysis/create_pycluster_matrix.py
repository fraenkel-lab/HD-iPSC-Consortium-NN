import sys,math
from optparse import OptionParser
import numpy as np

def multcond(Mfi,a,b,outfile,log2,cat=None):
    header=a.split(',')
    out=open(outfile,'w')
    out.write('\t'.join(['genes']+header)+'\n')
    B=np.loadtxt(b,dtype=str)
    M=np.loadtxt(Mfi,dtype=float)
    print len(B),len(M)
    for i,row in enumerate(M):
        if log2: row=np.log2(row+1)
        row_name=B[i]
        out.write('\t'.join([row_name]+["%.3f"%x for x in row])+'\n')
    out.close()

def main():
    usage = "usage: %prog [opts] arg"
    parser = OptionParser(usage)
    parser.add_option("--log2", action="store_true", dest="log2")
    parser.add_option("--cat", type=int, dest="cat")
    parser.add_option("--outfile", dest="outfile")
    (opts, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    
    M=args[0]
    a=args[1]
    b=args[2]
    
    multcond(M,a,b,opts.outfile,opts.log2,cat=opts.cat)
    
if __name__=='__main__':
  main()
