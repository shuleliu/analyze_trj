#!/usr/bin/python

# extract off-diagonal coupling strength to make histogram

import sys

inputname=sys.argv[1]

off_upper=2.0
off_lower=-30.0
deltabin=0.1
nbin=int((off_upper-off_lower)/deltabin)

off_hist=[0.0]*nbin

nconf=0
with open(inputname,'r') as f:
     for line in f:
         if "OFF-DIAGONAL" in line:
            nxtline = f.next()
            args = nxtline.split()
            fst_cp = float(args[1])
            ibin = int((off_upper-fst_cp)/deltabin)
            if ibin<nbin and ibin>=0:
               off_hist[ibin]+=1.0

            nconf+=1
            if nconf%1000==0:
               print "Reading config %d" %nconf


# normalize
g = open('off_hist.dat','w')
for i in range(0,nbin):
    off_hist[i]/=float(nconf)
    print>>g, "%f  %f" %(-deltabin*float(i)+off_upper,off_hist[i])

g.close()
             
