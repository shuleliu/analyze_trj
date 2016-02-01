#!/usr/bin/python

#calculate the number of states histogram

import sys

inputname=sys.argv[1]

# cutoff for the number of states
ncut=100
state_hist=[0.0]*ncut

nconf=0
with open(inputname,'r') as f:
     for line in f:
         if "state(s)" in line:
            args = line.split()
            nstate=int(args[2])
            if nstate<ncut:  
               state_hist[nstate]+=1.0
            nconf+=1
            if nconf%1000==0:
               print "Reading config %d" %nconf 


#normalize
g = open('state_hist.dat','w')
for i in range(0,ncut):
    state_hist[i]/=float(nconf)
    print>>g, "%d  %f" %(i,state_hist[i])

g.close() 
