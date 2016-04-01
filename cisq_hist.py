#!/bin/python 

# compute the histogram of ci^2 and ci_max^2

import sys
import numpy as np

cisq = []
cimaxsq =[]
cisq_2nd = []

g = open('cisq_hist.dat','w')
h = open('cimaxsq_hist.dat','w')
h2 = open('ci2ndsq_hist.dat','w')

inputfile = sys.argv[1]

with open(inputfile,'r') as f:
     for line in f:
         if "EIGEN_VECTOR" in line:
            nxtline = f.next()
            args = nxtline.split()
            for item in args:
                cisq.append(float(item)*float(item))

            cimaxsq.append(float(args[0])*float(args[0]))
            if len(args)>1:
               cisq_2nd.append(float(args[1])*float(args[1]))

cisq = np.array(cisq)
cimaxsq = np.array(cimaxsq)
cisq_2nd = np.array(cisq_2nd)

binsq = np.linspace(0.0,1.05,num=100)

cisq_hist, bin_cisq = np.histogram(cisq,bins = binsq)
cimaxsq_hist, bin_cimaxsq = np.histogram(cimaxsq,bins = binsq)
cisq2nd_hist, bin_cisq2nd = np.histogram(cisq_2nd,bins = binsq)

for i in range(0,bin_cisq.size):
      print>>g, "%.6f  %.6f" %(bin_cisq[i],cisq_hist[i]) 
      print>>h, "%.6f  %.6f" %(bin_cimaxsq[i],cimaxsq_hist[i])
      print>>h2, "%.6f  %.6f" %(bin_cisq2nd[i],cisq2nd_hist[i])

g.close()
h.close()
h2.close() 
