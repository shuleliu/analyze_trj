#!/usr/bin/python

import math

nrun=6
ndata=160000
tstep=0.5e-6
nrecord=10

e_init=[0.0]*nrun
e_final=[0.0]*nrun
e_drift_rate=[0.0]*nrun

for i in range(0,nrun):
    print "Reading folder %d" %(i+1)
    fname="nve_"+str(i+1)+"/log.out"
    with open(fname,'r') as f:
         for line in f:
             if "Step" in line and "Temp" in line and "TotEng" in line:
                nxtline=f.next()
                args=nxtline.split()
                e_init[i] = float(args[4])
                for j in range(0,ndata):
                    nxtline=f.next()

                args=nxtline.split()
                e_final[i] = float(args[4])

rate_ave=0.0
err_ave=0.0
g=open("energy_drift.txt",'w')
print>>g, "# Run    e_init     e_final    e_drift      e_drift_rate (kcal/mol/ns)"
for i in range(0,nrun):
    e_drift=e_final[i] - e_init[i]
    e_drift_rate[i] = abs(e_drift)/(tstep*float(ndata*nrecord))
    rate_ave+=e_drift_rate[i]
    print>>g, "%d   %12.6f   %12.6f  %12.6f  %12.6f" \
              %(i+1,e_init[i],e_final[i],e_drift,e_drift_rate[i])

rate_ave/=float(nrun)
for i in range(0,nrun):
    err_ave+=(e_drift_rate[i]-rate_ave)**2

err_ave/=float(nrun)
err_ave=math.sqrt(err_ave)

print>>g,"Average Energy Drift Rate is: %8.3f kcal/mol/ns" %rate_ave
print>>g, "The error of energy drift is: %8.3f kcal/mol/ns" %err_ave

g.close()  
    
                   
         
