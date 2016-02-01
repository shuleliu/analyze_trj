#!/usr/bin/python

#calculate the rdf between CEC and water

import math

def PBC_dis(x,y,z,r_cec,boxlength):
    x12=r_cec[0]-x
    y12=r_cec[1]-y
    z12=r_cec[2]-z
    x12=x12-boxlength*float(round(x12/boxlength))
    y12=y12-boxlength*float(round(y12/boxlength))
    z12=z12-boxlength*float(round(z12/boxlength))
    r12=math.sqrt(x12*x12+y12*y12+z12*z12)

    return r12

nwindow=49
natom=1466
nwater=486
boxlength=24.3794
deltar=0.05
nbin=240

rhoOW=float(nwater)/(boxlength**3)
rhoHW=rhoOW*2.0
f=open('cec-ow-rdf.dat','w')
g=open('cec-hw-rdf.dat','w')
print>>f, "# rwindow   rrr    g(r)"
print>>g, "# rwindow   rrr    g(r)"

cec_step=[]
cec_coord=[]

match=False
for i in range(0,nwindow):
    print "Processing window %d" %i
    fname="wndow_"+str(i)+"/"
    trjname=fname+"bin_"+str(i)+".0.lammpstrj"
    cecname=fname+"cec_coor_"+str(i)+".dat"
    rwindow=0.25*float(i)
    rdfname=fname+"cec-ow-hw-rdf.dat"
    h=open(rdfname,'w')
    print>>h, "# rrr   OW-rdf   HW-rdf"
    rdf_OW=[0.0]*nbin
    rdf_HW=[0.0]*nbin
    with open(cecname,'r') as cn:
         for line in cn:
             if "#" not in line:
                args=line.split()
                tstep=int(args[0])
                cec_step.append(tstep)
                coord=[float(args[1]),float(args[2]),float(args[3])]
                cec_coord.append(coord)
                
    print "Finish reading the CEC coordinate file"
    
    ind_step=0
    r_cec=[0.0,0.0,0.0]
    nconf=0
    with open(trjname,'r') as trj:
         for line in trj:
             if "TIMESTEP" in line:
                nxtline=trj.next()
                istep=int(nxtline)
                if istep%10000==0:
                   print "Reading step %d" %istep
                if istep in cec_step:
                   ind_step=cec_step.index(istep)
                   match=True
                   nconf+=1
                elif istep-1 in cec_step:
                   ind_step=cec_step.index(istep-1)
                   match=True
                   nconf+=1

                r_cec=cec_coord[ind_step]
 
             if "id" in line and match:
                for j in range(0,natom):
                    nxtline=trj.next()
                    args=nxtline.split()
                    atm_type=int(args[2])
                    xxx=float(args[4]) 
                    yyy=float(args[5]) 
                    zzz=float(args[6])
                    if atm_type==1:
                       rrr=PBC_dis(xxx,yyy,zzz,r_cec,boxlength)
                       ibin=int(rrr/deltar)
                       if ibin<nbin:
                          rdf_OW[ibin]+=1.0 
                    elif atm_type==2:
                       rrr=PBC_dis(xxx,yyy,zzz,r_cec,boxlength)
                       ibin=int(rrr/deltar)
                       if ibin<nbin:
                          rdf_HW[ibin]+=1.0 
                
                match=False
#                ind_step=0
    
    for ibin in range(0,nbin):
        rrr=deltar*(float(ibin)+0.5) 
        rdf_OW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoOW*float(nconf)
        print>>f, "%.6f  %.6f  %.6f" %(rwindow,rrr,rdf_OW[ibin])
        rdf_HW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoHW*float(nconf)
        print>>g, "%.6f  %.6f  %.6f" %(rwindow,rrr,rdf_HW[ibin])
        print>>h, "%.6f  %.6f  %.6f" %(rrr,rdf_OW[ibin],rdf_HW[ibin])

    h.close() 
    del cec_step[:]
    del cec_coord[:]
    del rdf_OW[:]
    del rdf_HW[:]

f.close()
g.close()

