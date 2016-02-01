#!/usr/bin/python

#compute RDFs for carbonate 
#the code is supposed to read the lammps trajectory file
#for the window when carbonate is in its deprotonated state

import math
import sys

def PBC_dis(r0,rrr,boxlength):
    x12=rrr[0]-r0[0]
    y12=rrr[1]-r0[1]
    z12=rrr[2]-r0[2]
    x12=x12-boxlength*float(round(x12/boxlength))
    y12=y12-boxlength*float(round(y12/boxlength))
    z12=z12-boxlength*float(round(z12/boxlength))
    r12=math.sqrt(x12*x12+y12*y12+z12*z12)

    return r12

trjname=sys.argv[1]
outputname=sys.argv[2]
natom=1466
boxlength=24.3794
nwater=486
deltar=0.05
nbin=240

rhoOW=float(nwater)/(boxlength**3)
rhoHW=rhoOW*2.0

nconf=0
nbi=0 #check if the number of times HCO3 shows up is consistent with nconf

rdf_CA_OW=[0.0]*nbin
rdf_CA_HW=[0.0]*nbin
rdf_HH_OW=[0.0]*nbin
rdf_HH_HW=[0.0]*nbin
rdf_OH_OW=[0.0]*nbin
rdf_OH_HW=[0.0]*nbin
rdf_OA_OW=[0.0]*nbin
rdf_OA_HW=[0.0]*nbin

OW_coord=[]
HW_coord=[]
OA_coord=[]
CA_coord=[0.0,0.0,0.0]
OH_coord=[0.0,0.0,0.0]
HH_coord=[]

with open(trjname,'r') as trj:
     for line in trj:
         if "id" in line:
            nconf+=1
            print "Reading config %d" %nconf
#           if nconf<1000:
            for j in range(0,natom):
                nxtline=trj.next()
                args=nxtline.split()
                atm_type=int(args[2])
                xxx=float(args[4])
                yyy=float(args[5])
                zzz=float(args[6])
                if atm_type==1:
                   OW_coord.append([xxx,yyy,zzz])
                elif atm_type==2:
                   HW_coord.append([xxx,yyy,zzz])
                elif atm_type==10:
                   OA_coord.append([xxx,yyy,zzz])
                elif atm_type==9:
                   CA_coord=[xxx,yyy,zzz]
                   nbi+=1
                elif atm_type==4:
                   HH_coord.append([xxx,yyy,zzz])
                elif atm_type==3:
                   OH_coord=[xxx,yyy,zzz]
             
           
            for OW in OW_coord:

                rrr=PBC_dis(OW,CA_coord,boxlength)
                ibin=int(rrr/deltar)
                if ibin<nbin:
                   rdf_CA_OW[ibin]+=1.0

                rrr=PBC_dis(OW,OH_coord,boxlength)
                ibin=int(rrr/deltar)
                if ibin<nbin:
                   rdf_OH_OW[ibin]+=1.0

                for ic in range(0,3):
                    rrr=PBC_dis(OW,OA_coord[ic],boxlength)
                    ibin=int(rrr/deltar)
                    if ibin<nbin:
                       rdf_OA_OW[ibin]+=1.0

                for ic in range(0,3):
                    rrr=PBC_dis(OW,HH_coord[ic],boxlength)
                    ibin=int(rrr/deltar)
                    if ibin<nbin:
                       rdf_HH_OW[ibin]+=1.0
               
            for HW in HW_coord:

                rrr=PBC_dis(HW,CA_coord,boxlength)
                ibin=int(rrr/deltar)
                if ibin<nbin:
                   rdf_CA_HW[ibin]+=1.0

                rrr=PBC_dis(HW,OH_coord,boxlength)
                ibin=int(rrr/deltar)
                if ibin<nbin:
                   rdf_OH_HW[ibin]+=1.0

                for ic in range(0,3):
                    rrr=PBC_dis(HW,OA_coord[ic],boxlength)
                    ibin=int(rrr/deltar)
                    if ibin<nbin:
                       rdf_OA_HW[ibin]+=1.0

                for ic in range(0,3):
                    rrr=PBC_dis(HW,HH_coord[ic],boxlength)
                    ibin=int(rrr/deltar)
                    if ibin<nbin:
                       rdf_HH_HW[ibin]+=1.0
              
            del OW_coord[:]
            del HW_coord[:]
            del OA_coord[:]
            del HH_coord[:]
            CA_coord=[0.0,0.0,0.0]
            OH_coord=[0.0,0.0,0.0]


print "nconf = %d" %nconf
print "nbi = %d" %nbi

g=open(outputname,'w')
print>>g, "# rrr   CA-OW   CA-HW   HH-OW  HH-HW  OA-OW  OA-HW  OH-OW  OH-HW" 
for ibin in range(0,nbin):
    rrr=deltar*(float(ibin)+0.5)
    rdf_CA_OW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoOW*float(nconf)
    rdf_CA_HW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoHW*float(nconf)
    rdf_HH_OW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoOW*float(nconf)*3.0
    rdf_HH_HW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoHW*float(nconf)*3.0
    rdf_OA_OW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoOW*float(nconf)*3.0
    rdf_OA_HW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoHW*float(nconf)*3.0
    rdf_OH_OW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoOW*float(nconf)
    rdf_OH_HW[ibin]/=4.0*math.pi*rrr*rrr*deltar*rhoHW*float(nconf)
    print>>g, "%.6f   %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f" %(rrr,rdf_CA_OW[ibin], \
              rdf_CA_HW[ibin],rdf_HH_OW[ibin],rdf_HH_HW[ibin],rdf_OA_OW[ibin],rdf_OA_HW[ibin], \
              rdf_OH_OW[ibin],rdf_OH_HW[ibin])

g.close()

     
