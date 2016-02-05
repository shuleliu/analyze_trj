#!/usr/bin/python

# compute rdf between H2CO3 (in either protonated 
# or deprotonated state) and water from Lammps trajectory

import math
from math import pi,sqrt
import sys
import numpy as np

def PBC_dis(r0,rrr,boxlength):
    x12=rrr[0]-r0[0]
    y12=rrr[1]-r0[1]
    z12=rrr[2]-r0[2]
    x12=x12-boxlength*float(round(x12/boxlength))
    y12=y12-boxlength*float(round(y12/boxlength))
    z12=z12-boxlength*float(round(z12/boxlength))
    r12=sqrt(x12*x12+y12*y12+z12*z12)

    return r12

trjname=sys.argv[1]
# outputname=sys.argv[2]
natom=1466
boxlength=24.3794
nwater=486
deltar=0.05
nbin=240

#rhoOW=float(nwater)/(boxlength**3)
#rhoHW=rhoOW*2.0

nconf=0

C_type=[5,9]     # carbon atom type
OC_type=[7,11]   # carboxyl oxygen atom type
HC_type=[8,12]   # carboxyl hydrogen atom type
OB_type=[6,10]   # carbonyl oxygen atom type
OW_type=[1]
HW_type=[2]
 
nbi=0 #check if the number of times HCO3 shows up is consistent with nconf

OW_coord=[]
HW_coord=[]
OC_coord=[]
OB_coord=[]
C_coord=[]
HC_coord=[]

nOW=0
nHW=0
nOB=0
nOC=0
nHC=0
nC=0

rOW_OC=[]
rOW_OB=[]
rOW_C=[]
rOW_HC=[]

rHW_OC=[]
rHW_OB=[]
rHW_C=[]
rHW_HC=[]

with open(trjname,'r') as trj:
     for line in trj:
         if "id" in line:
            nconf+=1
            if nconf%1000==0:
               print "Reading config %d" %nconf
#           if nconf<1000:
            for j in range(0,natom):
                nxtline=trj.next()
                args=nxtline.split()
                atm_type=int(args[2])
                xxx=float(args[4])
                yyy=float(args[5])
                zzz=float(args[6])
                if atm_type in OW_type:
                   OW_coord.append([xxx,yyy,zzz])
                elif atm_type in HW_type:
                   HW_coord.append([xxx,yyy,zzz])
                elif atm_type in C_type:
                   C_coord.append([xxx,yyy,zzz])
                elif atm_type in OC_type:
                   OC_coord.append([xxx,yyy,zzz])
                elif atm_type in OB_type:
                   OB_coord.append([xxx,yyy,zzz])
                elif atm_type in HC_type:
                   HC_coord.append([xxx,yyy,zzz])
             
            nOW+=len(OW_coord) 
            nHW+=len(HW_coord) 
            nOB+=len(OB_coord) 
            nOC+=len(OC_coord) 
            nHC+=len(HC_coord) 
            nC+=len(C_coord)
 
            for OW in OW_coord:
                for item in C_coord:
                    rrr = PBC_dis(OW,item,boxlength)
                    rOW_C.append(rrr)
                for item in OC_coord:
                    rrr = PBC_dis(OW,item,boxlength)
                    rOW_OC.append(rrr)
                for item in OB_coord:
                    rrr = PBC_dis(OW,item,boxlength)
                    rOW_OB.append(rrr)
                for item in HC_coord:
                    rrr = PBC_dis(OW,item,boxlength)
                    rOW_HC.append(rrr)

            for HW in HW_coord:
                for item in C_coord:
                    rrr = PBC_dis(HW,item,boxlength)
                    rHW_C.append(rrr)
                for item in OC_coord:
                    rrr = PBC_dis(HW,item,boxlength)
                    rHW_OC.append(rrr)
                for item in OB_coord:
                    rrr = PBC_dis(HW,item,boxlength)
                    rHW_OB.append(rrr)
                for item in HC_coord:
                    rrr = PBC_dis(HW,item,boxlength)
                    rHW_HC.append(rrr)

            OW_coord=[]
            HW_coord=[]
            OC_coord=[]
            OB_coord=[]
            C_coord=[]
            HC_coord=[]


print "nconf = %d" %nconf
#print "nbi = %d" %nbi

#print 'nHW %d' %nHW
#print 'nOW %d' %nOW
#print 'nOC %d' %nOC
#print 'nOB %d' %nOB
#print 'nHC %d' %nHC
#print 'nC %d' %nC

rhoHW=float(nHW)/float(nconf)/(boxlength**3)
rhoOW=float(nOW)/float(nconf)/(boxlength**3)
rhoOC=float(nOC)/float(nconf)/(boxlength**3)
rhoOB=float(nOB)/float(nconf)/(boxlength**3)
rhoHC=float(nHC)/float(nconf)/(boxlength**3)
rhoC=float(nC)/float(nconf)/(boxlength**3)

# set up bins for doing histogram
bin_rdf, deltar = np.linspace(0.2,8.0,num=390,retstep=True)

#print "bin_rdf shape is "
#print  bin_rdf.shape

# set up the prefactor for normalization
prefactor = 4.0*pi*np.square(bin_rdf)*deltar*float(nconf)*boxlength**3

#print "prefactor shape is"
#print prefactor.shape
 
# delete the first element of prefactor such that its dimension is consistent with rdf array
prefactor = np.delete(prefactor,0,0)

#print "prefactor shape is"
#print prefactor.shape
 
rOW_OC=np.array(rOW_OC)
rdf_OW_OC, bin_OW_OC = np.histogram(rOW_OC,bins=bin_rdf)
rdf_OW_OC = np.divide(rdf_OW_OC,prefactor)/rhoOW/rhoOC

rOW_C=np.array(rOW_C)
rdf_OW_C, bin_OW_C = np.histogram(rOW_C,bins=bin_rdf)
rdf_OW_C = np.divide(rdf_OW_C,prefactor)/rhoOW/rhoC

rOW_OB=np.array(rOW_OB)
rdf_OW_OB, bin_OW_OB = np.histogram(rOW_OB,bins=bin_rdf)
rdf_OW_OB = np.divide(rdf_OW_OB,prefactor)/rhoOW/rhoOB

rOW_HC=np.array(rOW_HC)
rdf_OW_HC, bin_OW_HC = np.histogram(rOW_HC,bins=bin_rdf)
rdf_OW_HC = np.divide(rdf_OW_HC,prefactor)/rhoOW/rhoHC

rHW_OC=np.array(rHW_OC)
rdf_HW_OC, bin_HW_OC = np.histogram(rHW_OC,bins=bin_rdf)
rdf_HW_OC = np.divide(rdf_HW_OC,prefactor)/rhoHW/rhoOC

rHW_C=np.array(rHW_C)
rdf_HW_C, bin_HW_C = np.histogram(rHW_C,bins=bin_rdf)
rdf_HW_C = np.divide(rdf_HW_C,prefactor)/rhoHW/rhoC

rHW_OB=np.array(rHW_OB)
rdf_HW_OB, bin_HW_OB = np.histogram(rHW_OB,bins=bin_rdf)
rdf_HW_OB = np.divide(rdf_HW_OB,prefactor)/rhoHW/rhoOB

rHW_HC=np.array(rHW_HC)
rdf_HW_HC, bin_HW_HC = np.histogram(rHW_HC,bins=bin_rdf)
rdf_HW_HC = np.divide(rdf_HW_HC,prefactor)/rhoHW/rhoHC


g=open("OW_rdf.dat",'w')
h=open("HW_rdf.dat",'w')

print>>g, "# rrr OW_OC  OW_OC_Coor OW_OB  OW_OB_Coor  OW_HC  OW_HC_Coor OW_C  OW_C_Coor"
print>>h, "# rrr HW_OC  HW_OC_Coor HW_OB  HW_OB_Coor  HW_HC  HW_HC_Coor HW_C  HW_C_Coor"

coor_OW_OC=0.0
coor_OW_OB=0.0
coor_OW_HC=0.0
coor_OW_C=0.0

coor_HW_OC=0.0
coor_HW_OB=0.0
coor_HW_HC=0.0
coor_HW_C=0.0

for i in range(0,rdf_OW_OC.size):
    volume = 4.0*pi*bin_rdf[i]*bin_rdf[i]*deltar
    coor_OW_OC+=rdf_OW_OC[i]*volume*rhoOW
    coor_OW_OB+=rdf_OW_OB[i]*volume*rhoOW
    coor_OW_HC+=rdf_OW_HC[i]*volume*rhoOW
    coor_OW_C+=rdf_OW_C[i]*volume*rhoOW

    coor_HW_OC+=rdf_HW_OC[i]*volume*rhoHW
    coor_HW_OB+=rdf_HW_OB[i]*volume*rhoHW
    coor_HW_HC+=rdf_HW_HC[i]*volume*rhoHW
    coor_HW_C+=rdf_HW_C[i]*volume*rhoHW

    print>>g,"%.6f   %.6f  %.6f  %.6f  %.6f  %.6f  %.6f %.6f %.6f" %(bin_rdf[i], \
              rdf_OW_OC[i],coor_OW_OC,rdf_OW_OB[i],coor_OW_OB,rdf_OW_HC[i],coor_OW_HC, \
              rdf_OW_C[i],coor_OW_C)
    print>>h,"%.6f   %.6f  %.6f  %.6f  %.6f  %.6f  %.6f %.6f %.6f" %(bin_rdf[i], \
              rdf_HW_OC[i],coor_HW_OC,rdf_HW_OB[i],coor_HW_OB,rdf_HW_HC[i],coor_HW_HC, \
              rdf_HW_C[i],coor_HW_C)
 
g.close()
h.close()

     
