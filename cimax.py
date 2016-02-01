#!/usr/bin/python

#calculate 2-D distribution of ci_max^2
#by reading the evb.out file

nwindow=49
nequil=100000   # number of steps for equilibration
nbin=200
deltabin=1.0/float(nbin)

pp=open("cisq_rrr_dist.dat",'w')
print>>pp, "#  rrr   ci_max_sq   probability" 
#timestep=0 
for i in range(0,nwindow):
    inputname="wndow_"+str(i)+"/bin_"+str(i)+".0.out"
    ci_name="wndow_"+str(i)+"/ci_max_"+str(i)+".dat"
    cec_name="wndow_"+str(i)+"/cec_coor_"+str(i)+".dat"

    g=open(ci_name,'w')
    h=open(cec_name,'w')
    print>>g, "# iconf    ci_max   ci_max_sq"
    print>>h, "# iconf    cec_x   cec_y  cec_z"

    timestep=0
    ci=[]
    ci_max=0.0
    ci_max_sq=0.0
    cisq=[0.0]*(nbin+1)
    rrr=0.25*float(i)
    try:
        print "Openning %s ..." %inputname
        with open(inputname,'r') as f:
             for line in f:
                 if "TIMESTEP" in line:
                    args=line.split()
                    timestep=int(args[1])
                    if timestep%10000==0:
                       print "Reading step %d" %timestep
                 elif "EIGEN_VECTOR" in line and timestep>nequil:
                  #  if timestep%1000==0:
                  #     print "Analyzing step %d" %timestep
                    nxtline=f.next()
                    args=nxtline.split()
                    for item in args:
                        ci.append(float(item)) 
                    
                    ci_max=max(ci)
                    ci_max_sq=ci_max*ci_max
                    print>>g, "%d  %.6e   %.6e" %(timestep,ci_max,ci_max_sq)
                    ibin=int(ci_max_sq/deltabin)
                    cisq[ibin]+=1.0
                    del ci[:]
                    ci_max=0.0
                    ci_max_sq=0.0
                 elif "CEC_COORDINATE" in line and timestep>nequil:
                    nxtline=f.next()
                    args=nxtline.split()
                    print>>h, "%d  %.6f  %.6f  %.6f" %(timestep,float(args[0]), \
                           float(args[1]),float(args[2]))
                    
        for j in range(0,nbin+1):
            cisq[j]/=float(timestep - nequil)
            print>>pp, "%.6f  %.6f  %.6f" %(rrr,deltabin*float(j),cisq[j])

        print>>pp, ""
       
               
    
    except IOError:
        print "Cannot open %s" %inputname
 
    g.close()
    h.close()

pp.close()
 
