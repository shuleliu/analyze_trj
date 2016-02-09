// rewrite the python code in cpp
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

double PBC_dis(double *r1, double *r2, double boxlength)
{
 double x12,y12,z12;
 double r12;

 x12 = r1[0] - r2[0] - double(round((r1[0] - r2[0])/boxlength))*boxlength;
 y12 = r1[1] - r2[1] - double(round((r1[1] - r2[1])/boxlength))*boxlength;
 z12 = r1[2] - r2[2] - double(round((r1[2] - r2[2])/boxlength))*boxlength;

 r12=sqrt(x12*x12 + y12*y12 + z12*z12);
 
 return r12;
}


int main(int argc, char **argv)
{
  // trj file name
  char* trjname = argv[1];

  int natom=1466;
  double boxlength=24.3794;
  int nwater=486;
  double deltar=0.05;
  int nbin=240;

  int nconf=0;
  int nskip=10000; // number of configurations to skip

  const double pi=3.141592653589793; 
  int ibin;
  double volume,rrr,deltavol;
 
  vector<double> xOW,yOW,zOW;
  vector<double> xHW,yHW,zHW;
  vector<double> xC,yC,zC;
  vector<double> xOC,yOC,zOC;  
  vector<double> xOB,yOB,zOB;  
  vector<double> xHC,yHC,zHC; 

  int atom_type,atom_id,mol_id;
  double q,xxx,yyy,zzz;
  double fx,fy,fz;

  int nOW,nHW,nOC,nOB,nHC,nC;
  int nOW_tot,nHW_tot,nOC_tot,nOB_tot,nHC_tot,nC_tot;

  double *r1,*r2;
  double r12;

  double *rdfOW_OC,*rdfOW_OB,*rdfOW_HC,*rdfOW_C;
  double *rdfHW_OC,*rdfHW_OB,*rdfHW_HC,*rdfHW_C;

  double  coord_OW_OC,coord_OW_OB,coord_OW_HC,coord_OW_C;
  double  coord_HW_OC,coord_HW_OB,coord_HW_HC,coord_HW_C;

  double rhoOW,rhoHW,rhoOC,rhoOB,rhoHC,rhoC;

  printf("Reading trajector file: %s\n", trjname);

  FILE *trjfile;
  trjfile = fopen(trjname,"r");
  if (trjfile==NULL) perror ("Error opening trajectory file");

  nOW_tot=0;
  nHW_tot=0;
  nOC_tot=0;
  nOB_tot=0;
  nHC_tot=0;
  nC_tot=0;

  r1=new double[3];
  r2=new double[3];

  rdfOW_OC = new double[nbin];
  rdfOW_OB = new double[nbin];
  rdfOW_HC = new double[nbin];
  rdfOW_C = new double[nbin];

  rdfHW_OC = new double[nbin];
  rdfHW_OB = new double[nbin];
  rdfHW_HC = new double[nbin];
  rdfHW_C = new double[nbin];

  volume = pow(boxlength,3.0);

  char arr[100];
  while(fgets(arr,100,trjfile)!= NULL)
  {
   if (strstr(arr,"id")!= NULL) 
   {
     nconf=nconf+1;
     if(nconf%1000==0) printf("Reading configuration %d\n", nconf);

     // skip the equilibration steps
     if(nconf>nskip){
    
     for(int k=0;k<natom;k=k+1)
     {
      fscanf(trjfile,"%d %d %d %lf %lf %lf %lf %lf %lf %lf\n",&atom_id,&mol_id,&atom_type,&q,&xxx,&yyy,&zzz,&fx,&fy,&fz);
//      printf("atom_type = %d\n",atom_type);

      if(atom_type==5 || atom_type==9)
      {
        //C
          xC.push_back(xxx);
          yC.push_back(yyy);
          zC.push_back(zzz);
 //         printf("carbon detected");
       }
       else if(atom_type==7 || atom_type==11)
       {
        //OC
          xOC.push_back(xxx);
          yOC.push_back(yyy);
          zOC.push_back(zzz);
       }
       else if(atom_type==8 || atom_type==12)
       {
        //HC
          xHC.push_back(xxx);
          yHC.push_back(yyy);
          zHC.push_back(zzz);
       }
       else if(atom_type==6 || atom_type==10)
       {
        //OB
          xOB.push_back(xxx);
          yOB.push_back(yyy);
          zOB.push_back(zzz);
       }
       else if(atom_type==1)
       {
          xOW.push_back(xxx);
          yOW.push_back(yyy);
          zOW.push_back(zzz);
       }
       else if(atom_type==2)
       {
          xHW.push_back(xxx);
          yHW.push_back(yyy);
          zHW.push_back(zzz);
       }
      /*
      switch(atom_type)
      {
        //C
        case 5:
        case 9:
          xC.push_back(xxx);
          yC.push_back(yyy);
          zC.push_back(zzz);
          printf("carbon detected");
          break;
        //OC
        case 7:
        case 11:
          xOC.push_back(xxx);
          yOC.push_back(yyy);
          zOC.push_back(zzz);
          break;
        //HC
        case 8:
        case 12:
          xHC.push_back(xxx);
          yHC.push_back(yyy);
          zHC.push_back(zzz);
          break;
        //OB
        case 6:
        case 10:
          xOB.push_back(xxx);
          yOB.push_back(yyy);
          zOB.push_back(zzz);
          break;
        case 1:
          xOW.push_back(xxx);
          yOW.push_back(yyy);
          zOW.push_back(zzz);
          break;
        case 2:
          xHW.push_back(xxx);
          yHW.push_back(yyy);
          zHW.push_back(zzz);
          break;
          
      }
      */
     }
 
    // count number of particles
    nOW=xOW.size();
    nHW=xHW.size(); 
    nOC=xOC.size(); 
    nOB=xOB.size(); 
    nHC=xHC.size(); 
    nC=xC.size();

    nOW_tot=nOW_tot + nOW;
    nHW_tot=nHW_tot + nHW;
    nOC_tot=nOC_tot + nOC;
    nOB_tot=nOB_tot + nOB;
    nHC_tot=nHC_tot + nHC;
    nC_tot=nC_tot + nC;
    
    //accumulate rdf statistics
    for(int i=0;i<nOW;i=i+1)
    {
      r1[0]=xOW[i];
      r1[1]=yOW[i];
      r1[2]=zOW[i];
      
      // OW-OC 
      for(int j=0;j<nOC;j=j+1)
      {
       r2[0]=xOC[j];
       r2[1]=yOC[j];
       r2[2]=zOC[j];

       r12 = PBC_dis(r1,r2,boxlength);

//       printf("r1 = %f  %f  %f\n",r1[0],r1[1],r1[2]);
//       printf("r2 = %f  %f  %f\n",r2[0],r2[1],r2[2]);
//       printf("boxlength = %f\n",boxlength);
//       printf("r12 = %f\n",r12);

       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfOW_OC[ibin]+=1.0;
                 
      }

      // OW-OB 
      for(int j=0;j<nOB;j=j+1)
      {
       r2[0]=xOB[j];
       r2[1]=yOB[j];
       r2[2]=zOB[j];

       r12 = PBC_dis(r1,r2,boxlength);
       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfOW_OB[ibin]+=1.0;
                 
      }

      // OW-HC 
      for(int j=0;j<nHC;j=j+1)
      {
       r2[0]=xHC[j];
       r2[1]=yHC[j];
       r2[2]=zHC[j];

       r12 = PBC_dis(r1,r2,boxlength);
       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfOW_HC[ibin]+=1.0;
                 
      }

      // OW-C 
      for(int j=0;j<nC;j=j+1)
      {
       r2[0]=xC[j];
       r2[1]=yC[j];
       r2[2]=zC[j];

       r12 = PBC_dis(r1,r2,boxlength);
       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfOW_C[ibin]+=1.0;
                 
      }

 
    }  
 
    for(int i=0;i<nHW;i=i+1)
    {
      r1[0]=xHW[i];
      r1[1]=yHW[i];
      r1[2]=zHW[i];
      
      // HW-OC 
      for(int j=0;j<nOC;j=j+1)
      {
       r2[0]=xOC[j];
       r2[1]=yOC[j];
       r2[2]=zOC[j];

       r12 = PBC_dis(r1,r2,boxlength);
       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfHW_OC[ibin]+=1.0;
                 
      }

      // HW-OB 
      for(int j=0;j<nOB;j=j+1)
      {
       r2[0]=xOB[j];
       r2[1]=yOB[j];
       r2[2]=zOB[j];

       r12 = PBC_dis(r1,r2,boxlength);
       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfHW_OB[ibin]+=1.0;
                 
      }

      // HW-HC 
      for(int j=0;j<nHC;j=j+1)
      {
       r2[0]=xHC[j];
       r2[1]=yHC[j];
       r2[2]=zHC[j];

       r12 = PBC_dis(r1,r2,boxlength);
       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfHW_HC[ibin]+=1.0;
                 
      }

      // HW-C 
      for(int j=0;j<nC;j=j+1)
      {
       r2[0]=xC[j];
       r2[1]=yC[j];
       r2[2]=zC[j];

       r12 = PBC_dis(r1,r2,boxlength);
       ibin = int(r12/deltar);
       if(ibin<nbin && ibin>=0) rdfHW_C[ibin]+=1.0;
                 
      }

 
    }

    xOW.clear();  
    yOW.clear();  
    zOW.clear();  
 
    xHW.clear();  
    yHW.clear();  
    zHW.clear();  
 
    xOC.clear();  
    yOC.clear();  
    zOC.clear();  
 
    xOB.clear();  
    yOB.clear();  
    zOB.clear();  
 
    xHC.clear();  
    yHC.clear();  
    zHC.clear();  
 
    xC.clear();  
    yC.clear();  
    zC.clear();  
 
   }
    }
  }
  
  fclose(trjfile);

  // subtract steps skipped
  nconf=nconf-nskip;

  rhoOW=double(nOW_tot)/double(nconf)/volume;
  rhoHW=double(nHW_tot)/double(nconf)/volume;
  rhoOC=double(nOC_tot)/double(nconf)/volume;
  rhoOB=double(nOB_tot)/double(nconf)/volume;
  rhoHC=double(nHC_tot)/double(nconf)/volume;
  rhoC=double(nC_tot)/double(nconf)/volume;

/*
  printf("rhoOW = %f\n",rhoOW);
  printf("rhoHW = %f\n",rhoHW);
  printf("rhoOC = %f\n",rhoOC);
  printf("rhoHC = %f\n",rhoHC);
  printf("rhoOB = %f\n",rhoOB);
  printf("rhoC = %f\n",rhoC);
*/

  FILE  *FP,*GP;

  FP = fopen("rdf_OW.dat","w"); 
  GP = fopen("rdf_HW.dat","w");

  fprintf(FP,"# rrr  OW-OC  Coor  OW-OB Coor OW-HC  Coor  OW-C  Coor\n"); 
  fprintf(GP,"# rrr  HW-OC  Coor  HW-OB Coor HW-HC  Coor  HW-C  Coor\n");

  coord_OW_OC=0.0;
  coord_OW_OB=0.0;
  coord_OW_HC=0.0;
  coord_OW_C=0.0;

  coord_HW_OC=0.0;
  coord_HW_OB=0.0;
  coord_HW_HC=0.0;
  coord_HW_C=0.0;

  for(ibin=0;ibin<nbin;ibin++)
  {
   rrr = deltar*(double(ibin)+0.5);
   deltavol = 4.0*pi*rrr*rrr*deltar;
   rdfOW_OC[ibin] = rdfOW_OC[ibin]/(deltavol*double(nconf)*rhoOW*rhoOC*volume);
   coord_OW_OC+=rdfOW_OC[ibin]*deltavol*rhoOW;
 
   rdfOW_OB[ibin] = rdfOW_OB[ibin]/(deltavol*double(nconf)*rhoOW*rhoOB*volume);
   coord_OW_OB+=rdfOW_OB[ibin]*deltavol*rhoOW;
 
   rdfOW_HC[ibin] = rdfOW_HC[ibin]/(deltavol*double(nconf)*rhoOW*rhoHC*volume);
   coord_OW_HC+=rdfOW_HC[ibin]*deltavol*rhoOW;
 
   rdfOW_C[ibin] = rdfOW_C[ibin]/(deltavol*double(nconf)*rhoOW*rhoC*volume);
   coord_OW_C+=rdfOW_C[ibin]*deltavol*rhoOW;

   fprintf(FP, "%f %f %f %f %f %f %f %f %f\n",rrr,rdfOW_OC[ibin],coord_OW_OC,rdfOW_OB[ibin], \
          coord_OW_OB,rdfOW_HC[ibin],coord_OW_HC,rdfOW_C[ibin],coord_OW_C);
 
   rdfHW_OC[ibin] = rdfHW_OC[ibin]/(deltavol*double(nconf)*rhoHW*rhoOC*volume);
   coord_HW_OC+=rdfHW_OC[ibin]*deltavol*rhoHW;
 
   rdfHW_OB[ibin] = rdfHW_OB[ibin]/(deltavol*double(nconf)*rhoHW*rhoOB*volume);
   coord_HW_OB+=rdfHW_OB[ibin]*deltavol*rhoHW;
 
   rdfHW_HC[ibin] = rdfHW_HC[ibin]/(deltavol*double(nconf)*rhoHW*rhoHC*volume);
   coord_HW_HC+=rdfHW_HC[ibin]*deltavol*rhoHW;
 
   rdfHW_C[ibin] = rdfHW_C[ibin]/(deltavol*double(nconf)*rhoHW*rhoC*volume);
   coord_HW_C+=rdfHW_C[ibin]*deltavol*rhoHW;

   fprintf(GP, "%f %f %f %f %f %f %f %f %f\n",rrr,rdfHW_OC[ibin],coord_HW_OC,rdfHW_OB[ibin], \
          coord_HW_OB,rdfHW_HC[ibin],coord_HW_HC,rdfHW_C[ibin],coord_HW_C);
 
  }
  fclose(FP);
  fclose(GP);

  delete [] rdfOW_OC; 
  delete [] rdfOW_OB; 
  delete [] rdfOW_HC; 
  delete [] rdfOW_C;
 
  delete [] rdfHW_OC; 
  delete [] rdfHW_OB; 
  delete [] rdfHW_HC; 
  delete [] rdfHW_C;

  return 0; 
}
   
      
   
