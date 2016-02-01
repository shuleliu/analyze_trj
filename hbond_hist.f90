program hbond_hist

!compute the h-bond angle distribution to debug energy drift
!note this is only for the case where bicarbonate is protonated

implicit none

integer :: i,j,k
integer :: nwindow
integer, parameter :: nwater=486
integer, parameter :: natom=1466

real*8, parameter :: boxlength=24.3794d0
real*8, parameter :: pi=3.141592653589793
real*8, parameter :: deltatheta=1.0d0

integer, parameter :: nbintheta=180

real*8 :: deltar
real*8 :: rwindow,rrr
integer :: nbin
integer :: nconf
integer :: ibin
character*10 :: iname
character*50 :: trjname,cecname,rdfname
character*80 :: tline

integer :: iatom,imol,atm_type
integer :: istep
real*8 :: x,y,z,q

real*8 :: PBC_dis
real*8 :: OP(3),H(3)
real*8,allocatable :: OW(:,:)
real*8 :: rmin,rOH,rOO
real*8 :: cos_theta,theta
real*8,allocatable :: theta_hist(:)

integer :: iwater
integer :: narg
character(len=50),allocatable :: args(:)

logical :: match

narg=iargc()
allocate(args(narg))
do i=1,narg
   call getarg(i,args(i))
end do

read(args(1),*) nconf
read(args(2),*) trjname

allocate(OW(nwater,3))
allocate(theta_hist(nbintheta))

theta_hist=0.d0

open(unit=20,file=trjname,action='read')
do i=1,nconf
   OP=0.d0
   H=0.d0
   OW=0.d0
   rOH=0.d0
   iwater=0
   read(20,*) tline
   read(20,*) istep
  
   if(mod(istep,10000).eq.0) then
      write(*,*) "Reading step ",istep
   end if

   do k=1,7
      read(20,*) tline
   end do

   do k=1,natom
      read(20,*) iatom,imol,atm_type,q,x,y,z
      !bicarbonate hydroxyl oxygen
      if(atm_type.eq.7) then
         OP(1)=x
         OP(2)=y
         OP(3)=z
      !bicarbonate oxygen
      else if(atm_type.eq.8) then
         H(1)=x
         H(2)=y
         H(3)=z
      ! water oxygen
      else if(atm_type.eq.1) then
         iwater=iwater+1
         OW(iwater,1)=x 
         OW(iwater,2)=y 
         OW(iwater,3)=z 
      end if
   end do

   rOH = PBC_dis(OP(1),OP(2),OP(3),H,boxlength)
   !find the closest water molecule
   !and compute the H-bond angle
   rmin = boxlength 
   do k=1,nwater
      rrr = PBC_dis(OW(k,1),OW(k,2),OW(k,3),H,boxlength)
      if(rrr.lt.rmin) then
        rmin = rrr
        rOO = PBC_dis(OW(k,1),OW(k,2),OW(k,3),OP,boxlength)
        cos_theta = rOH*rOH + rrr*rrr -rOO*rOO
        cos_theta = cos_theta/(2.d0*rOH*rrr)
        theta = acos(cos_theta)/pi*180.d0
      end if  
   end do
   ibin = int(theta/deltatheta)
   if((ibin.gt.0).and.(ibin.le.nbintheta)) then
       theta_hist(ibin) = theta_hist(ibin) + 1.d0      
   end if
end do

open(unit=21,file="hbond_angle.dat",status='replace')
do i=1,nbintheta
   write(21,'(f10.5, f14.6)') deltatheta*dble(i),theta_hist(i)/dble(nconf)
end do

close(20)
close(21)


end program hbond_hist

function PBC_dis(x,y,z,r_cec,boxlength)
real*8 :: x,y,z
real*8 :: r_cec(3)
real*8 :: boxlength
real*8 :: x12,y12,z12,r12
real*8 :: PBC_dis

x12=r_cec(1)-x
y12=r_cec(2)-y
z12=r_cec(3)-z

x12=x12-boxlength*dble(nint(x12/boxlength))
y12=y12-boxlength*dble(nint(y12/boxlength))
z12=z12-boxlength*dble(nint(z12/boxlength))

r12=sqrt(x12*x12+y12*y12+z12*z12)

PBC_dis=r12

end function
