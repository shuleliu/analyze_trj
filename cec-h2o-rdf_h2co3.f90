program cec_h2o_rdf
!the python script for computing cec-water rdf is too slow
!we rewrite it to fortran code.

integer :: i,j,k
integer :: nwindow
integer, parameter :: nwater=486
integer, parameter :: natom=1467

real*8, parameter :: boxlength=24.3794d0
real*8, parameter :: pi=3.141592653589793

real*8 :: deltar
real*8 :: rwindow,rrr
integer :: nbin
integer :: ncec,nconf
integer :: trueconf
integer :: ibin
character*10 :: iname
character*50 :: trjname,cecname,rdfname
character*80 :: tline

integer :: iatom,imol,atm_type
integer :: istep
real*8 :: x,y,z,q
real*8, allocatable :: rdf_OW(:),rdf_HW(:)

integer,allocatable :: iconf_cec(:)
real*8, allocatable :: cec_coord(:,:)
real*8 :: rhoOW,rhoHW
real*8 :: r_cec(3)
real*8 :: PBC_dis

integer :: narg
character(len=50),allocatable :: args(:)

logical :: match

narg=iargc()
allocate(args(narg))
do i=1,narg
   call getarg(i,args(i))
end do

read(args(1),*) nwindow
read(args(2),*) deltar
read(args(3),*) nbin
read(args(4),*) ncec   ! number of data points in CEC file
read(args(5),*) nconf  ! number of configs in the trj file

allocate(iconf_cec(ncec))
allocate(cec_coord(ncec,3))
allocate(rdf_OW(nbin))
allocate(rdf_HW(nbin))


open(unit=20,file='cec-ow_fortran.dat',status='replace')
open(unit=21,file='cec-hw_fortran.dat',status='replace')

write(20,*)  " rwindow  rrr   g(r)"
write(21,*)  " rwindow  rrr   g(r)"

rhoOW=dble(nwater)/(boxlength**3.d0)
rhoHW=rhoOW*2.d0

do i=1,nwindow
   write(*,*) "Processing window ",i
   cec_coord=0.d0
   iconf_cec=0
   write(iname,"(i4)") i
   trjname="wndow_"//trim(adjustl(iname))//"/bin_"//trim(adjustl(iname))//".0.lammpstrj" 
   cecname="wndow_"//trim(adjustl(iname))//"/cec_coor_"//trim(adjustl(iname))//".dat"
   rwindow=0.25d0*dble(i)
   open(unit=15,file=cecname,action='read')
   read(15,*) tline
   do j=1,ncec
      read(15,*)  iconf_cec(j),cec_coord(j,1:3)   
   end do
   close(15)

   write(*,*) "Finish reading the CEC coordinates" 
   
   rdf_OW=0.d0
   rdf_HW=0.d0
   trueconf=0
   open(unit=16,file=trjname,action='read')
   do j=1,nconf
      match=.false.
      r_cec=0.d0
      read(16,*) tline
      read(16,*) istep

      if(mod(istep,10000).eq.0) then
         write(*,*) "Reading step ",istep
      end if

      do k=1,ncec
         if(iconf_cec(k).eq.istep) then
            match=.true.
            trueconf=trueconf+1
            r_cec=cec_coord(k,:)
            exit
         else if(iconf_cec(k).eq.(istep-1)) then
            match=.true.
            trueconf=trueconf+1
            r_cec=cec_coord(k,:)
            exit
         end if
      end do
!      if(any(iconf_cec.eq.istep)) then
!         match=.true.
!         trueconf=trueconf+1
      !   where(iconf_cec.eq.istep)
      !         r_cec(1)=cec_coord(1) 
      !   endwhere
      !    forall(k=1:ncec)
      !          where(iconf_cec.eq.istep) r_cec=0.d0
      !    end forall
!      else if(any(iconf_cec.eq.(istep-1))) then
!         match=.true.
!         trueconf=trueconf+1
!         where(iconf_cec.eq.(istep-1))
!               r_cec=cec_coord(:) 
!         endwhere
!      end if
      do k=1,7
         read(16,*) tline 
      end do
      do k=1,natom
         read(16,*) iatom,imol,atm_type,q,x,y,z

         if((atm_type.eq.1).and.match) then
            rrr=PBC_dis(x,y,z,r_cec,boxlength)
            ibin=int(rrr/deltar)+1
            if(ibin.lt.nbin) then
               rdf_OW(ibin)=rdf_OW(ibin)+1.0
            end if
         else if((atm_type.eq.2).and.match) then
            rrr=PBC_dis(x,y,z,r_cec,boxlength)
            ibin=int(rrr/deltar)+1
            if(ibin.lt.nbin) then
               rdf_HW(ibin)=rdf_HW(ibin)+1.0
            end if
         end if

      end do

   end do 
   close(16)
 
    
   rdfname="wndow_"//trim(adjustl(iname))//"/cec-ow-hw_fortran.dat"
   open(unit=17,file=rdfname,status='replace')
   write(17,*) "#  rrr   ow   hw"
   do ibin=1,nbin
      rrr=deltar*(dble(ibin)-0.5d0)
      rdf_OW(ibin)=rdf_OW(ibin)/(4.d0*pi*rrr*rrr*deltar*rhoOW*dble(trueconf))
      rdf_HW(ibin)=rdf_HW(ibin)/(4.d0*pi*rrr*rrr*deltar*rhoHW*dble(trueconf))
      write(20,'(3f16.8)')  rwindow,rrr,rdf_OW(ibin)
      write(21,'(3f16.8)')  rwindow,rrr,rdf_HW(ibin)
      write(17,'(3f16.8)')  rrr,rdf_OW(ibin),rdf_HW(ibin)
   end do
   close(17)
   write(20,*)
   write(21,*) 
end do

close(20)
close(21)
end program cec_h2o_rdf

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
