program do_fft
implicit none

integer :: i
integer :: n
integer :: pn
integer :: l, err
logical :: d_flag
double precision,parameter :: Eh2eV=27.2114d0
double precision,parameter :: Pi=4.d0*atan(1.d0)
double precision,parameter :: c=137.d0
double precision,parameter :: kappa=1.d-4
double precision :: pf=(4.d0*Pi)/(3.d0*c*kappa)
double precision :: angle, damp, time, ts, omega
double precision,dimension(3) :: r0, rio
double precision,allocatable,dimension(:,:) :: r
double precision,allocatable,dimension(:,:) :: nf
double precision,allocatable,dimension(:) :: wsave

d_flag=.true.
!d_flag=.false.
!damp=2.d2
damp=50.0
open(unit=10,file='dip.dat')
!open(unit=30,file='nf.dat')
n=0
do
 read(10,*,iostat=err)
 if (err.lt.0) exit
 n=n+1
end do
rewind(10)
if (d_flag) then
 pn=2**ceiling(log(dble(10*n))/log(2.d0))
else
 pn=n
end if
allocate(r(3,pn))
allocate(nf(3,n))
allocate(wsave(2*pn+15))
!Initialize work array for FFT
!as long as pn does not change
!only have to call once
call dffti(pn,wsave)
l=pn/2
if (mod(pn,2).ne.0) l=(pn+1)/2
r0=0.d0
r=0.d0
nf=0.d0
do i=1,n
! read(30,*) time, nf(:,i)
 read(10,*) time, r(:,i)
 if (i==1) then
  ts=-time
 else if (i==2) then
  ts=ts+time
 end if
 r(:,i)=r(:,i)-nf(:,i)
 if (d_flag) then
  if (i==1) r0(:)=r(:,i)
  r(:,i)=r(:,i)-r0(:)
  r(:,i)=r(:,i)*exp(-time/damp)
 end if
end do
close(10)
close(30)
call dfftf(pn,r(1,:),wsave)
call dfftf(pn,r(2,:),wsave)
call dfftf(pn,r(3,:),wsave)
open(unit=20,file='fft.dat')
do i=2,l
 rio(:)=r(:,2*i-1)
 omega=2.d0*Pi*(dble(i-1)/(dble(pn-1)*ts))
 write(20,'(f12.6,4es19.10)') Eh2eV*omega, -omega*pf*sum(rio), -omega*pf*rio
end do
close(20)
deallocate(r)
deallocate(wsave)
end program
