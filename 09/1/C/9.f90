! This program 
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 02.10.2018
!
! INPUT: 
! 
! OUTPUT: 





module variables

!Parameters
integer*8,parameter     ::   N        =  200
real*8,parameter        ::   dH       =  1.0d0
real*8,parameter        ::   mp       =  10.0d0
real*8,parameter        ::   sigma    =  1.0d0
real*8,parameter        ::   L        =  28.0d0
real*8,parameter        ::   epslon   =  0.16021773             !Eletron-volts
real*8,parameter        ::   rc       =  (  2.0d0**(1.0d0/6.0d0)  )*sigma
real*8,parameter        ::   tn       =  10e3
real*8,parameter        ::   dt       =  0.001d0


!Global Variables

real*8		:: x(1:N),y(1:N),vx(1:N),vy(1:N),vi,vf,v2(1:N)
integer*8   :: i,Hn


end module variables

subroutine histograma(v2,N,vi,vf,dH,Hn)

integer*8   ::  N,Hn
real*8      ::  v2(1:N),vi,vf,dH

real*8      ::  x
integer*8   ::  i,k,H(0:Hn)
v2 = dsqrt(v2)
H = 0 

do i = 1,N
    k = floor(  (   v2(i) - vi    ) / dH  )
    H(k) = H(k) + 1
end do

do i = 1,Hn
    write(2,*) i,H(i)
end do








return
end subroutine histograma

program lista
use variables

open(1,file="config.dat")
open(2,file="histograma.dat")
open(3,file="mapa.dat")

do i = 1,N
    read(1,*)x(i),y(i),vx(i),vy(i)
    v2(i) = vx(i)*vx(i) + vy(i)*vy(i)
    write(3,*)x(i),y(i)
end do

vi =  0.0d0
vf =  30.0d0
Hn = (vf - vi) / dH

call histograma(v2,N,vi,vf,dH,Hn)








close(1)
close(2)
end program


