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
real*8,parameter        ::   kb       =  1.0d0
real*8,parameter        ::   dH       =  1.0d0
real*8,parameter        ::   mp       =  10.0d0
real*8,parameter        ::   sigma    =  1.0d0
real*8,parameter        ::   L        =  28.0d0
real*8,parameter        ::   epslon   =  0.16021773             !Eletron-volts
real*8,parameter        ::   rc       =  (  2.0d0**(1.0d0/6.0d0)  )*sigma
real*8,parameter        ::   tn       =  10e3
real*8,parameter        ::   dt       =  0.001d0


!Global Variables

real*8		:: x(1:N),y(1:N),vx(1:N),vy(1:N),vi,vf,v2(1:N),U,F,r
integer*8   :: i,Hn


end module variables


program lista
use variables

open(1,file="config.dat")
open(4,file="forca.dat")

do i = 1,N
    read(1,*)x(i),y(i),vx(i),vy(i)
end do


r = 0.95d0*sigma 
dr = (1.5d0*sigma - 0.95d0*sigma)/1000.0d0
print*,rc
r = 0.95d0   
do while (r < 1.5d0*sigma)

    U =  4.0d0*epslon*(      (sigma/r)**12   -   (sigma/r)**6    +   0.25d0      )
    F = 24.0d0*epslon*(2.0d0*(sigma/r)**13   -   (sigma/r)**7                    )/sigma
    if(r > rc) then
        U = 0d0
        F = 0d0
    end if

    write(4,*)r,U/epslon,100*F
    r = r + dr

end do

close(1)
close(2)
end program


