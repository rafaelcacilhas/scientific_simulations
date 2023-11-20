! This program calculates makes a delta r versus time graph and calculates its derivative
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 23.07.2018
!
! INPUT: Initial and final times
! 
! OUTPUT: Graph and derivatives





module variables

!Parameters

real*8      ::  e =  2.718281828459d0  
real*8      ::  t0 = 1.0E-7 
real*8      ::  tn = 1.0E-2
real*8      ::  epslon = 10.0E-4
real*8      ::  tau = 1.0E-4
real*8      ::  p = 0.7d0
integer*4   ::  N = 500




!Global Variables
integer*8   :: i
real*8		:: t,deltaR,y,y_1,y_2,y1_,y2_,h,dt                                   


end module variables

subroutine funcao(x,y,y_1,y_2,y1_,y2_,tau,p,e,h)

real*8      ::  x,y,y_1,y_2,y1_,y2_,h,tau,p
real*8      ::  e 
h = x / 10.0d0                                                                  !!!!! h = t*dt

y   = (1.0d0 - e**(  -  (x          /tau)             **p )   )

y_1 = (1.0d0 - e**(  -((x+h)       /tau)          **p )   )
y_2 = (1.0d0 - e**(  -((x+2.0d0*h)  /tau)     **p )   )

y1_ = (1.0d0 - e**(  -((x-h)        /tau)           **p )   )
y2_ = (1.0d0 - e**(  -((x-2.0d0*h)  /tau)     **p )   )


return
end subroutine funcao



program derivada
use variables

open(1,file="graph.dat")


t = t0
dt = (tn/t0)**(1.0d0/N)

do i = 1, N

    call funcao(t,y,y_1,y_2,y1_,y2_,tau,p,e,h)

    deltaR = epslon*(1.0d0 - e**(  -(t/tau)**p )   )
    df5 = epslon*(y2_ - 8.0d0*y1_ + 8.0d0*y_1 - y_2)    / (12.0d0*h)


    df5 = df5*t/deltaR

    write(1,*)t,deltaR,df5



t = dt*t
end do


end program


