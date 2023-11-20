! This program solves a partial diferencial equation for the wave equation
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 28.08.2018
!
! INPUT:   
! 
! OUTPUT:   


module variables


!Parameters
real*8,parameter      ::  pi      = 4.0d0*datan(1.0d0) 
real*8,parameter      :: v_luz    = 3.0d0*10e7  
real*8,parameter      ::   C      = 0.1d0  
real*8,parameter      ::   L      = 0.3d0  
real*8,parameter      ::   h      = 10.0d0
real*8,parameter      ::   k      = 0.1d0  
real*8,parameter      ::  L_fio   = 200.0d0     
integer*8,parameter   ::  Nmax    = (L_fio / h) 

!Global Variables

real*8		:: r,dx,dy,lambda,x,y,v(0:Nmax),f,f1,g,V_ini(0:Nmax),v1(0:Nmax),corrente(0:Nmax),corrente_ini(0:Nmax)  
real*8      :: k1t,v1_(0:Nmax),temp(0:Nmax),corrente1(0:Nmax),corrente1_(0:Nmax)
integer*8   :: i,j

end module variables



program onda
use variables

open(2,file="t0.dat")
open(3,file="t1.dat")
open(4,file="t2.dat")
open(5,file="t5.dat")

lambda = dsqrt( 1.0d0/(L*C)  )*k/h 


x = 0d0
t = 0d0
do i = 0,Nmax
    V_ini(i) = 110.0d0*dsin(x*pi/200.0d0)
    corrente_ini(i) = 5.5d0*dcos(x*pi/200.0d0)	
    x = x + h
end do

x = 0d0
do i = 0,Nmax

    v(i) = V_ini(i)
    corrente(i) = corrente_ini(i)   
    corrente1(i) = corrente_ini(i)
    write(2,*)x,v(i),corrente(i)
    x = x + h

    v1_(i) = v(i)
end do

corrente1(Nmax) = corrente(Nmax)
write(5,*)0,v(0),corrente(0)
write(4,*)0,v(0),corrente(0)
write(3,*)0,v(0),corrente(0)



!!!!!!!!!! Neste programa v(i) representa a voltagem no tempo atual, t; v1(i) representa no tempo t + dt; v1_(i) representa no tempo t - dt



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!t0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x = h
do i = 1,Nmax-1

    v1(i) = (1.0d0 - lambda**2)*V_ini(i) + ( 0.5d0*lambda**2)*(   V_ini(i+1) + V_ini(i-1) )
    corrente1(i) = (1.0d0 - lambda**2)*corrente_ini(i) + ( 0.5d0*lambda**2)*(   corrente_ini(i+1) + corrente_ini(i-1) )

    write(3,*)x,v1(i),corrente1(i)


    x = x + h
    v1_(i) = v_ini(i)
    v(i) = v1(i)
    corrente1_(i) = corrente_ini(i)
    corrente(i) = corrente1(i)

end do
t = t + k



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!t1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

x = h
do i = 1,Nmax-1

    v1(i)  = 2.0d0*(1.0d0 - lambda**2)*v(i) + ( lambda**2 )*(   v(i+1) + v(i-1) ) -   v1_(i)
    corrente1(i) = 2.0d0*(1.0d0 - lambda**2)*corrente(i) + ( lambda**2 )*( corrente(i+1) + corrente(i-1) ) - corrente1_(i)
    write(4,*)x,v1(i),corrente1(i)




    x = x + h
    v1_(i) = v(i)
    v(i) = v1(i)
    corrente1_(i) = corrente(i)
    corrente(i) = corrente1(i)

end do
t = t + k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!tn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



do j = 1,3000
x = h

    do i = 1,Nmax-1

        v1(i) = 2.0d0*(1.0d0 - lambda**2)*v(i) + ( lambda**2 )*(   v(i+1) + v(i-1) ) -   v1_(i)
        corrente1(i) = 2.0d0*(1.0d0 - lambda**2)*corrente(i) + ( lambda**2 )*( corrente(i+1) + corrente(i-1) ) - corrente1_(i)
        x = x + h
     
    v1_(i) = v(i)
    v(i) = v1(i)
    corrente1_(i) = corrente(i)
    corrente(i) = corrente1(i)
    end do

t = t + k

end do


x = 0d0
do i = 0,Nmax
    write(5,*)x,v1(i),corrente1(i)
    x  =  x + h
end do

write(4,*)nmax*h,v1(nmax),corrente1(nmax)
write(3,*)nmax*h,v1(nmax),corrente1(nmax)



close(2)
end program


