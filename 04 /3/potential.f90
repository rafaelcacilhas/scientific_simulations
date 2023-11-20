! This program solves a differential equation for the motion of a particle
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 05.08.2018
!
! INPUT:   
! 
! OUTPUT:   





module variables


!Parameters
real*8,parameter      ::   Energy = -0.8d0
real*8,parameter      ::   D =  1.0d0
real*8,parameter      ::   a =  1.0d0
real*8,parameter      ::   m = 1.0d0
real*8,parameter      ::  re = 0.0d0
real*8,parameter      ::  r0 = 0.0d0
real*8,parameter      ::  rf = 5.0d0
real*8,parameter      ::  t0 = 0.0d0
real*8,parameter      ::  tn = 10.0d0


!Global Variables

real*8		:: r,N,dr,V,F,x,velocidade,dt,exata,s,w0,phi,erro
real*8      :: k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v,x1,velocidade1


end module variables




program potential
use variables
open(1,file="VeF.dat")
open(2,file="dinamica.dat")

r = -1.0d0
t = 0d0
N = 300
dr = (rf - r0)/N
N = 2000          
dt = (tn - t0)/N



x = r0
V = -  D*(    1.0d0   -       ( 1.0d0 -   dexp(    -a*(x-re)   )     )**2       )       
velocidade = -dsqrt(2.0d0*dabs(Energy-V)/m )




do while (r < rf)
    V = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )       
    F = -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*(r-re)   )     ))*dexp(    -a*(r-re)   )            

    write(1,*)r,V,F

    r = r + dr
end do


w0 = dsqrt(-2.0d0*Energy*a*a/m)
phi =-dacos( (  (Energy+D)/D    )**0.5d0    )

t = 0d0


do while (t < tn)


    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*(x-re)   )     ))*dexp(    -a*(x-re)   )  
    k1x = dt*velocidade
    k1v = dt*F/m

    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*((x+0.5d0*k1x)-re)   )     ))*dexp(    -a*((x+0.5d0*k1x)-re)   )  
    k2x = dt*(velocidade + 0.5d0*k1x)
    k2v = dt*(F/m)

    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*((x+0.5d0*k2x)-re)   )     ))*dexp(    -a*((x+0.5d0*k2x)-re)   )  
    k3x = dt*(velocidade + 0.5d0*k2x)
    k3v = dt*(F/m)

    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*((x+1.0d0*k3x)-re)   )     ))*dexp(    -a*((x+1.0d0*k3x)-re)   )  
    k4x = dt*(velocidade + 1.0d0*k3x)
    k4v = dt*(F/m)


    x1 = x + (1.0/6.0d0)*(   k1x+ 2.0d0*k2x + 2.0d0*k3x + k4x   ) 
    velocidade1 = velocidade + (1.0/6.0d0)*(   k1v+ 2.0d0*k2v + 2.0d0*k3v + k4v   )





s = -(D/Energy)*(   1.0d0 - dsqrt(  (Energy+D)/D  )*dcos(w0*t + phi )   )
exata = r0 + dlog(s)/a
erro = dabs(exata - x)
    

    write(2,*)t,x,velocidade,erro


    x = x1
    velocidade = velocidade1


    t = t + dt
end do



end program


