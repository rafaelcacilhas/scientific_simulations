! This program solves a differential equation for the motion of a harmonic oscillator
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 17.08.2018
!
! INPUT:   
! 
! OUTPUT:   


module variables


!Parameters
real*8,parameter      ::  pi = 4.0d0*datan(1.0d0) 
real*8,parameter      ::   h = 0.1d0*pi
real*8,parameter      ::   q = 0.5d0  !0.5
real*8,parameter      ::   b = 0.9d0  !0.9
real*8,parameter      ::  w0 = (2.0d0/3.0d0)
real*8,parameter      ::  x0 = 0.0d0
real*8,parameter      ::  v0 = 2.0d0 !2.0
real*8,parameter      ::  t0 = 0.0d0
real*8                ::  tn = 10000.0d0*h


!Global Variables

real*8		:: r,N,dr,V,F,x,velocidade,dt,exata,s,phi,erro,teta,w
real*8      :: k1t,k2t,k3t,k4t,k1w,k2w,k3w,k4w,teta1,w1,ep,ek,et
integer*8   :: contador


end module variables




program oscilador
use variables

open(3,file="log")
    write(3,*) "Parametros:"
    write(3,*) "h = ",h
    write(3,*) "q = ",q
    write(3,*) "b = ",b
    write(3,*) "w0 = ",w0
    write(3,*) "x0 = ",x0
    write(3,*) "v0 = ",v0
    write(3,*) "t0 = ",t0
    write(3,*) "tn = ",tn
close(3)

open(2,file="dinamica.dat")


t = t0
dt = h
teta = x0
w = v0



do while (t < tn)


if( teta > pi) then
teta = teta - 2.0d0*pi
end if

if (teta < -1.0d0*pi) then
teta = teta + 2.0d0*pi
end if



    F   = (     -q*w                -1.0d0*dsin( teta            )    +      b*dcos(w0*   (t)             )          )
    k1t = dt*w
    k1w = dt*F

    F   = (     -q*(w+0.5d0*k1w)    -1.0d0*dsin( teta+0.5d0*k1t  )    +      b*dcos(w0*(t+ 0.5d0*dt)      )          )
    k2t = dt*(w + 0.5d0*k1w)
    k2w = dt*F


    F   = (     -q*(w+0.5d0*k2w)    -1.0d0*dsin(  teta+0.5d0*k2t )    +      b*dcos( w0*(t+ 0.5d0*dt)     )       )
    k3t = dt*(w + 0.5d0*k2w)
    k3w = dt*F


    F   = (     -q*(w+1.0d0*k3w)    -1.0d0*dsin(  teta+1.0d0*k3t )    +      b*dcos( w0*(t + dt)          )         )
    k4t = dt*(w + 1.0d0*k3w)
    k4w = dt*F



    teta1 = teta + (1.0/6.0d0)*(   k1t+ 2.0d0*k2t + 2.0d0*k3t + k4t   ) 
    w1    = w    + (1.0/6.0d0)*(   k1w+ 2.0d0*k2w + 2.0d0*k3w + k4w   )


    ep = -2.0d0*dcos(teta) + 2.0d0
    ek =  w*w
    et = ep + ek
    write(2,*)t,teta,w,ep,ek,et

    teta = teta1
    w = w1

     
    t = t + dt



end do

end program


