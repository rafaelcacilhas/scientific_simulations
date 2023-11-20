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
real*8,parameter      ::   E        = -0.8d0
real*8,parameter      ::   D        =  1.0d0
real*8,parameter      ::   a        =  1.0d0
real*8,parameter      ::   m        =  1.0d0
real*8,parameter      ::   re       =  0.0d0
real*8,parameter      ::   r0       =  0.0d0
real*8,parameter      ::   t0       =  0.0d0
real*8,parameter      ::   tn       =  100.0d0
real*8,parameter      ::   dt       =  0.01d0


!Global Variables

real*8		:: v0,U
                            


end module variables

subroutine verlet_original(r0,v0,t0,E,D,a,m,re,tn,dt)

real*8      ::  r0,t0,v0,E,D,a,m,re,tn,dt
real*8      ::  r,v,t,r1,temp,w0,phi,s,exata,erro,Ec,Ep,Et

w0 = dsqrt(-2.0d0*E*a*a/m)
phi =-dacos( (  (E+D)/D    )**0.5d0    )

r = r0
v = v0
t = t0
s = -(D/E)*(   1.0d0 - dsqrt(  (E+D)/D  )*dcos(w0*t + phi )   )
exata = r0 + dlog(s)/a
erro = dabs(exata - r)

Ec = 0.5d0*m*v*v
Ep = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )     +1
Et = Ec + Ep


write(2,*)t,r,v,erro,Ec,Ep,Et
r1 = r

t = dt
r = r + v*dt
v  = v + a*dt

s = -(D/E)*(   1.0d0 - dsqrt(  (E+D)/D  )*dcos(w0*t + phi )   )
exata = r0 + dlog(s)/a
erro = dabs(exata - r)

Ec = 0.5d0*m*v*v
Ep = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )  +1    
Et = Ec + Ep


write(2,*)t,r,v,erro,Ec,Ep,Et

t = 2.0d0*dt



do while (t < tn)
F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*(r-re)   )     ))*dexp(    -a*(r-re)   )  

temp = r
r = 2.0d0*r - r1 + (F/m)*dt*dt
v = (r-r1)/(2.0d0*dt)
r1 = temp


s = -(D/E)*(   1.0d0 - dsqrt(  (E+D)/D  )*dcos(w0*t + phi )   )
exata = r0 + dlog(s)/a
erro = dabs(exata - r)

Ec = 0.5d0*m*v*v
Ep = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )     +1
Et = Ec + Ep

write(2,*)t,r,v,erro,Ec,Ep,Et


t = t + dt
end do



return
end subroutine verlet_original




subroutine verlet_velocity(r0,v0,t0,E,D,a,m,re,tn,dt)

real*8      ::  r0,t0,v0,E,D,a,m,re,tn,dt
real*8      ::  r,v,t,v_meio,w0,phi,s,exata,erro

w0 = dsqrt(-2.0d0*E*a*a/m)
phi =-dacos( (  (E+D)/D    )**0.5d0    )

r = r0
v = v0
t = t0

s = -(D/E)*(   1.0d0 - dsqrt(  (E+D)/D  )*dcos(w0*t + phi )   )
exata = r0 + dlog(s)/a
erro = dabs(exata - r)

Ec = 0.5d0*m*v*v
Ep = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )   +1
Et = Ec + Ep


write(3,*)t,r,v,erro,Ec,Ep,Et

do while (t < tn-dt)                    !Por algum motivo se eu deixar t < tn o gráfico gerado fica em uma escala estranha. Como uma diferença de um ponto não é relevante, preferi mudar o programa.

F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*(r-re)   )     ))*dexp(    -a*(r-re)   )  

v_meio = v + (F/m)*dt/2.0d0
r = r + v_meio*dt
F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*(r-re)   )     ))*dexp(    -a*(r-re)   )  
v = v_meio +  (F/m)*dt/2.0d0

t = t + dt

s = -(D/E)*(   1.0d0 - dsqrt(  (E+D)/D  )*dcos(w0*t + phi )   )
exata = r0 + dlog(s)/a
erro = dabs(exata - r)

Ec = 0.5d0*m*v*v
Ep = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )   +1
Et = Ec + Ep


write(3,*)t,r,v,erro,Ec,Ep,Et

end do



return
end subroutine verlet_velocity

subroutine rougekutta(r0,v0,t0,E,D,a,m,re,tn,dt)

real*8      ::  r0,t0,v0,E,D,a,m,re,tn,dt
real*8      ::  t,F,r,k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v,r1,velocidade,velocidade1,w0,phi,s,exata,erro

w0 = dsqrt(-2.0d0*E*a*a/m)
phi =-dacos( (  (E+D)/D    )**0.5d0    )

r = r0
velocidade = v0
t = t0


do while (t < tn)

    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*(r-re)   )     ))*dexp(    -a*(r-re)   )  
    k1x = dt*velocidade
    k1v = dt*F/m

    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*((r+0.5d0*k1x)-re)   )     ))*dexp(    -a*((r+0.5d0*k1x)-re)   )  
    k2x = dt*(velocidade + 0.5d0*k1v)
    k2v = dt*(F/m)

    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*((r+0.5d0*k2x)-re)   )     ))*dexp(    -a*((r+0.5d0*k2x)-re)   )  
    k3x = dt*(velocidade + 0.5d0*k2v)
    k3v = dt*(F/m)

    F =  -2.0d0*a*D*(              ( 1.0d0 -   dexp(    -a*((r+1.0d0*k3x)-re)   )     ))*dexp(    -a*((r+1.0d0*k3x)-re)   )  
    k4x = dt*(velocidade + 1.0d0*k3v)
    k4v = dt*(F/m)


    r1 = r + (1.0/6.0d0)*(   k1x+ 2.0d0*k2x + 2.0d0*k3x + k4x   ) 
    velocidade1 = velocidade + (1.0/6.0d0)*(   k1v+ 2.0d0*k2v + 2.0d0*k3v + k4v   )


    s = -(D/E)*(   1.0d0 - dsqrt(  (E+D)/D  )*dcos(w0*t + phi )   )
    exata = r0 + dlog(s)/a
    erro = dabs(exata - r)

    Ec = 0.5d0*m*velocidade*velocidade
    Ep = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )   +1 
    Et = Ec + Ep


    write(1,*)t,r,velocidade,erro,Ec,Ep,Et

    r = r1
    velocidade = velocidade1


    t = t + dt
end do

return
end subroutine rougekutta

program morse
use variables

open(1,file="rk.dat")
open(2,file="verletO.dat")
open(3,file="verletV.dat")

r = r0
U = -        D*(    1.0d0   - ( 1.0d0 -   dexp(    -a*(r-re)   )     )**2       )       
v0 = -1.0d0*dsqrt(   (2.0d0/m)*dabs(E-U) )



call rougekutta(r0,v0,t0,E,D,a,m,re,tn,dt)
call verlet_original(r0,v0,t0,E,D,a,m,re,tn,dt)
call verlet_velocity(r0,v0,t0,E,D,a,m,re,tn,dt)
close(1)
close(2)
close(3)

end program


