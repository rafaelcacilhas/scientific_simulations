! This program 
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 03.10.2018
!
! INPUT: 
! 
! OUTPUT: 





module variables

! Unidades em 10^6 km e em 1 ano

!Parameters
real*8,parameter      ::   G        =  6.67408E-20 
real*8,parameter      ::   Mc       =  2.2E14 
real*8,parameter      ::   Ms       =  2.0E30 
real*8,parameter      ::   x0       =  8.80e7
real*8,parameter      ::   y0       =  0.0d0
real*8,parameter      ::   v0x      =  0.0d0
real*8,parameter      ::   v0y      =  -54.62d0
real*8,parameter      ::   t0       =  0.0d0
real*8,parameter      ::   tn       =  200.000d0
real*8,parameter      ::   dt       =  0.001d0


!Global Variables

real*8		:: x,y
                            


end module variables


subroutine verlet_velocity(G,Mc,Ms,x0,y0,v0x,v0y,t0,tn,dt)

real*8                ::   G,Mc,Ms,x0,y0,v0x,v0y,t0,tn,dt,Gms
real*8                ::   x,y,vx,vy,vx_meio,vy_meio,t,d,d2,v2,Ec,Ep,Et,L,escala
real*8,parameter      ::   escalaP  =  1e-6                                              !converte km para 10^6km 
real*8,parameter      ::   escalaV  =  31.563d0                                           !converte km/s para 10^6km/ano
real*8,parameter      ::   escalaT  =  3.1563e7                                           !converte s para ano



x  = x0*escalap
y  = y0*escalap
vx = v0x*escalav
vy = v0y*escalav
t  = t0
Gms = 132977435.859d0



do while (t < tn)                    

    d2 = x*x + y*y
    v2 = vx**2 + vy**2
    Ec = 0.5d0*v2
    Ep = -GMs/dsqrt(d2)
    Et = Ec + Ep

    L = Mc*dabs(x*vy - y*vx)
    write(3,*)t,x,y,vx,vy,Ec,Ep,Et,L




    Fx =  -GMs*x/(d2**1.5d0)
    Fy =  -GMs*y/(d2**1.5d0)


    vx_meio = vx + Fx*dt/2.0d0
    vy_meio = vy + Fy*dt/2.0d0

    x = x + vx_meio*dt
    y = y + vy_meio*dt
    d2 = x*x + y*y

    Fx =  -GMs*x/(d2**1.5d0)
    Fy =  -GMs*y/(d2**1.5d0)

    vx = vx_meio + Fx*dt/2.0d0
    vy = vy_meio + Fy*dt/2.0d0

    t = t + dt


end do

return
end subroutine verlet_velocity


program cometa
use variables

open(3,file="saida.dat")

call verlet_velocity(G,Mc,Ms,x0,y0,v0x,v0y,t0,tn,dt)




close(3)
end program


