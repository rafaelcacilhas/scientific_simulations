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
real*8,parameter      ::   k        =  1.0d0
real*8,parameter      ::   sigma    =  1.0d0
integer*8,parameter   ::   N        =  128
integer*8,parameter   ::   Ns       =  28
real*8,parameter      ::   m        =  1.0d0
real*8,parameter      ::   kb       =  1.0d0
real*8,parameter      ::   d        =  1.0d0
real*8,parameter      ::   t0       =  0.0d0
real*8,parameter      ::   tn       =  1000.00d0
real*8,parameter      ::   dt       =  0.01d0


!Global Variables

real*8		:: v0,U
                            


end module variables




subroutine verlet_velocity(k,sigma,N,Ns,kb,d,m,t0,tn,dt)

real*8      ::  k,sigma,kb,d,m,t0,tn,dt
integer*8   ::  N,Ns
real*8      ::  U(0:N-1),V(0:N-1),temp,P,t,Ec,Ep,Et
integer*8   ::  i,j,vizinhoE,vizinhoD

open(9,file="inicial1.dat")

do i = 0,N-1
    read(9,*)U(i),V(i)
end do



! F = k ( uk1_ - 2uk + uk1)


do while (t < tn)                    

Ec    = 0d0
Ep    = 0d0
Et    = 0d0
Temp  = 0d0
P     = 0d0

do i = 0,N-1

vizinhoE = i-1
vizinhoD = i+1

if(i == 0)   vizinhoE = N-1
if(i == N-1) vizinhoD = 0



    F = k*( u(vizinhoE) - 2.0d0*u(i)   + u(vizinhoD)   ) 
    v_meio = v(i) + (F/m)*dt/2.0d0
    u(i)   = u(i) + v_meio*dt

    F = k*( u(vizinhoE) - 2.0d0*u(i)   + u(vizinhoD)   ) 
    v(i) = v_meio +  (F/m)*dt/2.0d0


    Ep = Ep + 0.5d0*k*((   u(i)  -   u(vizinhoD)    )**2    )
    Ec = Ec + 0.5d0*m*v(i)*v(i)
    Et = Ec + Ep


    Temp = (2.0d0/(N*kb))*Ec 
    P    = P + m*v(i)

end do
write(3,*)t,Ec,Ep,Et,Temp,P

print*,t
t = t + dt
end do



return
end subroutine verlet_velocity

program mola
use variables

open(3,file="data.dat")



call verlet_velocity(k,sigma,N,Ns,kb,d,m,t0,tn,dt)

close(3)

end program


