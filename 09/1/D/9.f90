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
integer*8,parameter     ::   d        =  2
real*8,parameter        ::   kb       =  1.0d0
real*8,parameter        ::   dH       =  1.0d0
real*8,parameter        ::   mp       =  10.0d0
real*8,parameter        ::   sigma    =  1.0d0
real*8,parameter        ::   L        =  28.0d0*sigma
real*8,parameter        ::   epslon   =  0.16021773             !Eletron-volts
real*8,parameter        ::   rc       =  (  2.0d0**(1.0d0/6.0d0)  )*sigma
real*8,parameter        ::   tn       =  10e-4
real*8,parameter        ::   dt       =  0.001d0


!Global Variables

real*8		:: x(1:N),y(1:N),vx(1:N),vy(1:N),vi,vf,v2(1:N),U,F,Fx(1:N,1:N),Fy(1:N,1:N),r,Q,eta,v_meio,dr,distancia
real*8      ::  theta,Teq,forcas(1:N,1:N),Angulo
integer*8   :: i,Hn,particulai,particulaj


end module variables


subroutine inicial(x,y,vx,vy,N,r,dr,U,F,epslon,sigma,rc)

integer*8   ::  N
real*8      ::  x(1:N),y(1:N),vx(1:N),vy(1:N),r,dr,U,F,epslon,sigma,rc  

open(1,file="config.dat")
open(4,file="forca.dat")

do i = 1,N  
    read(1,*)x(i),y(i),vx(i),vy(i)
end do


r = 0.95d0*sigma 
dr = (1.5d0*sigma - 0.95d0*sigma)/1000.0d0

do while (r < 1.5d0*sigma)

    U =  4.0d0*epslon*(      (sigma/r)**12   -   (sigma/r)**6    +   0.25d0      )
    F = 24.0d0*epslon*(2.0d0*(sigma/r)**13   -   (sigma/r)**7                    )/sigma
    if(r > rc) then
        U = 0d0
        F = 0d0
    end if

    write(4,*)r,U/epslon,100.0*F
    r = r + dr



end do

close(1)
close(4)
end subroutine inicial


subroutine forca(x,y,epslon,sigma,rc,N,Fx,Fy,Angulo)
integer*8       ::  N,particulai,particulaJ,sinalx,sinaly
real*8          ::  x(1:N),y(1:N),distancia,theta,epslon,sigma,rc,Fx(1:N,1:N),Fy(1:N,1:N),dx,dy,forcas(1:N,1:N),Angulo

do particulai = 1,N
do particulaj = 1,N

    if(particulai == particulaj) then
        Angulo = 0d0
        Fx(particulai,particulaj) = 0d0
        Fy(particulai,particulaj) = 0d0
        goto 100
    end if

            dx = (  x(particulai) - x(particulaj)  )
            dy = (  y(particulai) - y(particulaj)  )
            distancia = dsqrt(dx*dx + dy*dy)
            Angulo = datan( ( y(particulaj) - y(particulai) ) / (  x(particulaj) - x(particulai) )   ) 

            if(dx > 0)  sinalx = -1
            if(dx <= 0) sinalx =  1
            if(dy > 0)  sinaly = -1
            if(dy <= 0) sinaly =  1

            if(distancia >= rc) then
                Fx(particulai,particulaj) = 0d0
                Fy(particulai,particulaj) = 0d0
            else if(distancia < rc) then
                Fx(particulai,particulaj) = (2.0d0*(sigma/distancia)**13   -   (sigma/distancia)**7    )/sigma
                Fx(particulai,particulaj) = -24.0d0*epslon*dcos(Angulo)*sinalx*Fx(particulai,particulaj) 

                Fy(particulai,particulaj) = (2.0d0*(sigma/distancia)**13   -   (sigma/distancia)**7    )/sigma
                Fy(particulai,particulaj) = -24.0d0*epslon*dabs(dsin(Angulo))*sinaly*Fy(particulai,particulaj) 
            end if


    100 continue
end do 
end do

end subroutine forca


subroutine verlet1(x,y,vx,vy,eta,Q,dt,mp,N,kb,Teq,particulai)

integer*8   ::  N,i,particulai
real*8      ::  x(1:N),y(1:N),vx(1:N),vy(1:N),eta,Q,dt,mp,kb,Teq,Ec

Ec = 0d0
do i = 1,N
    v2 = vx(i)*vx(i) + vy(i)*vy(i)
    Ec = Ec + mp*v2/2.0d0
end do

eta = eta + (   dt/(2.0d0*Q)    )*( (Ec) - (    (d*N + 1.0d0)   /   2.0d0   )*kb*Teq    )


vx(particulai) = vx(particulai) + (dt/2.0d0)*(  (Fx/mp) - (eta*vx(particulai) )   )  
vy(particulai) = vy(particulai) + (dt/2.0d0)*(  (Fy/mp) - (eta*vy(particulai) )   )  




x(particulai) = x(particulai) + vx(particulai)*dt 
y(particulai) = y(particulai) + vy(particulai)*dt


eta = eta + (   dt/(2.0d0*Q)    )*( (Ec) - (    (d*N + 1.0d0)   /   2.0d0   )*kb*Teq    )

end subroutine






subroutine verlet2(vx,particulai,dt,Fx,mp,eta,Fy,vy,N)
integer*8       ::  particulai,N
real*8          ::  vx(1:N),vy(1:N),dt,Fx,Fy,mp,eta

vx(particulai) = (  vx(particulai) + (dt/2.0d0)*(Fx/mp) ) /(1.0d0 + (dt/2.0d0)*eta )   
vy(particulai) = (  vy(particulai) + (dt/2.0d0)*(Fy/mp) ) /(1.0d0 + (dt/2.0d0)*eta )   

end subroutine



program lista
use variables
Teq = 100.0d0
eta = 0d0
t = 0d0

call inicial(x,y,vx,vy,N,r,dr,U,F,epslon,sigma,rc)

do while (t < tn)


        call forca(x,y,epslon,sigma,rc,N,Fx,Fy,Angulo)

        call verlet1(x,y,vx,vy,eta,Q,dt,mp,N,kb,Teq,particulai)

        call forca(x,y,epslon,sigma,rc,N,Fx,Fy,Angulo)





t = t + dt

end do
end program


