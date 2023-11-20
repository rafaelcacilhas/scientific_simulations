! This program solves three coupled differential equations
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 04.08.2018
!
! INPUT: Equation parameters
! 
! OUTPUT: Solution graphs





module variables



!Constants
real*8      ::  e       = 2.718281828459d0  


!Parameters
real*8      ::  a =  0.2d0
real*8      ::  b =  0.8d0
real*8      ::  m =  0.01d0
integer*8   ::  I0 = 10
integer*8   ::  R0 = 0
integer*8   ::  Nt = 40e3
integer*8   ::  t0 = 0
integer*8   ::  tn = 700
integer*8   ::  N = 7000


!Global Variables

real*8		:: t = 0,wi1(1:3),wi(1:3),h,y_exato,erro,w,fwi,k1(1:3),k2(1:3),k3(1:3),k4(1:3),S0,b1,S_exato,I_exato,R_exato


end module variables






program rougekutta
use variables


open(3,file="log")
    write(3,*) "Parametros:"
    write(3,*) "a = ",a
    write(3,*) "b = ",b
    write(3,*) "m = ",m
    write(3,*) "I0 = ",I0
    write(3,*) "R0 = ",R0
    write(3,*) "Nt = ",Nt
    write(3,*) "t0 = ",t0
    write(3,*) "tn = ",tn
    write(3,*) "N = ",N
close(3)


S0 = Nt - I0 - R0
b1 = b / (Nt**2)
S_exato = Nt*dsqrt(a/b)
I_exato = Nt*(1.0d0-dsqrt(a/b))/(1.0d0 + (a/m))
R_exato = Nt*(1.0d0-dsqrt(a/b))/(1.0d0 + (m/a))




open(1,file="sirs.dat")
dt = real(tn - t0)/N
wi(1) = S0
wi(2) = I0
wi(3) = R0


do while (t < tn)



    k1(1) = m*wi(3)              -  b1*wi(1)*wi(1)*wi(2)
    k1(2) = b1*wi(1)*wi(1)*wi(2) -  a*wi(2)
    k1(3) = a*wi(2)              -  m*wi(3)
    k1 = dt*k1

    k2(1) = m   *(  wi(3) + 0.5d0*k1(3)    ) - b1 *( wi(1)+ 0.5d0*k1(1) ) * (wi(1)+ 0.5d0*k1(1) ) * (   wi(2)+ 0.5d0*k1(2)   )
    k2(2) =-a   *(  wi(2) + 0.5d0*k1(2)    ) + b1 *( wi(1)+ 0.5d0*k1(1) ) * (wi(1)+ 0.5d0*k1(1) ) * (   wi(2)+ 0.5d0*k1(2)   ) 
    k2(3) =-m   *(  wi(3) + 0.5d0*k1(3)    ) + a  *(  wi(2) + 0.5d0*k1(2)        ) 
    k2 = dt*k2

    k3(1) = m   *(  wi(3) + 0.5d0*k2(3)    ) - b1 *( wi(1)+ 0.5d0*k2(1) ) * (wi(1)+ 0.5d0*k2(1) ) * (   wi(2)+ 0.5d0*k2(2)   )
    k3(2) =-a   *(  wi(2) + 0.5d0*k2(2)    ) + b1 *( wi(1)+ 0.5d0*k2(1) ) * (wi(1)+ 0.5d0*k2(1) ) * (   wi(2)+ 0.5d0*k2(2)   )
    k3(3) =-m   *(  wi(3) + 0.5d0*k2(3)    ) + a  *( wi(2)+ 0.5d0*k2(2)          ) 
    k3 = dt*k3

    k4(1) = m   *(  wi(3) + 1.0d0*k3(3)    ) - b1 *( wi(1)+ 1.0d0*k3(1) ) * (wi(1)+ 1.0d0*k3(1) ) * (   wi(2)+ 1.0d0*k3(2)   )
    k4(2) =-a   *(  wi(2) + 1.0d0*k3(2)    ) + b1 *( wi(1)+ 1.0d0*k3(1) ) * (wi(1)+ 1.0d0*k3(1) ) * (   wi(2)+ 1.0d0*k3(2)   )
    k4(3) =-m   *(  wi(3) + 1.0d0*k3(3)    ) + a  *( wi(2)+ 1.0d0*k3(2)          ) 
    k4 = dt*k4

    do i =1,3
        wi1(i) = wi(i) + (1.0/6.0d0)*(k1(i)+ 2.0d0*k2(i) + 2.0d0*k3(i) + k4(i)) 
    end do


    write(1,*)t,wi(1),wi(2),wi(3)


    wi = wi1
    t = t + dt
end do

print*,S_exato,I_exato,R_exato


end program


