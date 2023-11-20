! This program solves a differential equation using Euler's Method, which has a global error of h
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 04.08.2018
!
! INPUT:   
! 
! OUTPUT:   





module variables


!Parameters
real*8      ::  k =  4300.0d0
real*8      ::  r =  0.03d0
real*8      ::  N =  730d0
real*8      ::  y0 = 198d0
real*8      ::  t0 = 0d0
real*8      ::  tn = 365d0


!Global Variables
real*8		:: t = 0d0,soma = 0d0,wi1,wi,h,y_exato,erro,erro2,aux


end module variables






program pop
use variables

open(3,file="log")
    write(3,*) "Parametros:"
    write(3,*) "y0 = ",y0
    write(3,*) "k = ",k
    write(3,*) "r = ",r
    write(3,*) "N = ",N
    write(3,*) "t0 = ",t0
    write(3,*) "tn = ",tn
close(3)


open(1,file="730.dat")
dt = (tn - t0) / N

wi = y0


do while (t < tn)

aux = -r*t
y_exato = k*y0/(y0 + (k-y0)*dexp(aux)  )

wi1 = wi + dt* (    r*wi* ( 1.0d0 - wi/k)   )

erro = y_exato - wi
erro = dabs(erro)



write(1,*)t,y_exato,wi,erro


wi = wi1
t = t + dt
end do




end program




