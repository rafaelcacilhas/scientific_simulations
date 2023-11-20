! This program solves 
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 11.08.2018
!
! INPUT:   
! 
! OUTPUT:   





module variables


!Parameters

real*8,parameter      ::  e      = 2.718281828459045235
real*8,parameter      ::  y0     = -3.0d0
real*8,parameter      ::  tmax   = 20.0d0
real*8,parameter      ::  x0     =  2.0d0

real*8,parameter      ::  a      =  1.0d0
real*8,parameter      ::  b      =  1.0d0
real*8,parameter      ::  c      =  4.0d0
real*8,parameter      ::  d      = -2.0d0
real*8,parameter      ::  h      =  0.2d0


!Global Variables


real*8      :: xn,xn1,yn,yn1,x_exato,y_exato,t,x,y
integer*8   :: i,j,N
character(len=10) :: file_id
character(len=50) :: file_name

end module variables




program metodo
use variables

open(3,file="log")
    write(3,*) "Parametros:"
    write(3,*) "y0 = ",y0
    write(3,*) "tmax = ",tmax
    write(3,*) "x0 = ",x0
    write(3,*) "a = ",a
    write(3,*) "b = ",b
    write(3,*) "c = ",c
    write(3,*) "d = ",d
    write(3,*) "h = ",h
close(3)



! Saída de dados:
write(file_id, '(i0)') nome
file_name = '' // trim(adjustl(file_id)) // '.dat'
open(1,file = "matriz.dat")



yn = y0
xn = x0
wi = 0d0
t = 0d0
N = tmax/h


do i = 1,N

    xn1 = xn*(a*h + 1.0d0) + yn*(h*b)
    yn1 = yn*(d*h + 1.0d0) + xn*(h*c)

    x_exato = dexp(2.0d0*t) +       dexp(-3.0d0*t)
    y_exato = dexp(2.0d0*t) - 4.0d0*dexp(-3.0d0*t)

    write(1,*)t,xn,yn,x_exato,y_exato

    t = t + h
    yn = yn1
    xn = xn1
    

end do


close(1)

end program

