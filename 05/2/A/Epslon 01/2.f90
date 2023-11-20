! This program solves 
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 17.08.2018
!
! INPUT:   
! 
! OUTPUT:   





module variables


!Parameters

real*8,parameter      ::  e      = 2.718281828459045235
real*8,parameter      ::  y0     = 1.0d0
real*8,parameter      ::  tau    = 1.0d0
real*8,parameter      ::  epslon = 0.1d0
real*8,parameter      ::  tmax   = 25.0d0
real*8,parameter      ::  x0     = 0.1d0




!Global Variables


real*8      :: xn,xn1,yn,yn1,x_exato,y_exato,t,h,wi,wi1,x,y
integer*8   :: i,j,N
character(len=10) :: file_id
character(len=50) :: file_name

end module variables




program metodo
use variables



! Saída de dados:
write(file_id, '(i0)') nome
file_name = '' // trim(adjustl(file_id)) // '.dat'
open(1,file = "tau2.dat")



yn = y0
wi = 0d0
t = 0d0
h = 2.0d0*tau - epslon
N = tmax/h


do i = 1,N

    g   = -yn  / tau
    yn1 =  yn + h*g
    y_exato = y0*dexp(-t/tau)

    write(1,*)t,yn,y_exato

    t = t + h
    yn = yn1

    

end do


close(1)

end program

