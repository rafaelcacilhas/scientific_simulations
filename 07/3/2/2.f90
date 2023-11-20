! This program calculates a fourier series
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 02.09.2018
!
! INPUT:   
! 
! OUTPUT:   


module variables


!Parameters
real*8,parameter      ::  pi      = 4.0d0*datan(1.0d0) 
real*8,parameter      ::  T       = 3.255208333333333E-004
integer*8,parameter   ::  kmax    = 1024


!Global Variables

real*8		:: fk,tk,F,w0
complex(8)  :: soma,expoente,fator,Fn,complexo
integer*8   :: i,j,k,n

end module variables





program fourier
use variables

open(9,file="log")
    write(3,*) "Parametros:"
    write(3,*) "T = ",T
    write(3,*) "A = ",A
    write(3,*) "B = ",B
    write(3,*) "t1 = ",t1
    write(3,*) "dt = ",dt
    write(3,*) "M = ",M
close(9)


open(22,file="saida2.dat")



do n = 0,100


    open(2,file="serie2.dat")
    Fn=0d0
	do k=0,kmax-1
       read(2,*)tk,fk
       complexo = cmplx(    dcos(-2d0*pi*n*k/kmax), dsin(-2d0*pi*n*k/kmax)  )
	   Fn= Fn +complexo*fk
	end do
    close(2)


    write(22,*)n,cdabs(Fn)
end do


close(22)


end program

	
