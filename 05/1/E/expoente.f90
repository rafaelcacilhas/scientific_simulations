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

real*8,parameter      ::  e = 2.718281828459045235
integer*4,parameter   ::  N = 120

!Global Variables


real*8      :: xn,xn1,lambda,r,x0,y0,yn,yn1,deltaX,fr
integer*8   :: i,nome,j,jmax,contador
character(len=10) :: file_id
character(len=50) :: file_name

end module variables




program expoente
use variables

fr = 0d0
lambda = 0d0
r = 0.91d0
xn = 0.5000d0
yn = 0.5001d0
fr = 0d0
xn1 = 0d0
yn1 = 0d0
deltax = 0d0
contador = 0

open(1,file="data.dat")





    do i = 1,N-1

    xn1 = 4.0d0*r*xn*(1.0d0 - xn)
    fr = dabs(4.0d0*r*(1.0d0 - 2.0d0*xn) )

    yn1 = 4.0d0*r*yn*(1.0d0 - yn)
    fr = dabs(4.0d0*r*(1.0d0 - 2.0d0*xn) )

    if ( i > 10) then                           !!!Tempo de "relaxação"
    lambda = lambda + ( dlog(fr) / N    )
    end if

    deltaX = dabs(xn - yn)

    xn = xn1
    yn = yn1
    write(1,*)i,xn,yn,deltax


    end do







print*,lambda


end program

    
