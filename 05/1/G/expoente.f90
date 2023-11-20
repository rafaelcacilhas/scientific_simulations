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


real*8      :: xn,xn1,lambda,r,x0,y0,yn,yn1,deltaX,fr,dif,melhor_dif,melhor_r
integer*8   :: i,nome,j,jmax
character(len=10) :: file_id
character(len=50) :: file_name

end module variables




program expoente
use variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! r perto de 0.75 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fr = 0d0
lambda = 0d0
r = 0.700d0
dr = 0.001
melhor_dif = 1.0d0
dif = 1.0d0

jmax = (0.80d0 - r)/dr

open(1,file="data.dat")

do j = 1,jmax

    xn = 0.50d0
    lambda = 0d0
    fr = 0d0
    xn1 = 0d0



    do i = 1,N-1

    xn1 = 4.0d0*r*xn*(1.0d0 - xn)
    fr = dabs(4.0d0*r*(1.0d0 - 2.0d0*xn) )



    if ( i > 10) then                           !!!Tempo de "relaxação"
    lambda = lambda + ( dlog(fr) / N    )

    dif = dabs(lambda)

    if(dif < melhor_dif) then
        melhor_dif = dif
        melhor_r = r
    end if

    end if

    xn = xn1

    end do

r = r + dr
end do

print*,melhor_r,melhor_dif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! r perto de 0.85 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fr = 0d0
lambda = 0d0
r = 0.850d0
dr = 0.001
melhor_dif = 1.0d0
dif = 1.0d0

jmax = (0.875d0 - r)/dr

open(1,file="data.dat")

do j = 1,jmax

    xn = 0.50d0
    lambda = 0d0
    fr = 0d0
    xn1 = 0d0



    do i = 1,N-1

    xn1 = 4.0d0*r*xn*(1.0d0 - xn)
    fr = dabs(4.0d0*r*(1.0d0 - 2.0d0*xn) )



    if ( i > 10) then                           !!!Tempo de "relaxação"
    lambda = lambda + ( dlog(fr) / N    )

    dif = dabs(lambda)

    if(dif < melhor_dif) then
        melhor_dif = dif
        melhor_r = r
    end if

    end if

    xn = xn1

    end do

r = r + dr
end do

print*,melhor_r,melhor_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! r perto de 0.886 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fr = 0d0
lambda = 0d0
r = 0.875d0
dr = 0.001
melhor_dif = 1.0d0
dif = 1.0d0

jmax = (0.90d0 - r)/dr

open(1,file="data.dat")

do j = 1,jmax

    xn = 0.50d0
    lambda = 0d0
    fr = 0d0
    xn1 = 0d0



    do i = 1,N-1

    xn1 = 4.0d0*r*xn*(1.0d0 - xn)
    fr = dabs(4.0d0*r*(1.0d0 - 2.0d0*xn) )



    if ( i > 10) then                           !!!Tempo de "relaxação"
    lambda = lambda + ( dlog(fr) / N    )

    dif = dabs(lambda)

    if(dif < melhor_dif) then
        melhor_dif = dif
        melhor_r = r
    end if

    end if

    xn = xn1

    end do

r = r + dr
end do

print*,melhor_r,melhor_dif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! r perto de 0.896 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fr = 0d0
lambda = 0d0
r = 0.887d0
dr = 0.001
melhor_dif = 1.0d0
dif = 1.0d0

jmax = (0.90d0 - r)/dr

open(1,file="data.dat")

do j = 1,jmax

    xn = 0.50d0
    lambda = 0d0
    fr = 0d0
    xn1 = 0d0



    do i = 1,N-1

    xn1 = 4.0d0*r*xn*(1.0d0 - xn)
    fr = dabs(4.0d0*r*(1.0d0 - 2.0d0*xn) )



    if ( i > 10) then                           !!!Tempo de "relaxação"
    lambda = lambda + ( dlog(fr) / N    )

    dif = dabs(lambda)

    if(dif < melhor_dif) then
        melhor_dif = dif
        melhor_r = r
    end if

    end if

    xn = xn1

    end do

r = r + dr
end do

print*,melhor_r,melhor_dif



end program

    
