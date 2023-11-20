! This program calculates a fourier series
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 21.08.2018
!
! INPUT:   
! 
! OUTPUT:   


module variables


!Parameters
real*8,parameter      ::  pi      = 4.0d0*datan(1.0d0) 
real*8,parameter      ::  T       = 2.0d0*pi  
real*8,parameter      ::  A       = 7.0d0 
real*8,parameter      ::  B       = 0.0d0  
real*8,parameter      ::  t1      = 0.0d0  
real*8,parameter      ::  dt      = 0.05d0  
integer*8,parameter   ::  M       = 55


!Global Variables

real*8		:: a0,an(1:M),bn(1:M),w0,wn,soma,t_ini,t_f,h,integral,tempo,funcao,soma_a,soma_b
real*8      :: bn_exato(1:M)
integer*8   :: i,j,k,npassos,n

end module variables






program fourier
use variables

soma = 0d0
open(9,file="log")
    write(3,*) "Parametros:"
    write(3,*) "T = ",T
    write(3,*) "A = ",A
    write(3,*) "B = ",B
    write(3,*) "t1 = ",t1
    write(3,*) "dt = ",dt
    write(3,*) "M = ",M
close(9)


open(3,file="original.dat")
open(2,file="55.dat")

t_ini = -0.5d0*T
t_f   = 0.5d0*T

bn_exato = 0d0
do i = 1,M,2
bn_exato(i) = 4.0d0*A /(i*pi)
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! a0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integral = 0d0
tempo = t_ini                            
npassos = (t_f-t_ini)/dt

do k =1,npassos/2 +2

    if ( tempo < t1 ) then
        y  = -A + B
        y1 = -A + B
        y2 = -A + B

    else if ( tempo >= t1) then

        y  = A + B
        y1 = A + B
        y2 = A + B


    end if

integral = integral + (y + 4.0d0*y1 + y2 )*dt/3.0d0
tempo = tempo + 2.0d0*dt

end do
a0 = (2.0d0/T)*integral


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! an   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do n = 1,M

wn = n*2.0d0*pi/T

integral = 0d0
tempo = t_ini                   


    do k =1,npassos/2 + 2

        if ( tempo < t1 ) then
            y  = (-A + B)*dcos( wn*(tempo           ) )
            y1 = (-A + B)*dcos( wn*(tempo +       dt) )
            y2 = (-A + B)*dcos( wn*(tempo + 2.0d0*dt) )


        else if ( tempo >= t1) then

            y  = (A + B)*dcos( wn* tempo             )
            y1 = (A + B)*dcos( wn*(tempo +       dt) )
            y2 = (A + B)*dcos( wn*(tempo + 2.0d0*dt) )


        end if

    integral = integral + (y + 4*y1 + y2)*dt/3.0d0
    tempo = tempo + 2*dt
    end do


an(n) = (2.0d0/T)*integral

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! bn   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do n = 1,M

wn = n*2.0d0*pi/T

integral = 0d0
tempo = t_ini                   


    do k =1,npassos/2 + 2

        if ( tempo < t1 ) then
            y  = (-A + B)*dsin( wn*(tempo           ) )
            y1 = (-A + B)*dsin( wn*(tempo +       dt) )
            y2 = (-A + B)*dsin( wn*(tempo + 2.0d0*dt) )


        else if ( tempo >= t1) then

            y  = (A + B)*dsin( wn* tempo             )
            y1 = (A + B)*dsin( wn*(tempo +       dt) )
            y2 = (A + B)*dsin( wn*(tempo + 2.0d0*dt) )


        end if

    integral = integral + (y + 4*y1 + y2)*dt/3.0d0
    tempo = tempo + 2*dt
    end do


bn(n) = (2.0d0/T)*integral

end do

tempo = -4.0d0*pi



an = 0d0

do k = 1,4*npassos + 1
soma_a = 0d0
soma_b = 0d0

    do i = 1,M
    wn = i*2.0d0*pi/T
    soma_a = soma_a + an(i)      *dcos(wn*tempo)
    soma_b = soma_b + bn(i)      *dsin(wn*tempo)
    end do
    funcao = 0.5d0*a0 + soma_a + soma_b

    write(2,*)tempo,funcao
tempo = tempo + dt
end do


funcao = 0d0
tempo = -4.0d0*pi

do while (tempo < 4.0d0*pi)

if ( tempo < 0) then
funcao = -A


else if (tempo >= 0) then
funcao = A

end if

write(3,*) tempo,funcao
tempo = tempo + dt


end do

close(2)
end program


