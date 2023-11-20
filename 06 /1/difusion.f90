! This program solves a partial diferencial equation
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 17.08.2018
!
! INPUT:   
! 
! OUTPUT:   


module variables


!Parameters
real*8,parameter      ::  pi      = 4.0d0*datan(1.0d0) 
real*8,parameter      ::   L      = 10.0d0  
real*8,parameter      ::   a      = 1.0d0 
real*8,parameter      ::   kbT    = 4.2109d0        ! T = 305K 
real*8,parameter      ::   n      = 0.7644d0  
real*8,parameter      ::   h      = 0.05d0 
real*8,parameter      ::   b      = 0.1d0 
real*8,parameter      ::   xp     = 0.6d0*L 
real*8,parameter      ::   tmax   = 50.00d0 
real*8,parameter      ::   Cp     = 1.00d0 
integer*8,parameter   ::   Nmax   = L / h 


!Global Variables

real*8		:: r,dt,lambda,C0(0:nmax),C(0:nmax),Cvelho(0:nmax),Salvar,ds
real*8      :: k1t
integer*8   :: i,j,Tempo,nome
character(len=10) :: file_id
character(len=50) :: file_name


end module variables




program difusion
use variables

open(3,file="log")
    write(3,*) "Parametros:"
    write(3,*) "L = ",L
    write(3,*) "a = ",a
    write(3,*) "kbT = ",kbT
    write(3,*) "n = ",n
    write(3,*) "h = ",h
    write(3,*) "b = ",b
    write(3,*) "xp = ",xp
close(3)

open(4,file="inicial.dat")

do i = 0,Nmax
C(i) = Cp*dexp(    -b*(i*h-xp)**2   )
write(4,*)i*h,C(i)
end do

close(4)

Salvar = 5.0d0
ds = tmax 





t = dt


D = kbT/(6.0d0*pi*n*a)              ! * 10^-15 m²/s
lambda = 0.1d0
dt = h*h*lambda/D
Tempo = 2000 + tmax / dt
print*,D

C(0) = C(1)
C(nmax) = C(nmax-1)

!Tempo = 1
do j = 1,Tempo


    Cvelho(1) = C(1)
    C(1) = (1.0d0 - 2.0d0*lambda)*C(1) + lambda*(   C(2) + C(0)    )
    C(0) = C(1)

    do i = 2,Nmax-1
        Cvelho(i) = C(i)                                                                        !Como o elemento i usa o elemento i-1, que já foi atualizado, precisamos de uma variavel temporaria
        C(i) = (1.0d0 - 2.0d0*lambda)*C(i) + lambda*(   C(i+1) + Cvelho(i-1)    )



    end do



        if(t >= Salvar) then

            nome = Salvar
            write(file_id, '(i0)') nome
            file_name = '' // trim(adjustl(file_id)) // '.dat'
            open(2,file = trim(file_name))

            write(2,*)0,C(1)
            do i = 1,Nmax-1
                write(2,*)i*h,C(i)
            end do
            write(2,*)Nmax*h,C(Nmax-1)

            close(2)
            Salvar = Salvar + dS

        end if




    C(0) = C(1)
    C(nmax) = C(nmax-1)



t = t + dt
end do
end program


