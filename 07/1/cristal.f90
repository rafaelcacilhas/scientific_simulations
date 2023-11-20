! This program calculates atoms positions in a crystal
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
real*8,parameter      ::  x_ini   = -pi
real*8,parameter      ::  x_fin   =  pi
real*8,parameter      ::  dx      = 0.1d0
integer*8,parameter   ::  N       = 85  




!Global Variables

complex(8)  :: soma,exp_complexo
real*8      :: posicaox(1:N),posicaoy(1:N),S,posicao1,posicao2,q(1:2),rl(1:2),rm(1:2),dist(1:2)
real*8      :: x,y,expoente,dy,temp
integer*8   :: l,m, i,j,Nx


end module variables


subroutine modullo(z,tamanho)
complex(8)  ::  z
real*8      ::  tamanho

tamanho = dsqrt(real(z,8)**2 + aimag(z)**2)

end subroutine



program cristal
use variables

Nx = (x_fin - x_ini) / dx
dy = dx

open(1,file="posicoes.dat")
open(2,file="mapa.dat")

soma = 0d0


do i = 1, N
read(1,*)posicaox(i),posicaoy(i)
end do

l = 3
m = 4

x = x_ini
y = x_ini

!Dado x e y, temos:

do i = 1,Nx
x = x_ini
do j = 1,Nx

soma = 0d0

q(1) = x
q(2) = y

    do l = 1,N

    rl(1) = posicaox(l)
    rl(2) = posicaoy(l)

        do m = 1,N


        rm(1) = posicaox(m)
        rm(2) = posicaoy(m)

        dist = rm - rl
        expoente = q(1)*dist(1) + q(2)*dist(2)
        exp_complexo = dcmplx(0.0d0,-1.0d0*expoente)
        soma = soma + cdexp(exp_complexo)

        end do
    end do


call modullo(soma,S)
S = S / N
write(2,*)x,y,S


x = x + dx
end do
y = y + dy
end do






















close(4)

end program











subroutine rascunho

complex(8)  :: z,q
real*8		:: z1,z2,parte_real,parte_im, tamanho

z1 = 3.0d0
z2 = 4.0d0

z = dcmplx(   z1 , z2  )
parte_real = real(z)
parte_im   = aimag(z)
q = cdexp(z)


end subroutine rascunho
