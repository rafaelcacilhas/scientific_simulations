! This program solves a second order differential equation
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 06.08.2018
!
! INPUT:   
! 
! OUTPUT:   





module variables

!Parameters
real*8   ::  g   =-4.0d0
real*8   ::  k   = 1.04d0
real*8   ::  L   = 30d0
real*8   ::  N   = 99
real*8   ::  t0  = 305d0
real*8   ::  q   = 1.50d0


integer*8,parameter ::  Ntamanho = 99

!Global Variables

real*8		:: dx,u(1:Ntamanho),M(1:Ntamanho,1:Ntamanho),w(1:Ntamanho),Lower(1:Ntamanho,1:Ntamanho)
real*8      :: x,xi,y(1:Ntamanho),funcao,soma,rascunho(1:Ntamanho),u0,u100,temp,gammaL,gammaH,exata,c1,c2
real*8      ::  Upper(1:Ntamanho,1:Ntamanho),teste(1:Ntamanho,1:Ntamanho),alfaL,alfaH,betaL,betaH
integer*8   :: linha,coluna,i,j

end module variables


subroutine LU(a,L,U,Ntamanho)

integer*8   ::  i,j,k,Ntamanho
real*8      ::  a(1:Ntamanho,1:Ntamanho),L(1:Ntamanho,1:Ntamanho),U(1:Ntamanho,1:Ntamanho),uij,lij,soma,teste(1:Ntamanho,1:Ntamanho)



l = 0d0
u = 0d0

do i = 1,Ntamanho
l(i,i) = 1.0d0
end do


do j = 1,Ntamanho

    u(1,j) = a(1,j)

    do i = 2,j
        soma = 0d0
        do k = 1,i-1
            soma = soma + l(i,k)*u(k,j)

        end do

        U(i,j) = a(i,j) - soma


    end do


    do i = j+1,Ntamanho
        soma = 0d0
        do k = 1,j-1
            soma = soma + L(i,k)*U(k,j)
        end do

        L(i,j) = (  a(i,j) - soma   ) /U(j,j)

    end do  


end do
    
end subroutine



program calor
use variables
open(1,file="calor.dat")
open(3,file="log")
    write(3,*) "Parametros:"
    write(3,*) "g = ",g
    write(3,*) "k = ",k
    write(3,*) "L = ",L
    write(3,*) "N = ",N
    write(3,*) "t0 = ",t0
    write(3,*) "tn = ",tn
close(3)

Lower = 0d0
Upper = 0d0
x = 0d0

dx = L/(N+1.0d0)
u0 = t0


betaL = 0d0
alfaL = 1.0d0
gammaL = t0

alfaH = 0d0
betaH = 1.0d0
gammaH = g

funcao = -q/k



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definir W, M e LU    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

w(1)  = dx*dx*funcao  - u0
w(99) = dx*dx*funcao  - g*dx

do linha = 2,98
w(linha) = dx*dx*funcao
end do


M = 0d0
M(1,1)   = -2.0d0
M(1,2)   =  1.0d0
M(99,98) =  1.0d0
M(99,99) = -1.0d0 


do linha = 2,Ntamanho-1
do coluna = 1,Ntamanho

    if (linha == coluna) then
        M(linha,(coluna-1)) = 1.0d0
        M(linha,(coluna+1)) = 1.0d0
        M(linha,coluna) = -2.0d0
    end if

end do
end do


call LU(M,Lower,Upper,Ntamanho)
teste = matmul(Lower,Upper)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Verificação !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do linha = 1,Ntamanho
do coluna = 1,Ntamanho
    if(teste(linha,coluna) /= M(linha,coluna)   )  then

!        print*,"Erro no elemento (",linha,",",coluna,"):"
!        print*,"Decomposição LU: ",teste(linha,coluna),"Matriz M:",M(linha,coluna) 

    end if
end do
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Sistema      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!L.y = w

y(1) = w(1)
do linha = 2,Ntamanho
    soma = w(linha)
        do coluna = 1,linha-1
            soma = soma - lower(linha,coluna)*y(coluna)
        end do
    y(linha) = soma
end do

!Upper.u = y

do linha = Ntamanho,1,-1
    soma = y(linha)
        do coluna = 99,linha+1,-1
            soma = soma - Upper(linha,coluna)*u(coluna)
        end do
    u(linha) = soma / Upper(linha,linha)
end do





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Saída de dados      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c1 = 30.0d0*q/k
c2 = 305d0

    write(1,*)0,u0,u0
do i = 1,99
    exata = -(q/k)*(L/(Ntamanho+1.0d0) )*i*(L/(Ntamanho+1.0d0) )*i + c1*(L/(Ntamanho+1.0d0) )*i + c2
    write(1,*)(L/(Ntamanho+1.0d0) )*i,u(i),exata
end do


u0 = (gammaL*dx - betaL*u(1)) / (alfaL*dx - betaL)
u100 = (gammaH*dx + betah*u(99)) / (alfaH*dx + betaH)
print*,u0,u100

end program


