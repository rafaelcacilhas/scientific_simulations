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
real*8,parameter   ::  g   = -4.0d0
real*8,parameter   ::  k   = 1.04d0
real*8,parameter   ::  L   = 30d0
real*8,parameter   ::  t0  = 305d0
real*8,parameter   ::  tn1 = 390d0
integer*8,parameter ::  Ntamanho = 99

!Global Variables

real*8		:: dx,u(1:Ntamanho),M(1:Ntamanho,1:Ntamanho),w(1:Ntamanho),Lower(1:Ntamanho,1:Ntamanho)
real*8      :: x,xi,y(1:Ntamanho),funcao,soma,u0,u100,temp,exato,erro
real*8      ::  Upper(1:Ntamanho,1:Ntamanho),teste(1:Ntamanho,1:Ntamanho)
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

open(3,file="log")
    write(3,*) "Parametros:"
    write(3,*) "g = ",g
    write(3,*) "k = ",k
    write(3,*) "L = ",L
    write(3,*) "N = ",N
    write(3,*) "t0 = ",t0
    write(3,*) "tn = ",tn
close(3)



open(1,file="calor.dat")
dx = L/(N+1.0d0)
x = 0d0
u0 = t0
u100 = tn1
Lower = 0d0
Upper = 0d0

funcao = 0d0
w(1) = h*h*funcao - u0
w(99) = h*h*funcao - u100

do linha = 2,98
w(linha) = h*h*funcao
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definir M e LU    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
M = 0d0
M(1,1)   = -2.0d0
M(1,2)   =  1.0d0
M(99,98) = 1.0d0
M(99,99) = -2.0d0

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

        print*,"Erro no elemento (",linha,",",coluna,"):"
        print*,"Decomposição LU: ",teste(linha,coluna),"Matriz M:",M(linha,coluna) 

    end if
end do
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Sistema      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Erro e Resultados !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    write(1,*)0,u0,0
do i = 1,99
    exato = ((tn1-t0)/(Ntamanho+1.0d0))*i + t0
    erro = dabs(exato - u(i))
    write(1,*)(L/(Ntamanho+1.0d0) )*i,u(i),erro


end do
    write(1,*)L,tn1,0

end program


