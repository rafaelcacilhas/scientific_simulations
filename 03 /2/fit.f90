! This program solves a linear system to interpolate data points
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 29.08.2018
!
! INPUT:    
! 
! OUTPUT:   





module variables

!Parameters


!Global Variables

real*8      :: X(1:13,1:13),U(1:13,1:13),L(1:13,1:13), a(1:13),y(1:13),x_aux(1:13),u_aux(1:13),l_aux(1:13),beta(1:13),temp = 0d0
real*8      :: teste(1:13,1:13)
integer*4   :: i,j,k,linha,coluna



end module variables

subroutine lu(x,l,u)

real*8      ::  X(1:13,1:13),U(1:13,1:13),L(1:13,1:13),temp
integer*4   ::i,j,linha,coluna,k

u = 0d0
temp =0d0
do coluna=1,13,1
        u(1,coluna) = x(1,coluna) 
end do
l(1,1) = 1.0d0
l(2,1) = x(2,1) / u(1,1) 




do i=2,13,1

    do j=1,13,1

        if (i <= j) then

            do k=1,(i-1)
                temp = temp + l(i,k)*u(k,j)
            end do
           
            u(i,j) = x(i,j) - temp

        end if

temp = 0d0


    if(i >= j+1) then

        do k=1,(j-1)
            temp = temp + l(i,k)*u(k,j)
        end do
        
        l(i,j) = (x(i,j) - temp)/u(j,j)   
    end if     

    if(i < j+1) then
        l(i,j) = 0d0
    end if

    if(j == i) then
        l(i,j) = 1.0d0
    
    end if


temp = 0d0


    end do




end do



end subroutine



program matriz
use variables


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Leitura de dados !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(1,file="planck.dat")

do i=1,13
    read(1,*)x_aux(i),y(i)
end do

close(1)


do coluna = 1,13,1

    do linha = 1,13,1
        x(linha,coluna) = x_aux(linha)**(13-coluna)     
    end do

end do 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Solução do sistema !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call lu(x,l,u)

!!!! y = L(U.a)      ->  U.a = Beta  

beta = 0d0
do linha=1,13
    soma = 0d0
    do coluna=1,linha-1
        soma = soma - L(linha,coluna)*beta(coluna)                      !Note que beta(coluna) = 0, a não ser que ela receba outro valor na linha abaixo
    end do
    beta(linha) = y(linha) + soma
end do


a = 0d0
do linha=1,13
    soma = 0d0
    do coluna=1,linha-1
        soma = soma - U(linha,coluna)*a(coluna)                      !Note que beta(coluna) = 0, a não ser que ela receba outro valor na linha abaixo
    end do
    a(linha) = beta(linha) + soma
end do

print*,a

end program


