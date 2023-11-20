! This program calculates the coefficients to a Lorentz fit
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 30.08.2018
!
! INPUT:    
! 
! OUTPUT:   





module variables

!Parameters

real*8      ::  w0       = 5.0d0;
real*8      ::  gama2    = 1.0d0;
real*8      ::  A        = 16.0d0;

!Global Variables

real*8      :: x(1:16),y(1:16),Jac(1:3,1:3),yq, z,qui, mi, delta,H(3,3),q(1:3),f(1:3),U(1:3,1:3),L(1:3,1:3),yes(1:3,1:3)
real*8      :: funcao,y_sub,x_sub,df1(1:3),df2(1:3),df3(1:3),JacI(1:3,1:3)
integer*4   :: i,j,linha,coluna


end module variables

subroutine lu(x,l,u)

real*8      ::  X(1:3,1:3),U(1:3,1:3),L(1:3,1:3),temp
integer*4   ::i,j,linha,coluna,k,N

N = 3               !!! Tamanho da matriz
u = 0d0
temp =0d0
do coluna=1,N,1
        u(1,coluna) = x(1,coluna) 
end do
l(1,1) = 1.0d0
l(2,1) = x(2,1) / u(1,1) 




do i=2,N,1

    do j=1,N,1

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
open(1,file="lorentziana.dat")
do i=1,16
    read(1,*)x(i),y(i)
end do

!Chutes iniciais
q(1) = A                
q(2) = w0                 
q(3) = gama2       



jac = 0.d0    
chi = 0.d0    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Erro   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	do i=1,16	
        yq = q(1) / (      (   x(i) - q(2)  )**2 + q(3)  ) 
        delta = yq  -   y(i)   

		qui = qui + delta**2.d0	
	enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Matriz Jacobiana !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i = 1,1
    !!Vamos definir duas variaveis auxiliares, porque elas vão aparecer bastante:
    y_sub = 1.0d0 / (     (   x(i) - q(2)  )**2 + q(3) ) 
    x_sub = x(i) - q(2)

    f(1) =          y_sub               *(A*y_sub -y(i))        
    f(2) = 2.0d0*A*(y_sub**2)*x_sub     *(A*y_sub -y(i))   
    f(3) =       A*(y_sub**2)           *(A*y_sub -y(i))

    df1(1)    =   f(1)*f(1)                     ;     df1(2)    =  f(1)*f(2)* f(2)/q(1)                               
    df1(3)    =  f(1)*f(3)* f(3)/q(1)           ;     df2(1)    = f(2)*f(1)*  2*(y_sub**2)*x_sub  ;      
    df2(2)    =  f(2)*f(2)* 2*A*(y_sub**2)*(4*x_sub*y_sub +1.0d0)  ;  df2(3)    =  f(2)*f(3)* (-4)*x_sub*A*y_sub**3;   
    df3(1)    = f(3)*f(1)*  (y_sub**2)          ;     df3(2)    =  f(3)*f(2)* 4*A*(y_sub**3)*x_sub                    ;  
    df3(3)    =  f(3)*f(3)* (-2)*A*y_sub**3;   





    do coluna   = 1,3

        Jac(1,coluna) = Jac(1,coluna) +  df1(coluna)
        Jac(2,coluna) = Jac(2,coluna) +  df2(coluna)
        Jac(3,coluna) = Jac(3,coluna) +  df3(coluna)

    end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Achar inversa !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call lu(Jac,L,U)



yes(1,1) = 1.0d0  ;                     yes(1,2)  =   0d0 ;   yes(1,3)  =   0d0;

yes(2,1) = -L(2,1);                     yes(2,2)  =   1d0 ;   yes(2,3)  =   0d0;
yes(3,1) = -L(3,2)*yes(2,1) - L(3,1);  yes(3,2) =   -L(3,2);    yes(3,3)  =   1d0;

JacI(3,1) = yes(3,1) /U(3,3)                                   ;   JacI(3,2) = yes(3,2)/U(3,3); JacI(3,3) = yes(3,3)/U(3,3);

JacI(2,1) = (   yes(2,1) - U(2,3)*JacI(3,1)    ) / U(2,2)  ;   
JacI(2,2) = (   yes(2,2) - U(2,3)*JacI(3,2)    ) / U(2,2) ; 
JacI(2,3) = (   yes(2,3) - U(2,3)*JacI(3,3)    ) / U(2,2);

JacI(1,1) = (   yes(1,1) - U(1,3)*JacI(3,1) - U(1,2)*JacI(2,1)      )   /U(1,1);  
JacI(1,2) = (   yes(1,2) - U(1,3)*JacI(3,2) - U(1,2)*JacI(2,2)      )   /U(1,1); 
JacI(1,3) = (   yes(1,3) - U(1,3)*JacI(3,3) - U(1,2)*JacI(2,3)      )   /U(1,1); 





print*,matmul(Jac,JacI)             !Está correto até 15 casas decimais





end do
end program


