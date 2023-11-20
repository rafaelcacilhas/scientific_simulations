! This program solves the Wheatstone bridge
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 28.08.2018
!
! INPUT:    Six resistances values and the external potential
! 
! OUTPUT:   Current in all circuits resistor 





module variables

!Parameters
real*8      ::  V0 =    1.5d0
real*8      ::  r1 =  100.0d0
real*8      ::  r2 =  100.0d0
real*8      ::  r3 =  150.0d0
real*8      ::  rx =  120.0d0
real*8      ::  ra = 1000.0d0
real*8      ::  rs =   10.0d0


!Global Variables

real*8      :: R(1:3,1:3),U(1:3,1:3),L(1:3,1:3), i(3),V(1:3),y(1:3)
real*8		:: i4,i5,ia
real*8      :: u1(1:3,1:3),l1(1:3,1:3),soma = 0d0
integer*4   :: linha,coluna

end module variables






program matriz
use variables



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definir as matrizes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

R(1,1)  =  rs   ;   R(1,2)  = r1            ;   R(1,3)  = r2            ;
R(2,1)  = -rx   ;   R(2,2)  = (r1+rx+ra)    ;   R(2,3)  = -ra           ;
R(3,1)  = -r3   ;   R(3,2)  =  -ra          ;   R(3,3)  = (r2+r3+ra)    ;




U(1,1) =   R(1,1)   ;   U(1,2) =  R(1,2)                         ;   U(1,3) =  R(1,3)      ;
U(2,1) =   0d0      ;   U(2,2) = R(2,2) - R(2,1)*R(1,2)/R(1,1)   ;   U(2,3) =  R(2,3) - R(2,1)*R(1,3)/R(1,1)  ;
U(3,1) =   0d0      ;   U(3,2) =   0d0                           ;   
U(3,3)=R(3,3) -1.0d0*(R(3,1)*R(1,3)/R(1,1) +  ( R(3,2) - R(3,1)*R(1,2)/R(1,1)   )*(     R(2,3) - R(2,1)*R(1,3)/R(1,1)    )/U(2,2))



	
L(1,1) =  1.0d0            ;   L(1,2) =   0d0                                         ;   L(1,3) =  0d0  ;
L(2,1) =  R(2,1)/R(1,1)    ;   L(2,2) =   1.0d0                                       ;   L(2,3) =  0d0   ;
L(3,1) =  R(3,1)/R(1,1)    ;   L(3,2) = ( R(3,2) - R(3,1)*R(1,2)/R(1,1)  )   / U(2,2) ;   L(3,3) = 1.0d0  ;


V(1) = V0;
V(2) = 0d0;
V(3) = 0d0;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Solução do sistema  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y = 0d0
do linha=1,3
    soma = 0d0
    do coluna=1,linha-1
        soma = soma - L(linha,coluna)*y(coluna)
    end do  
    y(linha) = V(linha) + soma
end do
print*,y



i(3) = (y(3)                                )/U(3,3)
i(2) = (y(2) - U(2,3)*i(3)                  )/U(2,2)
i(1) = (y(1) - U(1,3)*i(3) - U(1,2)*i(2)    )/U(1,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Outras correntes  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

i4 = i(1) - i(2)
i5 = i(1) - i(3)
ia = i(3) - i(2) 
print*, "Correntes 1, 2 e 3: ",i(1),i(2),i(3)
print*,"Corrente no amperímetro: ",ia



end program


