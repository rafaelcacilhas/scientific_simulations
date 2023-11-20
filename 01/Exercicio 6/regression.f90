! This program calculates a linear regression from a data set
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 2017.18.16
!
! INPUT: Data set with x and y
! 
! OUTPUT: angular and linear coefficients for a linear regression


!xmgrace results: y = 1.6382 * x + 0.029209  

module variables

!Parameters
real*8      ::  electron_charge = 1.60217662E-19

integer*4   ::  N = 15


!Global Variables
integer*4   ::  x(15),A,C,i
real*8      ::  y(15),B,D,c_angular,c_linear,c_angular_2,erro_ab,erro_rel


end module variables


program regression
use variables
A = 0                                   !Essas variaveis serão definidas por uma soma, portanto é necessario zera-las
B = 0d0
C = 0
D = 0d0

open(1,file="millikan.dat")

do i=1,15
    read(1,*)x(i),y(i)
end do


do i=1,15
    A = A + x(i)
    B = B + y(i)
    C = C + (   x(i)*x(i)   )
    D = D + (   x(i)*y(i)   )   

end do


c_angular = (N*D - A*B) /   (N*C - A*A)
c_linear =  (B*C - A*D) /   (N*C - A*A)

print*, c_angular, " x + ",c_linear
print*,""


c_angular_2 = (N*D - A*A*D/C) / (N*C - A*A)
print*,"Impondo b = 0, temos que a = ",c_angular_2



qm = 10E18*electron_charge
erro_ab = qm - c_angular_2
erro_rel = 100*(-qm + c_angular_2) / qm
print*,""
print*, "Valor de a esperado: ",qm
print*,"Erro relativo: ",erro_rel, "%"




end program
