! This program calculates the difference in a subtraction when the terms are truncated and the "exact" solution
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 2017.18.15
!
! INPUT: two numbers
! 
! OUTPUT: the absolute and relative differences between the inputs.

module variables

!Parameters
real*8		::	x1 = 22.0d0/7.0d0
real*8		::	x2 = 60.0d0/19.0d0

!Global Variables
real*8		::	delta,delta5,erro_ab,erro_rel,x1T,x2T               !x1T e x2T representam x1 Truncado e x2 Truncado

end module variables


program error
use variables



x1 = 22.0d0/7.0d0
x2 = 60.0d0/19.0d0



temp = 10e3 * x1
x1T = (temp*1d0 / 10e3)

temp = 10e3 * x2
x2T = (temp*1d0 / 10e3)






delta = x1 - x2
delta5 = (x1T - x2T) 

erro_ab = delta - delta5
    if(erro_ab < 0.0) then
        erro_ab = erro_ab*(-1.0d0)
    end if

erro_rel = (delta - delta5)/delta
    if(erro_rel < 0.0) then
        erro_rel = erro_rel*(-1.0d0)
    end if


write(*,*)"Delta exato: ",delta
write(*,*)"Delta truncado: ",delta5
write(*,*)"Erro absolulto: ",erro_ab
write(*,*) "Erro relativo: ",erro_rel

end program


