! This program calculates functions derivates through several methods
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 23.07.2018
!
! INPUT: Initial and final x
! 
! OUTPUT: First derivatives using 3 methods and second derivatives using 2 methods





module variables

!Parameters
real*8      ::  h1 = 0.10d0
real*8      ::  h2 = 0.25d0
real*8      ::  e =  2.718281828459d0  
real*8      ::  xi = 0.0d0 
real*8      ::  xf = 5.0d0
real*8      ::  dx = 0.01d0



!Global Variables
integer*8   :: i, npassos
real*8		:: df2,df3,df5,erro2,erro3,erro5,df_exato,x,erro_temp,d2f3,d2f5,d2f_exato
real*8      :: y,y_1,y_2,y1_,y2_                                !y_1 significa y(i+1), ao passo que y1_ significa y(i-1). O mesmo vale para y_2 e y2_	



end module variables


subroutine funcao(x,y,y_1,y_2,y1_,y2_,h)

real*8      ::  x,y,y_1,y_2,y1_,y2_,h
real*8      ::  e =  2.718281828459d0  

y = x*e**x

y_1 = (x+h)*e**(x+h)
y_2 = (x+2.0d0*h)*e**(x+2.0d0*h)

y1_ = (x-h)*e**(x-h)
y2_ = (x-2.0d0*h)*e**(x-2.0d0*h)

return
end subroutine funcao





program derivada
use variables

npassos = (xf - xi) / dx

open(1,file='h1.dat')
open(2,file='h2.dat')
open(3,file='erro.dat')

open(4,file='d2h1.dat')
open(5,file='d2h2.dat')
open(6,file='d2erro.dat')




x = xi
do i = 1, npassos
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIMEIRA DERIVADA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!! h = 0.10 !!
    call funcao(x,y,y_1,y_2,y1_,y2_,h1)                                                 !Recebe x e h e retorna y,y(x+h),y(x+2h),y(x-h) e y(x-2h)
    df2 = ( y_1 - y  )                           / h1
    df3 = (y_1 - y1_)                            / (2.0d0*h1)
    df5 = (y2_ - 8.0d0*y1_ + 8.0d0*y_1 - y_2)    / (12.0d0*h1)

    df_exato = e**x + x*e**x

    erro2 = df2 - df_exato
    erro3 = df3 - df_exato
    erro5 = df5 - df_exato

    if(erro2 < 0) then
        erro2 = erro2*(-1.0d0)
    end if
    if(erro3 < 0) then
        erro3 = erro3*(-1.0d0)
    end if
    if(erro5 < 0) then
        erro5 = erro5*(-1.0d0)
    end if

    erro_temp = erro5
    write(1,*)x,df2,df3,df5,erro2,erro3,erro5,y



!!! h = 0.25 !!
    call funcao(x,y,y_1,y_2,y1_,y2_,h2)                                                 !Recebe x e h e retorna y,y(x+h),y(x+2h),y(x-h) e y(x-2h)
    df2 = ( y_1 - y  )                           / h2
    df3 = (y_1 - y1_)                            / (2.0d0*h2)
    df5 = (y2_ - 8.0d0*y1_ + 8.0d0*y_1 - y_2)    / (12.0d0*h2)

    df_exato = e**x + x*e**x

    erro2 = df2 - df_exato
    erro3 = df3 - df_exato
    erro5 = df5 - df_exato

    if(erro2 < 0) then
        erro2 = erro2*(-1.0d0)
    end if
    if(erro3 < 0) then
        erro3 = erro3*(-1.0d0)
    end if
    if(erro5 < 0) then
        erro5 = erro5*(-1.0d0)
    end if

    write(2,*)x,df2,df3,df5,erro2,erro3,erro5,y

    write(3,*)x,erro_temp,erro5



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SEGUNDA DERIVADA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!   Teoricamente poderiamos fazer este cálculo junto com o calculo da primeira derivada, o que pouparia o número de calculos feito pelo programa. No entando, como o programa é extremamente leve e 
!ser de fácil compreensão é um objetivo importante, iremos separar o calculo da segunda derivada em uma seção propria.

    d2f_exato = 2.0d0*e**x + x*e**x




    call funcao(x,y,y_1,y_2,y1_,y2_,h1)   
   
    d2f3 = (y_1 -2.0d0*y  + y1_) / (h1*h1)
    d2f5 = (-1.0d0*y2_ + 16.0d0*y1_ - 30.0d0*y + 16.0d0*y_1 - 1.0d0*y_2) / (12.0d0*h1*h1)      

    erro3_2 = d2f3 - d2f_exato
    erro5_2 = d2f5 - d2f_exato

    if(erro3_2 < 0) then
        erro3_2 = erro3_2*(-1.0d0)
    end if
    if(erro5_2 < 0) then
        erro5_2 = erro5_2*(-1.0d0)
    end if

    write(4,*)x,d2f3,d2f5,erro3_2,erro5_2



    erro_temp = erro5_2


    call funcao(x,y,y_1,y_2,y1_,y2_,h2)  

    d2f3 = (y_1 -2.0d0*y  + y1_) / (h2*h2)
    d2f5 = (-1.0d0*y2_ + 16.0d0*y1_ - 30.0d0*y + 16.0d0*y_1 - 1.0d0*y_2) / (12.0d0*h2*h2)      

    erro3_2 = d2f3 - d2f_exato
    erro5_2 = d2f5 - d2f_exato

    if(erro3_2 < 0) then
        erro3_2 = erro3_2*(-1.0d0)
    end if
    if(erro5_2 < 0) then
        erro5_2 = erro5_2*(-1.0d0)
    end if

    write(5,*)x,d2f3,d2f5,erro3_2,erro5_2


    write(6,*)x,erro_temp,erro5_2


x = x + dx
end do










end program


