! This program integrates functions
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update:  23.07.2018
!
! INPUT: initial and final x
! 
! OUTPUT: integral results and error





module variables

!Parameters

real*8      ::  xi        = 0.0d0 
real*8      ::  xf        = 3.14159265359/2.0d0




!Global Variables
integer*8   :: i,npassos,k
real*8		:: f,f1,x,y,y1,integral,h,valor_exato,erro_100,erro_1000,y2,y_1,y1_,temp



end module variables


subroutine funcao(x,y,y1,y2,h)

real*8      ::  x,y,y1,y2,h


y = dsin(x)
y1 = dsin(x+h)
y2 = dsin(x+2*h)



return
end subroutine funcao




program derivada
use variables

valor_exato = 1.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! n = 100 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111
npassos = 100
h = (xf-xi)/npassos

integral = 0d0
do i =0,(npassos-1)

    call funcao(x,y,y1,y2,h)
    integral = integral + (y1+y)*h/2.0d0


x = x + h
end do
erro_100 = integral - valor_exato
temp = integral



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! n = 1000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integral = 0d0
x = 0d0
npassos = 1000
h = (xf-xi)/npassos


do i =0,(npassos-1)

    call funcao(x,y,y1,y2,h)
    integral = integral + (y1+y)*h/2.0d0


x = x + h

end do
erro_1000 = integral - valor_exato

print*,"Método do trapézio: "
print*,"Valor da integral: "
print*,"N = 100: ",temp,"N = 1000: ",integral
print*,"Erro para 100 passos: ",erro_100,"Erro para 1000 passos: ",erro_1000






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simpsons !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integral = 0d0
x = xi
npassos = 2*100                                     !Como o algoritmo faz npassos/2 iterações, definiremos npassos = 2*100 
h = (xf-xi)/npassos

do k =1,npassos/2

    y  = dsin(x                    )
    y1 = dsin(x    +         h     )
    y2 = dsin(x    +   2.0d0*h     )


    integral = integral + (y + 4*y1 + y2)*h/3.0d0


x = x + 2*h
end do
erro_100 = integral - valor_exato
temp = integral


integral = 0d0
x = xi
npassos = 2*1000                                        !Como o algoritmo faz npassos/2 iterações, definiremos npassos = 2*1000 
h = (xf-xi)/npassos
do k =1,npassos/2

    y  = dsin(x                    )
    y1 = dsin(x    +         h     )
    y2 = dsin(x    +   2.0d0*h     )


    integral = integral + (y + 4*y1 + y2)*h/3.0d0


x = x + 2*h
end do
erro_1000 = integral - valor_exato


print*,"Metodo de Simpsons: "
print*,"Valor da integral: "
print*,"N = 100: ",temp,"N = 1000: ",integral
print*,"Erro para 100 passos: ",erro_100,"Erro para 1000 passos: ",erro_1000











end program


