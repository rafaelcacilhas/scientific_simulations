! This program calculates roots using the Newton method
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 23.07.2018
!
! INPUT: The initial X value, x0 
! 
! OUTPUT: Root.

!numero de iterações: 
!x0 = -2 -> 4
!x0 = -3 -> 5
!x0 = -4 -> 7




module variables

!Parameters

real*8      ::  e =  2.718281828459d0  
real*8      ::  xi = -40.0, xf = 40.0
real*8      ::  x0 = -2.0d0
real*8      ::  deltaX = 1.0
real*8      ::  epslon = 10.0E-5
real*8      ::  erroM = 10e-5



!Global Variables
integer*8   :: i, npassos,j,teste1,teste2,nraizes,contador
real*8		:: x,xn,xn1,erro,derivada,funcao,raiz,raizes(1:5),erros(1:5),seeds(1:5)


end module variables



program newton
use variables

npassos = (xf - xi)/ deltaX
erro = 1.0
raizes = 9999                                                               !Defini este valor para que ao analilsar o array raizes o valor padrão 0 não seja confundido com uma raiz.
nraizes = 0


open(2,file="2.dat")
open(3,file="3.dat")
open(4,file="4.dat")


do i=1,npassos,1
x = x0
contador = 0

    do
    !write(2,*)contador,x
    contador = contador + 1

    raiz = x

    funcao =    (    1.0d0   +   (1.0d0+x**2)*dsin(x/5.0d0)    )   /   (1.0d0+x**2)
    derivada =  (    ((1.0d0 + x**2)**2 )*dcos(x/5.0d0) - 10.0d0*x  )   /           (   5.0d0*(1.0d0 + x**2)**2    )

    x = x - funcao/derivada

    erro = (raiz-x)/x
    if(erro < 0) then
        erro = -1.0d0*erro
    end if







    if ( erro < erroM) then     

        do j = 1,nraizes,1                                  !Checa se a raiz encontrada ja foi encontrada anteriormente
            teste1 = raiz
            teste2 = raizes(j)
            !print*,teste1,teste2
            if(teste1 == teste2) then
                goto 32
             endif
        end do
        

        if(x < xi .or. x > xf) then                         !Checa se a raiz encontrada está no intervalo de interesse
            goto 32
        end if
        nraizes = nraizes + 1  
        raizes(nraizes) = x
        erros(nraizes) = erro
        seeds(nraizes) = x0
        print*,"Raiz: ",raizes(nraizes),"Erro: ", erros(nraizes), "X inicial: ", seeds(nraizes), "numero de iterações: ",contador


        32 continue
        exit
    end if


    
    end do


x0 = x0 + deltaX
end do

end program


