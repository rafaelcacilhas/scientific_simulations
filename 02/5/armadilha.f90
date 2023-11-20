! This program
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 23.07.2018
!
! INPUT: 
!   
! OUTPUT: Number of molecules





module variables

!Constants

real*8      ::  pi      = 4.0d0*datan(1.0d0) 
real*8      ::  e       = 2.718281828459
real*8      ::  kb      = 1.380648520E-23
real*8      ::  kbEV    = 11604.517650070278
real*8      ::  Na      = 6.02E23 
real*8      ::  eV      = 1.60217610E-22            !10^-19 do eletron volt e 10^-3 do mili


!Parameters

real*8      ::  T1      = 10.0d0                   !temperatura em Kelvin 
real*8      ::  Ne      = 0.3d0
real*8      ::  dx      = 0.1d0
integer*4   ::  n       = 2


!Global Variables
integer*8   :: i,npassos,k
real*8		:: C1,energia,g,v,m,de,C0,integral=0d0,integralI=0d0,a,u,x,integralP=0d0,y,valor_exato


end module variables






program armadilha
use variables

!!!!!!!!!!!!!! Mediremos energia em mili eletronvolts        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


energia = 0d0
de = 0.01d0                                                                                   
npontos = 800



open(1,file="graph.dat")

C1 = 3.0E13 / (dSqrt(pi) *   (kb*T1)**(3.0d0/2.0d0)  )
C0=(Sqrt(eV))



do i=1,npontos
    g = C0*C1*dSqrt(energia)*e**( (-1d0 * eV / (kb*T1)  ) *    energia     )                                  !mili eletron volts 
    write(1,*)energia,g



    energia = energia + de
end do










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Integrar de 0 ate a
a = 10.0d0
h = 10e-4
energia = 0d0
integralP = 0d0
npassos = ( (a - 0.0d0) / h)

valor_exato = 15000000000000.d0

do k =0,npassos/2


    g = C0*C1*dSqrt(energia)*e**( (-1d0 * eV / (kb*T1)  ) *   energia     )
    g1 = C0*C1*dSqrt(energia+h)*e**( (-1d0 * eV / (kb*T1)  ) *   (energia+h)     )
    g2 = C0*C1*dSqrt(energia+2.0d0*h)*e**( (-1d0 * eV / (kb*T1)  ) *    (energia+2.0d0*h)     )

    integralP = integralP + (g + 4*g1 + g2)*h*ev/(3.0d0)                                            !eV precisa entrar aqui devido a mudança de variaveis

energia = energia + 2*h    
end do
print*,energia,g

print*,integralP, valor_exato,-integralP + valor_exato


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Integrar de 1/a ate 0
a = 1d0/a
h = a*10e-4
energia = 0d0
integralI = 0d0
npassos = ( (a - 0.0d0) / h)
y = 1d0/(energia*ev*ev)             ! y deve ser dividido por ev^2 para compensar os evs do lado de fora de y na equação. Divide uma vez pra eliminar o de fora e divide denovo para re-escalar o y.
y = h

do k =0,npassos/2


    g = C0*C1*dSqrt(y)*e**( (-1d0 * eV / (kb*T1)  ) *   y     )
    g1 = C0*C1*dSqrt(energia+h)*e**( (-1d0 * eV / (kb*T1)  ) *   (energia+h)     )
    g2 = C0*C1*dSqrt(energia+2.0d0*h)*e**( (-1d0 * eV / (kb*T1)  ) *    (energia+2.0d0*h)     )

    integralI = integralI + (g + 4*g1 + g2)*h*ev/(3.0d0)                                            !eV precisa entrar aqui devido a mudança de variaveis

energia = energia + 2*h    
end do


print*,integralI+IntegralP, -(integralP+IntegralI) + valor_exato





end program

