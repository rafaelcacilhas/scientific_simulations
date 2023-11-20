! This program makes a P vs V graph for two different sets of P and V values
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 2017.18.15
!
! INPUT: Volume range
! 
! OUTPUT: Volume versus Pressure data file


module variables

!Parameters
integer*4   ::  nmol = 2
integer*4   ::  pontos = 5000
integer*16	::	Na = 6022E20
real*8      ::  kb = 1.380648520E-23
real*8      ::  a = 1.3373d0
real*8      ::  b = 0.0320d0

real*8      ::  Vi = 0.05d0
real*8      ::  Vf = 0.2d0

real*8      ::  Vli = 1d0
real*8      ::  Vlf = 5d0


!Definitions
real*8      ::  a1,b1,R,Tc,N


!Global Variables
integer*4   ::  i,imax,j
real*8		::	sum8 = 0.0,T,p,V,dV,p1,V1,Vc,pc,T1
real*8      ::  Pl,Vl,Tl,dVl




end module variables



program pressure
use variables




!Definitions
a1 = a/Na
b1 = b/Na
R = kb*Na
Tc = 8*a1/(27.0*kb*b1)
N = nmol*Na



open (unit = 1, file = "Tc2.dat")
open (unit = 2, file = "Tc.dat")
open (unit = 3, file = "3Tc2.dat")

open (unit = 4, file = "temp2Tc2.dat")
open (unit = 5, file = "temp2Tc.dat")
open (unit = 6, file = "temp23Tc2.dat")




dV = (Vf-Vi)/(pontos)
V = Vi
dVl = (Vlf-Vli)/(pontos)
Vl = Vli


do i=1,pontos


        T = Tc/2
        p = (   N*kb*T/(V-N*b1)   ) - (   N*N*a1/(V*V)  ) 
        write(1,*) V,p
    
        T = Tc    
        p = (   N*kb*T/(V-N*b1)   ) - (   N*N*a1/(V*V)  ) 
        write(2,*) V,p
        
        T = 3*Tc/2
        p = (   N*kb*T/(V-N*b1)   ) - (   N*N*a1/(V*V)  ) 
        write(3,*) V,p


        V = V + dv
end do

close(1)
close(2)
close(3)


do i=1,pontos


        Tl = 0.5d0
        Pl = (   Tl/(   3d0*Vl-1d0  )   ) - (   3d0/(8d0*Vl*Vl)  ) 
        write(4,*) Vl,pl
    
        Tl = 1d0
        Pl = (   Tl/(   3d0*Vl-1d0  )   ) - (   3d0/(8d0*Vl*Vl)  ) 
        write(5,*) Vl,pl
        
        Tl = 1.5d0
        Pl = (   Tl/(   3d0*Vl-1d0  )   ) - (   3d0/(8d0*Vl*Vl)  ) 
        write(6,*) Vl,pl


        Vl = Vl + dVl
end do

close(4)
close(5)
close(6)




end program


