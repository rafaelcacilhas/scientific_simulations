! This program uses the partition function to calculate energy values
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 2017.18.15
!
! INPUT: Data set with k and Omega(k)
! 
! OUTPUT: 



subroutine soma(A,B,Z)
real*8 A,B,Z

Z = max(A,B) + log(     1   +   exp(    min(A,B)    -   max(A,B)    )   )

end subroutine soma




program partition



!Constants
real*8      ::  e       = 2.718281828459d0  
real*8      ::  kb      = 1.38064852E-23   
real*8      ::  T       = 300d0                     !Temperature in Kelvin
integer*4   ::  N       = 32





!Global Variables

integer*4   ::  i,k,Energy
real*8      ::  A,B,C,O(1024),Z,Z_temp
real*8      ::  beta 
real*16     ::  peso

beta =  1d0
Z = 0d0


open(1,file="LNgE.txt")

do i = 0,10
    read(1,*)k,O(k)

    Energy = 4*k -2*N

    Z = Z +  e**( O(k) - beta*Energy)


    print*,k,Z

end do

call soma(Z,( O(k) - beta*Energy),Z)                    !Conceitualmente errado pois o log da soma não é a soma dos logs







end program


