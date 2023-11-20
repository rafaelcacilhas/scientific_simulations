! This program solves 
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 11.08.2018
!
! INPUT:   
! 
! OUTPUT:   





module variables


!Parameters

real*8,parameter      ::  e = 2.718281828459045235
integer*4             ::  N = 40

!Global Variables


real*8      :: xn,xn1,lambda,r,x0,exato
integer*8   :: i,nome,j
character(len=10) :: file_id
character(len=50) :: file_name

end module variables




program expoente
use variables

x0 = 0.1d0
xn = x0
r = 0.33d0




! Saída de dados:
nome = 10*x0
write(file_id, '(i0)') nome
file_name = '0' // trim(adjustl(file_id)) // '.dat'
open(1,file = trim(file_name))
print*,file_name


do j = 1,5

do i = 1,N

xn1 = 4.0d0*r*xn*(1d0-xn)
xn  = xn1



end do
exato = 1.0d0 - (1.0d0/(4.0d0*r))
write(1,*)1.0d0/r,xn1,exato
r = r + 0.1d0
end do

close(1)

end program

    
