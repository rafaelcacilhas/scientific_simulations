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
integer*4,parameter   ::  N = 100

!Global Variables


real*8      :: xn,xn1,lambda,r,valor(4),x0
integer*8   :: i,nome,j
character(len=10) :: file_id
character(len=50) :: file_name

end module variables




program expoente
use variables



x0 = 0.1d0
xn = x0

valor(1) = 0.85d0 
valor(2) = 0.87d0
valor(3) = 0.89d0
valor(4) = 0.96d0

do j = 1,4
!Saída de dados
nome = 100*valor(j)
write(file_id, '(i0)') nome
file_name = '0' // trim(adjustl(file_id)) // '.dat'
open(1,file = trim(file_name))
print*,file_name



do i = 1,N
r = valor(j)
        xn1 = 4.0d0*r*xn*(1d0-xn)
        xn  = xn1
write(1,*)i,xn1
end do

close(1)
end do



end program

    
