! This program calculates sum with two different precisions
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 2017.18.15
!
! INPUT: Sum initial and final term
! 
! OUTPUT: Sum 


module variables

!Parameters
integer*8	::	ki = 1
integer*8	::	kf = 100


!Global Variables

integer*8	::	i
real*8		::	sum8 = 0.0;
real*4      ::  sum4 = 0.0;

end module variables




program sum
use variables


write(*,*) "Conta feita começando pelos menores termos: "

do i=ki,kf

	sum4 = sum4 + (i**5)*(-1)**i
	sum8 = sum8 + (i**5)*(-1)**i

end do

write(*,*)"Real*4: ", sum4
write(*,*)"Real*8: ", sum8


write(*,*) " "
write(*,*) "Conta feita começando pelos maiores termos: "



do i=kf,ki,-1

	sum4 = sum4 + (i**5)*(-1)**i
	sum8 = sum8 + (i**5)*(-1)**i

end do

write(*,*)"Real*4: ", sum4
write(*,*)"Real*8: ", sum8




end program


