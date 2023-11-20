! This program 
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 29.10.2018
!
! INPUT:    
! 
! OUTPUT:   





module variables

!Parameters


!Global Variables

real*8                  ::  x(1:2**20),Ci,C0,x_medio,sigma,soma,tau
integer*4               ::  i,j,linha,coluna,s,k,nk



end module variables


program random
use variables


open(1,file="random_seq.dat")
open(2,file="saida.dat")
open(3,file="tau.dat")

soma = 0d0

k = 20
Nk = 2**k

do i=1,Nk
    read(1,*)   x(i)
    write(2,*)i,x(i)
    soma = soma + x(i)
end do

x_medio = soma / Nk


soma = 0d0
do i=1,Nk
    soma = soma + ( x(i) - x_medio)**2
end do
sigma = soma / Nk

print*,x_medio,"+",sigma




do k=1,2**20
    soma  = 0d0
    sigma = 0d0
	tau   = 0d0
	Nk    = 2**k

    do i=1,Nk
        soma = soma + x(i)
    end do
    x_medio = soma / Nk

	do i=1,Nk
		sigma = sigma + (   x(i)    -   x_medio )**2
	end do
	sigma = sigma / Nk



	do i=0,30
		Ci=0
		do j=1,Nk-i
			Ci = Ci + (x(j) - x_medio)*(x(j+i) - x_medio)				
		end do
        if(i==0) C0 = Ci

		tau = tau + Ci/(Nk-i)
	end do

	tau = 0.5*C0 + tau/sigma
	write(3,*) k, 2.0d0*tau
print*,k
enddo	





close(1)


end program


