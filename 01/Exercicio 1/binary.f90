! This program converts a binary number in a decimal one  and vice-versa
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 2017.18.15
!
! INPUT: binary number in a number.dat file (1 character per line)
! 
! OUTPUT: binary conversion on screen and decimal number conversion in a bin.dat file (1 character per line and decimals after the comma)





module variables

!Parameters
real*4      ::  numero = 1227.86d0
integer*8   ::  precisao = 10


!Global Variables
integer*8	::	bin(32),i,sign, expo(0:7),mantiss(1:23),c = 0,c1 = 127,numero_int,numero_intT,resto,binaI(1:20) = 0,binaD(1:20) = 0,j
integer*8   :: imax,temp(1:20) = 0
real*8		::	real_number,f = 0.0,dec


real*8      ::  dec_temp,expoente
integer*4   ::  decimal(20)
integer*16  ::  number_bin,sub

end module variables




program binary
use variables


open(1,file='number.dat')
read(1,*)bin          
close(1)


!do i = 1,16                                                                ! Este trecho de código foi uma tentativa de ler o numero em binario sem precisar de um arquivo externo. Não funcionou pois o fortran não aceita
                                                                            !  numeros inteiros da ordem de 10**20
!    expoente = 31-i
!    number_bin_temp = number_bin/10**expoente    

!    sub = 0
!    do j = 1,i
!    sub = sub - bin(j )*10**(i-j+1)
!    end do

!    number_bin_temp = number_bin_temp + sub
!    bin(i+1) = number_bin_temp

!    print*,number_bin_temp,sub

!end do



sign = (-1)**bin(1)


do i=0,7
	expo(i) = bin(9-i)
end do


do i=1,23
	mantiss(i) = bin(9+i)
end do

		


do i = 0,7
	c = c + expo(i)*(2**i)
end do


do i = 1,23
	f = f + mantiss(i)*(2.0**(-1.0*i) )
end do	






real_number = sign*(2.0**(c-c1))*(1.0+f)
write(*,*) "O numero binário dado escrito na forma decimal e: " ,real_number





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  PARTE 2:    Decimal para Binário    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


decimal = 0.d0

numero_int = numero                 !numero_int é definido como inteiro, portanto ele ignora a parte decimal
dec = numero - 1d0*(numero_int)

numero_intT = numero_int            !Variavel temporaria para que numero_int seja constante ao longo do programa

i = 1
do while (numero_intT >= 1)         !Calcula os numeros binarios
    resto = mod(numero_intT,2)
    binaI(i) = resto
    numero_intT = numero_intT/2
    i = i+1 
end do
imax = i-1

do i = 1,imax                       !O numero é criado invertido. Este do serve para inverter o array bina(i)
    temp(i) = binaI(imax+1-i)
end do
binaI = temp



dec_temp = dec                      !Variavel temporaria para que dec seja constante ao longo do programa

do i = 1,precisao                   !A variavel precisao define quantas casas "decimais" serao utilizadas. Ela pode ser variada  até 20 sem nenhum problema; acima disso é necessario aumentar o tamanho do array decimal(20)
    dec_temp = dec_temp * 2.0d0

    if(dec_temp >= 1.0d0) then
    decimal(i) = 1
    dec_temp = dec_temp - 1.0d0    

    else if (dec_temp < 1.0d0) then
    decimal(i) = 0
    end if


end do



open(2,file = "bin.dat")

do i = 1,imax
    write(2,*) binaI(i)
end do

write(2,*)","

do i = 1,precisao
    write(2,*) decimal(i)
end do





end program


