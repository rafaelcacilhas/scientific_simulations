program histograma
use Marsaglia

!Parameters

integer*8,parameter     ::   Mcs      =  100
integer*8,parameter     ::   ncaixas  =  30


!Global Variables

real*8		:: energia(1:Mcs),maior_energia,menor_energia,dh,H(1:ncaixas)
integer*8   :: i,k,s,media,caixa

open(1,file="energia30.dat")                              ! Arquivo de entrada
open(2,file="hist30.dat")                                 ! Arquivo de saída
maior_energia = -20.0d0
menor_energia = 999.0d0
H = 0d0

do k = 0,Mcs
    read(1,*)i,energia(k)
    if(energia(k) > maior_energia) maior_energia = energia(k)
    if(energia(k) < menor_energia) menor_energia = energia(k)
end do

dh = (maior_energia - menor_energia) / (1.0d0*ncaixas)

do k = 0,Mcs
    caixa = (energia(k) - menor_energia) / dh 
    H(caixa) = H(caixa) + 1
end do

do k = 1,ncaixas
    write(2,*)dh*k,H(k)
end do


end program histograma

!============================================================================
