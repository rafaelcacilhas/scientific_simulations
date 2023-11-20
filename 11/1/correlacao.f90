program correlacao
use Marsaglia

!Parameters
integer*8,parameter     ::   itotal   =  2**20                !numero de linhas no arquivo de entrada
integer*8,parameter     ::   Mcs      =  2**20 -1
integer*8,parameter     ::   nmedia   =  itotal/Mcs
integer*8,parameter     ::   imax     =  Mcs/5

!Global Variables

real*8		:: energia(1:Mcs),energiaMED(1:Mcs),energia_media,sigma2,soma,C(0:Mcs),tau,somatau,C_medio(0:Mcs),tau_medio
integer*8   :: i,k,s,media

open(1,file="random_seq.dat")                              ! Arquivo de entrada
open(2,file="ci.dat")                                 ! Arquivo de saída
 C_medio   = 0d0
 tau_medio = 0d0



            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Média e Desvio padrão !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            soma = 0d0
            do k = 0,Mcs
                read(1,*)energia(k)
                soma = soma + energia(k)
            end do

            energia_media = soma / (1.0d0*Mcs)

            soma = 0d0
            do k = 1,Mcs
                soma = soma + (   energia(k) - energia_media  )**2
            end do

            sigma2      =  soma / (1.0d0*Mcs)
            print*,energia_media,sigma2

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do i = 0,imax
            soma = 0d0
                do s = 0,Mcs-i
                    soma = soma + ( energia(s) - energia_media  )*(   energia(s+i) - energia_media    )
                end do

                C(i) =  soma / (   sigma2*1.0d0*(Mcs - i)  )
                somatau = somatau + C(i)
                write(2,*)i,C(i)/C(0)
            end do
            
            tau = somatau - 0.5d0*C(0)          !! O sinal de menos aqui é porque somamos 1C0 no somatório acima, então estamos subtraindo 0.5C0
            print*,2.0d0*tau


end program correlacao

!============================================================================
