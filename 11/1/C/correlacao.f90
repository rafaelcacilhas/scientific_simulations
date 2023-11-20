program correlacao


!Parameters
integer*8,parameter     ::   itotal   =  2**20                !numero de linhas no arquivo de entrada
integer*8,parameter     ::   Mcs      =  100
integer*8,parameter     ::   nmedia   =  itotal/Mcs
integer*8,parameter     ::   imax     =  Mcs/5

!Global Variables

real*8		:: energia(1:Mcs),energiaMED(1:Mcs),energia_media,sigma2,soma,C(0:Mcs),tau,energia2
integer*8   :: i,k,s,media

open(1,file="d=15.dat")                              ! Arquivo de entrada
open(2,file="ci15.dat")                                 ! Arquivo de saída
 C_medio   = 0d0
 tau_medio = 0d0
 energia2  = 0d0
 soma      = 0d0

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Média e Desvio padrão !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do k = 0,Mcs
                read(1,*)i,energia(k)
                soma = soma + energia(k)
                energia2 = energia2 + energia(k)**2
                !print*,soma
            end do

            energia_media = soma     / (1.0d0*Mcs)
            energia2      = energia2 / (1.0d0*Mcs) 

            soma = 0d0
            do k = 1,Mcs
                soma = soma + (   energia(k) - energia_media  )**2
            end do

            sigma2      =  soma / (1.0d0*Mcs)
 !           print*,energia_media/200.0d0

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

open(3,file="cv.dat")
write(3,*)-15,83.470415445129106
write(3,*)0,92.636588827426408  
write(3,*)15,76.569347040140656
write(3,*)30,111.71098394023124 

open(4,file="energia.dat")
write(4,*)-15,0.12375567595600806 
write(4,*)0,0.13234361967821229 
write(4,*)15,0.12307613507149399     
write(4,*)30,0.13638457599479981  



end program correlacao

!============================================================================
