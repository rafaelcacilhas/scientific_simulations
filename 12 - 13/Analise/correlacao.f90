program correlacao

!Parameters
integer*8,parameter     ::   terma    =  10e5
integer*8,parameter     ::   itotal   =  10e6                 !numero de linhas no arquivo de entrada
integer*8,parameter     ::   ndata    =  itotal - terma
integer*8,parameter     ::   nmedia   =  itotal
integer*8,parameter     ::   imax     =  ndata/10000
integer*8,parameter     ::   N        =  32*32

!Global Variables

real*8		:: energia(1:ndata),energia_media,sigma2,soma,C(0:ndata),tau,somatau,C_medio(0:ndata),tau_medio,mag_medio,C_mag(0:ndata)
real*8      :: sigma2_mag,magnetizacao(1:ndata),soma2,lixo1,lixo2,lixo3,histograma_energia(4*N),histograma_magnetizacao(-N:N)
real*8      :: somatau2
integer*8   :: i,k,s,media,l,l2
open(1,file="exercicio2.dat")                                        ! Arquivo de entrada
open(2,file="ci2.dat")                                      ! Arquivo de saída
open(5,file="energia2.dat")
open(7,file="mag2.dat")
 C_medio   = 0d0
 tau_medio = 0d0
histograma_magnetizacao = 0
histograma_energia      = 0

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Média e Desvio padrão !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            soma  = 0d0
            soma2 = 0d0
        
            do k = 1,terma
                read(1,*)lixo1,lixo2,lixo3
            end do

            do k = 1,ndata
                read(1,*)i,energia(k),magnetizacao(k)
                soma  = soma  + energia(k)
                soma2 = soma2 + magnetizacao(k)  


                l  = (   (2*N - energia(k)    ) / 4     )    +1
                l2 = (magnetizacao(k))
                histograma_magnetizacao( l2) = histograma_magnetizacao(l2) + 1         
                histograma_energia(l) = histograma_energia(l) + 1 
            end do

histograma_magnetizacao = histograma_magnetizacao / ndata
histograma_energia      = histograma_energia / ndata

            energia_media = soma  / (1.0d0*ndata)
            mag_medio     = soma2 / (1.0d0*ndata) 

            soma  = 0d0
            soma2 = 0d0
            do k = 1,ndata
                soma  = soma  + (   energia(k)      - energia_media  )**2
                soma2 = soma2 + (   magnetizacao(k) - mag_medio      )**2
            end do

            sigma2      =  soma  / (1.0d0*ndata)
            sigma2_mag  =  soma2 / (1.0d0*ndata)

            print*,"Energia: ",energia_media,dsqrt(sigma2)
            print*,"Mag: ",mag_medio,dsqrt(sigma2_mag)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            somatau  = 0d0
            somatau2 = 0d0

            do i = 0,imax
            soma  = 0d0
            soma2 = 0d0
                do s = 0,ndata-i
                    soma  = soma  + ( energia(s)      - energia_media  )*(   energia(s+i)      - energia_media    )
                    soma2 = soma2 + ( magnetizacao(s) - mag_medio      )*(   magnetizacao(s+i) - mag_medio        )
                end do

                C(i)     =  soma  / (   sigma2    *1.0d0*(ndata - i)  )
                C_mag(i) =  soma2 / (   sigma2_mag*1.0d0*(ndata - i)  )
                somatau  = somatau  + C(i)
                somatau2 = somatau2 + C_mag(i)

                write(2,*)i,C(i)/C(0)
            end do
            
            tau = somatau2 - 0.5d0*C_mag(0)          !! O sinal de menos aqui é porque somamos 1C0 no somatório acima, então estamos subtraindo 0.5C0
            print*,2.0d0*tau




do k = 1, N
    write(5,*)( -4.0d0*(k-1) + 2*N  ),histograma_energia(k)
end do

do k = -N,N
write(7,*)k,histograma_magnetizacao(k)
end do



end program correlacao

!============================================================================
