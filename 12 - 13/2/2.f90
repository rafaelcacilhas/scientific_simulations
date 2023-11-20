! Modulo e subrotinas adapatados por:
! Leandro G. Rizzi (lerizzi at ufv.br)
! (Maio de 2017)

	module Marsaglia
                logical, parameter :: seed_list = .FALSE.	! se "seed_list = .TRUE." le arquivos "input.dat" e "seeds.txt"
                logical, parameter :: get_seed  = .FALSE.	! se "get_seed  = .TRUE." le arquivo "seedin.rng"

! Random Number Generator [RNG]:
		type :: RNG_seed_state
			integer :: ij=1802, kl=9373     ! default values
			integer :: i, j
			double precision :: u(1:97), c, cd, cm
		end type

		type (RNG_seed_state) :: raset	! variavel que define o estado do RNG

		integer :: Nseed

	contains

!----------------------------------------------------------------------------

		subroutine rmaset()
! Initializing routine for ranmar or ravmar, must be called
! before generating pseudorandom numbers with ranmar (ravmar).
! (James version of Marsaglia and Zaman, FSU-SCRI-87-50)
! ranges: 0<=ij<=31328 AND 0<=kl<=30081.
		implicit double precision (a-h,o-z)
		integer :: np

		Nseed = 0
		if (seed_list) then
		        open(27,file="input.dat")
	        	read(27,*) Nseed
	        	close(27)
	        	if (Nseed.GT.3376) then
	        	        write(*,*) "ERROR: Nseed .GT. 3376"
	        	        stop
	        	end if
	        	open(32,file="seeds.txt")
	        	do np=1,Nseed
	        	        read(32,*) raset%ij, raset%kl
	        	end do
	        	close(32)
		end if

		if (get_seed) then
			open(78,file="seedin.rng")
			read(78,*) raset
			close(78)
			Nseed = -1
		else
			raset%i=MOD(raset%ij/177, 177)+2
			raset%j=MOD(raset%ij,     177)+2
			k=MOD(raset%kl/169, 178)+1
			m=MOD(raset%kl,     169)
			do ii=1,97
				s=0.d0
				t=0.5d0
				do jj=1,24
					n=MOD(MOD(raset%i*raset%j,179)*k, 179)
					raset%i=raset%j
					raset%j=k
					k=n
					m=MOD(53*m+1, 169)
					if (MOD(m*n,64).GE.32) s=s+t
					t=0.5d0*t
				end do
				raset%u(ii)=s
			end do

			raset%c =  362436.d0/16777216.d0
			raset%cd= 7654321.d0/16777216.d0
			raset%cm=16777213.d0/16777216.d0
!
			raset%i=97
			raset%j=33
		end if
		end subroutine rmaset

!----------------------------------------------------------------------------

		subroutine rmaget()
		open(78,file="seedou.rng")
		write(78,*) raset
		close(78)
		end subroutine rmaget

!----------------------------------------------------------------------------

		subroutine ranmar(xr)
! pseudo random number generator
! proposed by Marsaglia, Zaman and Tsang
		double precision :: uni, xr
		uni=raset%u(raset%i)-raset%u(raset%j)
		if (uni.LT.0.d0) uni=uni+1.d0
		raset%u(raset%i)=uni
		raset%i=raset%i-1
		if (raset%i.EQ.0)raset%i=97
		raset%j=raset%j-1
		if (raset%j.EQ.0) raset%j=97
		raset%c=raset%c-raset%cd
		if (raset%c.LT.0.d0) raset%c=raset%c+raset%cm
		uni=uni-raset%c
		if (uni.LT.0.d0) uni=uni+1.d0
		xr=uni
		end subroutine ranmar

        end module Marsaglia

!============================================================================
subroutine inicial(L,N,M)
use Marsaglia
integer*8       ::  L,N,M(1:L,1:L),i,j
real*8          ::  xr


do i = 1,L
    do j = 1,L

        call ranmar(xr)		! Utilizacao do RNG
        if(xr <= 0.5) M(i,j) = -1
        if(xr >  0.5) M(i,j) = +1

    end do
end do



end subroutine inicial

subroutine calcula_energia(L,M,energia_Total,J_troca)
integer*8       ::  L,M(1:L,1:L),vizinho_esquerda,vizinho_direita,vizinho_cima,vizinho_baixo,i,j
real*8          ::  energia_Total,J_troca,soma

energiaTotal = 0d0
soma = 0d0

do i = 1,L
do j = 1,L

    vizinho_esquerda = j-1
    vizinho_direita  = j+1
    vizinho_cima     = i-1
    vizinho_baixo    = i+1

    if( i == 1) then
        vizinho_cima     = L
    end if

    if( i == L) then
        vizinho_baixo    = 1
    end if

    if( j == 1) then
        vizinho_esquerda = L
    end if

    if( j == L) then
        vizinho_direita  = 1
    end if


    soma = soma + M(i,j)*(  M(vizinho_baixo,j)  +   M(i,vizinho_direita) )

end do
end do

energia_total = (-1.0d0)*J_troca*soma

end subroutine calcula_energia

subroutine tenta_mover(M,T,L,J_troca,e,kb,energia_total,mag,energia_antiga,LNgE)
use Marsaglia
integer*8       ::  L,M(1:L,1:L),i,j,M_novo,vizinho_esquerda,vizinho_direita,vizinho_cima,vizinho_baixo,magnetizacao,mag
integer*8       ::  energia_nova,energia_antiga
real*8          ::  T,J_troca,xr,deltaE,temp,e,kb,peso,probabilidade,energia_total,LNgE(0:1024)

call ranmar(xr)
i = L*xr + 1

call ranmar(xr)
j = L*xr + 1


    vizinho_esquerda = j -1
    vizinho_direita  = j +1
    vizinho_cima     = i -1
    vizinho_baixo    = i +1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Condições de contorno periódicas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( i  == 1) then
        vizinho_cima     = L
    end if

    if( i == L) then
        vizinho_baixo    = 1
    end if

    if( j  == 1) then
        vizinho_esquerda = L
    end if

    if( j  == L) then
        vizinho_direita  = 1
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


temp = M(vizinho_baixo,j) + M(vizinho_cima,j) + M(i,vizinho_esquerda) + M(i,vizinho_direita)      ! Magnetização total dos vizinhos
deltaE = 2.0d0*J_troca*M(i,j)*temp

energia_nova = energia_antiga + deltaE/4


peso = e**( (-1.0d0)*   (  LNgE(energia_nova) - LNgE(energia_antiga)    )   )
probabilidade = min(1.0d0,peso)


call ranmar(xr)
if ( xr < probabilidade) then
    M(i,j) = (-1)*M(i,j)
    energia_total = energia_total + deltaE
    energia_antiga = energia_nova
    if( M(i,j) > 0) mag = mag + 2
    if( M(i,j) < 0) mag = mag - 2

end if

end subroutine tenta_mover




program spin
use Marsaglia

!Parameters
integer*8,parameter     ::   L        =  32
integer*8,parameter     ::   N        =  L*L
integer*8,parameter     ::   Mcs      =  10e6
real*8,parameter        ::   J_troca  =  1.0d0
real*8,parameter        ::   e        =  2.718281828
real*8,parameter        ::   kb       =  1.0d0



!Global Variables

real*8		:: beta,T,xr,energia_Total,deltaE,eMax,Energia_salva(0:Mcs),magnetizacao(0:Mcs),LNgE(0:1024)
integer*8   :: i,Hn,npassos,k,escolhida,alfa,M(1:L,1:L),Histograma_energia(N),Histograma_mag(2),mag,k1,lixo
integer*8   :: energia_antiga

open(1,file="exercicio2.dat")
open(9,file="LNgE.txt")

do i = 0,1024
    read(9,*)lixo,LNgE(k)
end do

call rmaset()	

energiaTotal = 0d0                                                                                
M = 1
call calcula_energia(L,M,eMax,J_troca)


call inicial(L,N,M)
call calcula_energia(L,M,energia_Total,J_troca)
mag              = sum(M)
magnetizacao(0)  = mag
energia_salva(0) = energia_Total

energia_antiga = (energia_total - emax) / 4

do k  = 1,Mcs
    do k1 = 1,N 
        call tenta_mover(M,T,L,J_troca,e,kb,energia_total,mag,energia_antiga,LNgE)
    end do

    magnetizacao( k)   = mag
    energia_salva(k)   = energia_total

end do

magnetizacao    = magnetizacao   
energia_salva   = energia_salva  

do k = 0, Mcs
    write(1,*)k,  energia_salva(k),    magnetizacao(k)
end do

end program spin

!============================================================================
