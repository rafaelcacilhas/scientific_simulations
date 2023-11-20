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

subroutine inicial(sigma2,kb,T,m,sigma,u,phi,x,y,pi,r,vx,vy,N)
	integer*8            :: k,N
	real*8               :: r,vx(1:N),vy(1:N),phi,u,sigma,sigma2,x(1:N),y(1:N),kb,T,m,pi




end subroutine


subroutine histograma(x,y,vx,vy,N)
integer*8       ::  N,ncaixas,caixa
real*8          ::  x(1:N),y(1:N),vx(1:N),vy(1:N),H(1:10)



end subroutine histograma



	program direto
	use Marsaglia
	integer*8, parameter :: N     = 599
	integer*8, parameter :: M     = 10000000
	integer*8, parameter :: M1    = 5e4
	integer*8, parameter :: M2    = 5e6
	integer*8, parameter :: i01   = 534
	integer*8, parameter :: i02   = 535
	integer*8, parameter :: salvo = 1000
    integer*8            :: i,escolhida,pi(1:10,N),particula(0:N),Na,Nb,j,H(N),SomaH(1:10,N)
	real*8               :: sorteio,mediaNa,mediaNb

! Inicializacao do RNG
	call rmaset()	

	open(100,file="saida.dat")
    particula = 0
    pi = 0
    mediaNa = 0
    mediaNb = 0


do j = 1, M1

Na    = 0
Nb    = 0
H     = 0
SomaH = 0

!!! Condição Inicial !!!!
    do i = 1,i01    
        particula(i) = +1
        Na = Na + 1
    end do
    do i = i01+1, N
        particula(i) = -1
        Nb = Nb + 1
    end do

    H(Na) = H(Na) + 1 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    do i = 1,10000

        call ranmar(sorteio)
        sorteio = N*sorteio
        escolhida = dabs( sorteio )

        if( particula(escolhida) == 1 ) then                    !!!!!! Atualiza o número de particulas dependendo do sorteio
            Na = Na - 1
            Nb = Nb + 1
        else if (particula(escolhida) == -1) then
            Na = Na + 1
            Nb = Nb - 1
        end if

    

        H(Na) = H(Na) + 1                                       !!! Atualiza o histograma
        particula(escolhida) = particula(escolhida)*(-1)        !!! Transfere a particula



        if(mod(i,salvo) == 0) then
            SomaH((i/1000),Na) = SomaH((i/1000),Na) + 1         !!! Salva o resultado caso i seja multiplo de 1000
        end if


    end do


    mediaNa = mediaNa + Na
    mediaNb = mediaNb + Nb

!print*,Na,Nb
end do

mediaNa = real(mediaNa) / real(M1)
mediaNb = real(mediaNb) / real(M1)

pi = SomaH/M1


	stop
	end program direto

!============================================================================
