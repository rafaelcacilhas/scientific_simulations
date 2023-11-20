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

	program direto
	use Marsaglia
	integer*8, parameter :: ni = 1
	integer*8, parameter :: nf = 100
	integer*8, parameter :: nm = 50
	integer*8            :: k
	real*8               :: xr

! Inicializacao do RNG
	call rmaset()	

	write(*,*) "============================================================"
	write(*,*) "RNG init param.:", Nseed, raset%ij, raset%kl
	write(*,*) "============================================================"






	open(100,file="sequencia_direta.dat")

	do k=ni,nf
		call ranmar(xr)		! Utilizacao do RNG
		write(100,*) k, xr
	end do
	close(100)

	call rmaget()	! Salva estado do RNG no arquivo "seedin.rng"



	stop
	end program direto

!============================================================================
