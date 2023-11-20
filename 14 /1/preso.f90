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
subroutine inicial(L,P,pc,numero_coopera)
use Marsaglia
integer*8       ::  L,P(1:L),numero_coopera,j
real*8          ::  pc,xr

P = 0
numero_coopera = 0

do j = 1,L
    call ranmar(xr)		! Utilizacao do RNG
    if (xr < pc) then
        P(j) = 1
        numero_coopera = numero_coopera + 1
    end if
end do


end subroutine inicial



program preso
use Marsaglia

!Parameters
integer*8,parameter     ::   L        =  1000 !1000
integer*8,parameter     ::   M        =  1 !1000
integer*8,parameter     ::   N        =  400 !400
integer*8,parameter     ::   Nmax     =  500
integer*8,parameter     ::   z        =  16 ! 4,8,16
real*8,parameter        ::   T        =  1.6d0



!Global Variables

real*8		:: xr,ganho(1:L),melhor_ganho,pc,p_infinito,p_medio
integer*8   :: i,j,k,P(1:L),numero_coopera,tempo,vizinho,detento,melhor_vizinho,densidade,media,salva


open(1,file="z16.dat")
call rmaset()



do densidade = 8,8
    pc       = 0.1d0*densidade 
    p_medio  = 0d0
print*,z
    do media = 1,M

        call inicial(L,P,pc,numero_coopera)
        !print*,"Número de cooperadores: ",numero_coopera
        write(1,*)0,numero_coopera

        
        do tempo = 1,N

        ganho = 0d0
        maior_ganho = -1.0d0


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Determina os ganhos !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do detento = 1,L
                do k = -z/2,z/2
                    if (k == 0) goto 90

                    vizinho = detento + k
                    if ( vizinho < 1) vizinho = vizinho + L
                    if ( vizinho > L) vizinho = vizinho - L

                    ganho(detento) = ganho(detento) + P(detento)*P(vizinho) + T*( 1.0d0 - P(detento) )*P(vizinho) 


                90 continue        
                end do
            end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Determina melhores vizinhos !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do detento = 1,L
                melhor_vizinho = detento
                melhor_ganho   = 0d0

                do k = -z/2,z/2
                    if(k == 0) goto 92
                    vizinho = detento + k
                    if ( vizinho < 1) vizinho = vizinho + L
                    if ( vizinho > L) vizinho = vizinho - L

                    if( ganho(vizinho) > melhor_ganho  )  then
                        melhor_ganho = ganho(vizinho) 
                        melhor_vizinho = vizinho
                    end if

                92  continue    
                end do

                if( melhor_ganho > ganho(detento) .and. P(detento).ne. P(melhor_vizinho) ) then
                    P(detento) = P(melhor_vizinho)
                    if(P(detento) == 1) numero_coopera = numero_coopera  + 1
                    if(P(detento) == 0) numero_coopera = numero_coopera  - 1   
                end if
            


            end do


                write(1,*)tempo,numero_coopera





        end do

        p_infinito = numero_coopera / (1.0d0*L)
        p_medio    = p_medio + p_infinito

  
    end do


p_medio = p_medio /(1.0d0*M) 

print*,pc,p_medio



end do

end program preso

!============================================================================
