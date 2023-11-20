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


subroutine inicia(maior_n,Universo,x,y,nvizinhos,n,possiveis_lugares,x_possivel,y_possivel)
integer*8       ::  maior_n,Universo(0:200,0:200),x(1:Nmax),y(1:Nmax),nvizinhos(1:Nmax),n,possiveis_lugares
integer*8       ::  x_possivel(1:200),y_possivel(1:200)

maior_n = 0
Universo = 0
Universo(100,100) = 1
Universo(101,100) = 1
x(1) = 100; y(1) = 100
x(2) = 101; y(2) = 100

n = 2
maior_n = n




end subroutine inicia


subroutine probabilidade(possiveis_lugares,nvizinhos,Nmax,s,psi,kmais,kmenos,e,n)
use Marsaglia
integer*8   ::  Nmax,nvizinhos(1:Nmax),possiveis_lugares,n,escolhido
real*8      ::  s,kmais,kmenos(1:Nmax),norma,e,psi,soma,sorteio

norma   = 0d0

kmais   = e**(s)

do i = 1,n
    kmenos(i) = e**(    2.0d0*( 2.0d0 - nvizinhos(i)    )*psi   )   
    norma = norma + kmenos(i)
end do

norma = norma + possiveis_lugares*kmais

kmais   = kmais     / norma
kmenos  = kmenos    / norma



end subroutine probabilidade





subroutine decide(escolhido,possiveis_lugares,Universo,x_possivel,y_possivel,n,kmais,kmenos,nmax,sorteio,x,y,nvizinhos)
use Marsaglia
integer*8       ::  escolhido,possiveis_lugares,Universo(0:200,0:200),x_possivel(1:200),y_possivel(1:200),n,nmax,x(1:Nmax),y(1:Nmax)
integer*8       ::  nvizinhos(1:Nmax)
real*8          ::  kmais,kmenos(1:nmax),sorteio,prob,norma
norma = 0d0
prob  = 0d0

call ranmar(sorteio)		! Utilizacao do RNG



if( sorteio < possiveis_lugares*kmais ) then                                                         !   Cria partícula
!print*,"cria",sorteio,"menor que",possiveis_lugares*kmais,possiveis_lugares
    call ranmar(sorteio)		                                                                     !   Sorteia qual partícula será criada entre as (possiveis_lugares) possibilidades
    escolhido = sorteio*possiveis_lugares + 1

    n = n + 1
    x(n) = x_possivel(escolhido)
    y(n) = y_possivel(escolhido)    
    Universo(   x(n)    ,   y(n)    ) = 1

else if ( sorteio > possiveis_lugares*kmais) then                                                   !   Destrói partícula
!print*,"destroi",sorteio,"maior que",possiveis_lugares*kmais,possiveis_lugares
    do i = 1,N
        norma = norma + kmenos(i)
    end do
    kmenos = kmenos / norma


    call ranmar(sorteio)

    do i = 1,N
        prob = prob + kmenos(i)

        if(sorteio < prob) then
            escolhido = i
            goto 1002
        end if

    end do

    1002 continue

    x(escolhido)            = x(n)
    y(escolhido)            = y(n)
    nvizinhos(escolhido)    = nvizinhos(n)
    Universo(x(escolhido),y(escolhido)  )    = 0
    n = n-1     

end if

end subroutine decide


subroutine avalia_situacao(n,universo,x,y,nvizinhos,nMax,e,psi,possiveis_lugares,x_possivel,y_possivel)
integer*8       ::  i,n,nmax,Universo(0:200,0:200),x(1:Nmax),y(1:Nmax),nvizinhos(1:Nmax),possiveis_lugares
integer*8       ::  x_possivel(1:200),y_possivel(1:200)
real*8          ::  kmenos(1:nMax),e,psi

nvizinhos = 0
possiveis_lugares = 0
x_possivel = 0
y_possivel = 0

do i = 1,n

if( x(i)  == Nmax  )   goto 234
    if      ( Universo(   x(i) + 1,y(i)       )   ==  1) then
        nvizinhos(i) = nvizinhos(i) + 1
    else if ( Universo(   x(i) + 1,y(i)       )   ==  0) then
        possiveis_lugares = possiveis_lugares + 1

        x_possivel(possiveis_lugares) = x(i) + 1
        y_possivel(possiveis_lugares) = y(i)
!        print*,"Direita"
    end if
234 continue

if( x(i)  == 1  )   goto 123
    if      ( Universo(   x(i) - 1,y(i)       )   ==  1) then
        nvizinhos(i) = nvizinhos(i) + 1
    else if ( Universo(   x(i) - 1,y(i)       )   ==  0) then
        possiveis_lugares = possiveis_lugares + 1

        x_possivel(possiveis_lugares) = x(i) - 1
        y_possivel(possiveis_lugares) = y(i)
!        print*,"Esquerda"
    end if
123 continue

if( y(i)  == Nmax  )   goto 789
    if      ( Universo(   x(i)    ,y(i) + 1    )   ==  1) then
        nvizinhos(i) = nvizinhos(i) + 1
    else if ( Universo(   x(i)    ,y(i) + 1    )   ==  0) then
        possiveis_lugares = possiveis_lugares + 1

        x_possivel(possiveis_lugares) = x(i)
        y_possivel(possiveis_lugares) = y(i) + 1
!        print*,"Cima"
    end if
789 continue


if( y(i)  == 1  )   goto 987
    if      ( Universo(   x(i)    ,y(i) - 1    )   ==  1) then
        nvizinhos(i) = nvizinhos(i) + 1
    else if ( Universo(   x(i)    ,y(i) - 1    )   ==  0) then
        possiveis_lugares = possiveis_lugares + 1

        x_possivel(possiveis_lugares) = x(i)
        y_possivel(possiveis_lugares) = y(i) - 1
!        print*,"Baixo"
    end if
987 continue

end do



end subroutine avalia_situacao



program kmc
use Marsaglia

!Parameters
integer*8,parameter     ::  N0      =   100000
integer*8,parameter     ::  n_ini   =   2
integer*8,parameter     ::  Nmax    =   200
real*8,parameter        ::  Psi     =   1.0d0
real*8,parameter        ::  e       =   2.718281828459045235360287


!Global Variables

real*8		:: xr,s,P(0:Nmax),kmais,kmenos(1:Nmax),norma,sorteio,prob,P3
integer*8   :: i,j,k,maior_n,Numero(0:Nmax),n,simulacao,Universo(0:200,0:200),x(1:Nmax),y(1:Nmax),nvizinhos(1:Nmax),stopp
integer*8   :: saturacao,possiveis_lugares,x_possivel(1:200),y_possivel(1:200),escolhido,passo,maior


call rmaset()
stopp = 2
P = 0d0
Numero = 0
maior = 0


    saturacao = 25
    s = 0.1d0 + ( saturacao - 1.0d0)*0.1d0
  open(1,file="s25.dat")
  open(2,file="n.dat")

    do simulacao = 1,N0

        call inicia(maior_n,Universo,x,y,nvizinhos,n,possiveis_lugares,x_possivel,y_possivel)
        call avalia_situacao(n,universo,x,y,nvizinhos,nmax,e,psi,possiveis_lugares,x_possivel,y_possivel)



            do passo = 1,1000
            

                call probabilidade(possiveis_lugares,nvizinhos,Nmax,s,psi,kmais,kmenos,e,n)

                call decide(escolhido,possiveis_lugares,Universo,x_possivel,y_possivel,n,kmais,kmenos,nmax,sorteio,x,y,nvizinhos)

                call avalia_situacao(n,universo,x,y,nvizinhos,nmax,e,psi,possiveis_lugares,x_possivel,y_possivel)


            if  ( n > maior_n) maior_n = n        
            if  ( n == 0     ) goto 999
            if  ( n >= 200   ) goto 999

            end do
            999 continue



        do j = 1,maior_n
                Numero(j) = Numero(j) + 1
        end do

    if  ( maior_n > maior) maior = maior_n        
    end do


    do k = 1,maior
        P(k) = 1.0d0*Numero(maior) / (   1.0d0*Numero(k)    )
        write(1,*)k,P(k)
        write(2,*)k,Numero(k)
    end do

    P3 = 1.0d0*Numero(maior) / (1.0d0*N0)
    print*,s,P3

!do i = 1,200
!do j = 1,200

!if(Universo(i,j) == 1) write(1,*)i,j

!end do
!end do


end program kmc

!============================================================================
