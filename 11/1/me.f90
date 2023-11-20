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
subroutine inicial(N,sigma,L,rc,x,y)
use Marsaglia
integer*8       ::  N,i,j
real*8          ::  sigma,L,rc,x(1:N),y(1:N),xr,x_rascunho,y_rascunho,r,r_rascunho,deltax,deltay


do i = 1,N
1000 continue

    call ranmar(xr)		! Utilizacao do RNG
    x(i) = xr*L
    call ranmar(xr)		! Utilizacao do RNG
    y(i) = xr*L


    do j = 1,i-1                                                                        !!! Verifica se todas estão com distancias superiores a rc


        deltax = x(j) - x(i)
        deltay = y(j) - y(i)

        if(deltax >        L/2.0d0) deltax = deltax - L
        if(deltax < -1.0d0*L/2.0d0) deltax = deltax + L

        r = dsqrt(  deltax**2 + deltay**2   )

        if(r < rc)                   goto 1000

    end do

end do


end subroutine inicial

subroutine calcula_energia(x,y,N,rc,epslon,sigma,energiaTotal,L)
integer*8   ::  N,i,j
real*8      ::  x(1:N),y(1:N),energia,r,energiaTotal,rc,epslon,sigma,deltax,deltay,L

energiaTotal = 0d0

do i = 1,N

    do j = i+1,N                                                                        !!! Verifica se todas estão com distancias superiores a rc

        deltax = x(j) - x(i)
        deltay = y(j) - y(i)

        if(deltax >        L/2.0d0) deltax = deltax - L
        if(deltax < -1.0d0*L/2.0d0) deltax = deltax + L

        r = dsqrt(  deltax**2 + deltay**2   )

        if(r > rc) energia = 0d0
        if(r < rc) energia =  4.0d0*epslon*(      (sigma/r)**12   -   (sigma/r)**6    +   0.25d0      )

        energiaTotal = energiaTotal + energia

    end do



end do

end subroutine calcula_energia

subroutine tenta_mover(x,y,N,energianova,xnovo,ynovo,escolhida,rc,x_aux,y_aux)
use Marsaglia

integer*8               ::  N,escolhida,k
real*8                  ::  x(1:N),y(1:N),dx,dy,sorteio,xnovo,ynovo,deltax,deltay,energia,energianova
real*8                  ::  energiatotal,rc,probabilidade,random,x_aux(1:N),y_aux(1:N)
real*8                  ::   dmax
dmax = 1.5d0*rc


call ranmar(sorteio)		! Utilizacao do RNG
sorteio = sorteio*N + 1.0d0
escolhida = sorteio

call ranmar(dx)		! Utilizacao do RNG
dx = dmax*dx

call ranmar(dy)
dy = dmax*dy

x_aux = x  
y_aux = y
x_aux(escolhida) = x(escolhida) + dx
y_aux(escolhida) = y(escolhida) + dy


end subroutine tenta_mover

subroutine energiaLocal(x,y,N,x_aux,y_aux,escolhida,L,deltaE,rc,sigma,epslon)
integer*8       ::  N,escolhida,j
real*8          ::  x(1:N),x_aux(1:N),y(1:N),y_aux(1:N),deltax,deltay,L,energianova,energiavelha,deltaE
real*8          ::  energianovaTotal,energiavelhaTotal,rc,sigma,epslon

energianovaTotal  = 0d0
energiavelhaTotal = 0d0

do j = 1,N
if(j == escolhida) goto 321
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Energia antiga !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
        deltax = x(j) - x(escolhida)
        deltay = y(j) - y(escolhida)

        if(deltax >        L/2.0d0) deltax = deltax - L
        if(deltax < -1.0d0*L/2.0d0) deltax = deltax + L

        r = dsqrt(  deltax**2 + deltay**2   )

        if(r > rc) energiavelha = 0d0
        if(r < rc) energiavelha = 4.0d0*epslon*(      (sigma/r)**12   -   (sigma/r)**6    +   0.25d0      )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Energia nova    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
        deltax = x(j) - x_aux(escolhida)
        deltay = y(j) - y_aux(escolhida)

        if(deltax >        L/2.0d0) deltax = deltax - L
        if(deltax < -1.0d0*L/2.0d0) deltax = deltax + L

        r = dsqrt(  deltax**2 + deltay**2   )

        if(r > rc) energianova = 0d0
        if(r < rc) energianova = 4.0d0*epslon*(      (sigma/r)**12   -   (sigma/r)**6    +   0.25d0      )


        energiavelhaTotal = energiavelhaTotal + energiavelha
        energianovaTotal  = energianovaTotal  + energiaNova 

321 continue
end do

deltaE = energianovaTotal - energiavelhaTotal

end subroutine energiaLocal

program estatistica
use Marsaglia

!Parameters
integer*8,parameter     ::   N        =  200
integer*8,parameter     ::   Mcs      =  100
real*8,parameter        ::   dH       =  1.0d0
real*8,parameter        ::   mp       =  10.0d0
real*8,parameter        ::   sigma    =  1.0d0
real*8,parameter        ::   L        =  28.0d0*sigma
real*8,parameter        ::   epslon   =  0.16021773             !Eletron-volts
real*8,parameter        ::   rc       =  (  2.0d0**(1.0d0/6.0d0)  )*sigma
real*8,parameter        ::   tn       =  1.0
real*8,parameter        ::   dt       =  0.6d0
real*8,parameter        ::   e        =  2.718281828
real*8,parameter        ::   kb       =  1.38064852E-2
real*8,parameter        ::   T        =  (273.0d0 - 15.0d0)                    !30ºC em K
real*8,parameter        ::   beta     =  1.0d0/(kb*T)



!Global Variables

real*8		:: x(1:N),y(1:N),vx(1:N),vy(1:N),vi,vf,v2(1:N),xr, energiaTotal,energiaNova,random
real*8      :: xnovo,ynovo,deltaE,x_aux(1:N),y_aux(1:N),salvaEnergia(1:Mcs),energiaTotalAntiga
integer*8   :: i,Hn,npassos,k,escolhida

open(1,file="d=15.dat")
call rmaset()	
energiaTotal = 0d0


call inicial(N,sigma,L,rc,x,y)
call calcula_energia(x,y,N,rc,epslon,sigma,energiaTotal,L)
write(1,*)0,energiaTotal

npassos = N*Mcs
do k = 1,npassos
energiaNova  = 0d0
 
    call tenta_mover(x,y,N,energianova,xnovo,ynovo,escolhida,rc,x_aux,y_aux)

!    call calcula_energia(x_aux,y_aux,N,rc,epslon,sigma,energiaNova,L)                              !Energia Global (Caso queira a energia por este metodo, é necessario 
!    deltaE = energiaNova - energiaTotal                                                            !                descomentar esta linha  também   

    call energiaLocal(x,y,N,x_aux,y_aux,escolhida,L,deltaE,rc,sigma,epslon)                         !Energia Local

    peso = e**( -1.0d0*beta*deltaE        )
    probabilidade = min(1.0d0,peso)
    call ranmar(random)
    
    if(random < probabilidade) then
        x = x_aux
        y = y_aux
        energiaTotal = energiaTotal + deltaE
    end if



    


    if(mod(k,N) == 0) then
        salvaEnergia(  k/N   ) = energiaTotal
        print*,(k/N),salvaEnergia(k/N),deltaE
    end if

end do



do i = 1,Mcs
    write(1,*)i,salvaEnergia(i)

end do


end program estatistica

!============================================================================
