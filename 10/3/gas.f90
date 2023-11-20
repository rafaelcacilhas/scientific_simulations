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
	use Marsaglia
	integer*8            :: k,N
	real*8               :: r,vx(1:N),vy(1:N),phi,u,sigma,sigma2,x(1:N),y(1:N),kb,T,m,pi


    do i = 1,N
        call ranmar(u)		! Utilizacao do RNG
        call ranmar(phi)
        call ranmar(x(i))
        call ranmar(y(i))
        phi = 2.0d0*pi*phi

        x(i) = 100.0d0*x(i)
        y(i) = 100.0d0*y(i)

         r    = sigma*dsqrt( -2.0d0*dlog(1.0d0 - u)  )
        vx(i) = r    *dcos(phi) *0.01d0                         ! Velocidade em Ang/ps
        vy(i) = r    *dsin(phi) *0.01d0                         ! Velocidade em Ang/ps
        
       
    end do



end subroutine


subroutine histograma(x,y,vx,vy,N)
integer*8,parameter ::  ncaixas = 25
integer*8           ::  N,caixax,caixav,i,Hx(1:ncaixas),Hv(1:ncaixas)
real*8              ::  x(1:N),y(1:N),vx(1:N),vy(1:N),xmax,xmin,vmax,vmin,dHx,dHv

xmax = 0d0
xmin = 100d0
vmax = 0d0
vmin = 100d0


Hx = 0
Hv = 0

open(10,file="histogramax.dat")
open(11,file="histogramavx.dat")


do i = 1,N

    if (x(i)  > xmax) xmax =  x(i)
    if (x(i)  < xmin) xmin =  x(i)
    if (vx(i) > vmax) vmax = vx(i)
    if (vx(i) < vmin) vmin = vx(i)

end do

dHx = real(xmax - xmin)/real(ncaixas)
dHv = real(vmax - vmin)/real(ncaixas)


do i = 1,N

caixax = (x(i)  - xmin)/dHx
caixav = (vx(i) - vmin)/dHv

Hx(caixax) = Hx(caixax) + 1
Hv(caixav) = Hv(caixav) + 1

end do

do i = 1,ncaixas
write(10,*)xmax*i/real(ncaixas),Hx(i)
write(11,*)vmax*i/real(ncaixas),Hv(i)
end do

end subroutine histograma



	program gas
	use Marsaglia
	integer*8, parameter :: N     = 1000
	integer*8, parameter :: d     = 2
	integer*8, parameter :: L     = 100
	real*8, parameter    :: tmax  = 1000.0d0
	real*8, parameter    :: T     = 303
	real*8, parameter    :: m  = 1.0d0          !x10^-26
	real*8, parameter    :: kb = 1380.64852d0   !x10^-26
	real*8, parameter    :: pi = 3.14159265   
	real*8, parameter    :: dt = 1.0e-1   
	integer*8            :: k
	real*8               :: r,vx(1:N),vy(1:N),phi,u,sigma,sigma2,x(1:N),y(1:N),tempo,energia,mediavrms,mediav,v2

! Inicializacao do RNG
	call rmaset()	

	open(100,file="saida.dat")
	open(200,file="media.dat")

    x = 0d0
    tempo = 0d0
    energia = 0d0

sigma2 = kb*T/m
sigma = dsqrt(sigma2)


call inicial(sigma2,kb,T,m,sigma,u,phi,x,y,pi,r,vx,vy,N)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Inicial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Dinâmica !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call histograma(x,y,vx,vy,N)
mediavrms = 0d0
mediav    = 0d0

do while (tempo < tmax)
    energia = 0d0


    do i = 1,N

        x(i) = x(i) + vx(i)*dt
        y(i) = y(i) + vy(i)*dt

        if (x(i) >= L .or. y(i) >= L)   then

            call ranmar(u)		! Utilizacao do RNG
            call ranmar(phi)
            call ranmar(x(i))
            call ranmar(y(i))
            phi = 2.0d0*pi*phi

            x(i) = L*x(i)
            y(i) = L*y(i)

            r     = sigma*dsqrt( -2.0d0*dlog(1.0d0 - u)  )
            vx(i) = r    *dcos(phi) *0.01d0                         ! Velocidade em Ang/ps
            vy(i) = r    *dsin(phi) *0.01d0                         ! Velocidade em Ang/ps


        end if
    
       v2 = vx(i)*vx(i) + vy(i)*vy(i)
       energia = energia + m*v2/2.0d0
 
       mediavrms = mediavrms + v2
       mediav    = mediav    + dsqrt(v2)  

    end do

mediavrms = dsqrt(mediavrms /  real(N))
mediav    = mediav          /  real(N)

write(100,*)tempo,energia
write(200,*)tempo,mediavrms,mediav



tempo = tempo + dt

end do






	stop
	end program gas

!============================================================================
