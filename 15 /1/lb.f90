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


subroutine inicial(L,f,c,w,c0)
integer*8       ::  L,Lmin,Lmax,i
real*8          ::  f(0:L,0:L,0:8),c(1:2,0:8),w(0:8),c0


Lmin = 120
Lmax = 180

 c(1,0)   =   0d0;  c(1,1)   = c0   ;  c(1,2)   = 0d0   ;  c(1,3)   = -c0   ;  c(1,4)   =  0d0  ;  
 c(2,0)   =   0d0;  c(2,1)   = 0.0d0;  c(2,2)   = c0    ;  c(2,3)   =  0.0d0;  c(2,4)   = -c0   ;  

 c(1,5)   = c0  ;  c(1,6)   =-c0;  c(1,7)   =-c0;  c(1,8)   = c0;  
 c(2,5)   = c0  ;  c(2,6)   = c0;  c(2,7)   =-c0;  c(2,8)   =-c0;  


 w(0) = 4.0d0/9.0d0;    
 w(1) = 1.0d0/9.0d0 ;    w(2) = w(1);    w(3) = w(1);  w(4) = w(1)
 w(5) = 1.0d0/36.0d0;    w(6) = w(5);    w(7) = w(5);  w(8) = w(5)
    
f        = 0.0d0
f(:,:,0) = 0.4d0

do i = 120,180
    f(i,:,0)    =   0.8d0
end do



end subroutine inicial



subroutine parametros(L,f,pk,jk,c,uk)
integer*8       ::  L,i,j,k,componente
real*8          ::  f(0:L,0:L,0:8),pk(0:L,0:L),jk(1:2,0:L,0:L),c(1:2,0:8),uk(1:2,0:L,0:L),Forca

pk      = 0d0
jk      = 0d0
Forca   = 0d0


do k = 0,L
    do j = 0,L
        do i = 0,8

            pk(k,j)   = pk(k,j)            + f(k,j,i)
            jk(:,k,j) = jk(:,k,j)          + f(k,j,i)*c(:,i)
 
        end do


       uk(:,k,j) = (   jk(:,k,j) ) / pk(k,j)    

    end do
end do


end subroutine parametros





subroutine f_eq(L,f,w,c0,feq,pk,uk,c)
integer*8       ::  L,i,k1,k2
real*8          ::  f(0:L,0:L,0:8),feq(0:L,0:L,0:8),w(0:8),pk(0:L,0:L),c0,produto,uk(1:2,0:L,0:L),uk2,c(1:2,0:8),soma 

produto = 0d0

do k1 = 0,L
    do k2 = 0,L
    
        uk2 = sum(uk(:,k1,k2)*uk(:,k1,k2) )        


        do i = 0,8
	        produto = sum(  c(:,i)*uk(:,k1,k2) )
            feq(k1,k2,i) = w(i)*pk(k1,k2)*( 1.0d0 + 3.0d0*produto/ (c0**2) + 4.5d0*produto*produto / (c0**4) - 1.5d0*uk2 / (c0**2) )

        end do

    end do
end do


do k1=0,L
    do k2=0,l
        soma = soma + sum(feq(k1,k2,:))
    end do
end do



end subroutine f_eq




subroutine nova_f(L,f,feq,dt,tau,pk,f_novo)
integer*8       ::  L,i,j,k,c,viz1,viz2
real*8          ::  f(0:L,0:L,0:8),feq(0:L,0:L,0:8),f_novo(0:L,0:L,0:8),dt,tau,pk(0:L,0:L),uk(1:2,0:L,0:L),temp(0:L,0:L,0:8)

pk = 0d0




do i = 0,L
do j = 0,L
do c = 0,8

            if( c == 0) then

                viz1 = i
                viz2 = j

            end if

            if( c == 1) then

                viz1 = i + 1
                viz2 = j


            end if

            if( c == 2) then

                viz1 = i
                viz2 = j + 1

            end if

            if( c == 3) then

                viz1 = i
                viz2 = j -1

            end if

            if( c == 4) then

                viz1 = i - 1
                viz2 = j

            end if


            if( c == 5) then

                viz1 = i + 1
                viz2 = j + 1

            end if

            if( c == 6) then

                viz1 = i - 1
                viz2 = j + 1

            end if

            if( c == 7) then

                viz1 = i - 1
                viz2 = j - 1

            end if

            if( c == 8) then

                viz1 = i + 1
                viz2 = j - 1

            end if

            if ( viz1 > L   ) viz1 = viz1 - (L+1)
            if ( viz1 < 0   ) viz1 = viz1 + (L+1)
            if ( viz2 > L   ) viz2 = viz2 - (L+1)
            if ( viz2 < 0   ) viz2 = viz2 + (L+1)



    f_novo(viz1,viz2,c) = f(i,j,c) - (dt/tau)*( f(i,j,c) - feq(i,j,c)  )


end do  
end do
end do
    
f = f_novo



end subroutine nova_f



subroutine salva_perfil(L,pk,contador,uk)
integer*8       ::  L,i,nome,contador
real*8          ::  pk(0:L,0:L),uk(1:2,0:L,0:L)
character(len=10) :: file_id,file_id2
character(len=50) :: file_name,file_name2


    ! Saída de dados:
    nome = contador
    write(file_id, '(i0)') nome
    file_name = 'perfil' // trim(adjustl(file_id)) // '.dat'
    open(1,file = trim(file_name))

    do i = 0,200
        write(1,*)i,pk(i,100) 
    end do


    write(file_id, '(i0)') nome
    file_name = 'mapa' // trim(adjustl(file_id)) // '.dat'
    open(2,file = trim(file_name))


    do i = 0,200
        do j = 0,200
            write(2,*)i,j,pk(i,j)
        end do
    end do




    write(file_id, '(i0)') nome
    file_name = 'vetor' // trim(adjustl(file_id)) // '.dat'
    open(2,file = trim(file_name))

    do i = 0,200
        do j = 0,200
            write(2,*)i,j,i+(4.0d0*uk(1,i,j)),j+(4.0d0*uk(2,i,j))
        end do
    end do




end subroutine salva_perfil





program lb
use Marsaglia

!Parameters
integer*8,parameter     ::  L       =   200
integer*8,parameter     ::  vmax    =   20000
integer*8,parameter     ::  salva   =   50
real*8,parameter        ::  dt      =   1.0d0
real*8,parameter        ::  tau     =   1.0d0
real*8,parameter        ::  c0      =   1.0d0
real*8,parameter        ::  e       =   2.718281828459045235360287


!Global Variables

real*8		:: f(0:L,0:L,0:8),tempo,c(1:2,0:8),pk(0:L,0:L),jk(1:2,0:L,0:L),uk(1:2,0:L,0:L),w(0:8)
real*8      :: feq(0:L,0:L,0:8),f_novo(0:L,0:L,0:8)
integer*8   :: i,j,k,Lmin,Lmax,contador


call rmaset()



call inicial(L,f,c,w,c0)
tempo = 0d0

do contador = 1,vmax


    call    parametros(L,f,pk,jk,c,uk)
    
    if(contador < 500) then
        if( mod(contador,50) == 0   )       call    salva_perfil(L,pk,contador,uk)
    else if (contador < 1000) then
        if( mod(contador,100) == 0   )       call    salva_perfil(L,pk,contador,uk)
    else if (contador < 5000) then
        if( mod(contador,500) == 0   )       call    salva_perfil(L,pk,contador,uk)
    else if (contador < 20000) then
        if( mod(contador,2500) == 0   )       call    salva_perfil(L,pk,contador,uk)
    end if




    call    f_eq(L,f,w,c0,feq,pk,uk,c)
   
    call    nova_f(L,f,feq,dt,tau,pk,f_novo)



!    call    salva_mapa(L,pk)

tempo = tempo + dt
end do

1001 continue


end program lb

!============================================================================
