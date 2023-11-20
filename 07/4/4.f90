! This program solves a partial diferencial equation for the wave equation
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 28.08.2018
!
! INPUT:   
! 
! OUTPUT:   


module variables


!Parameters
real*8,parameter      ::  pi      = 4.0d0*datan(1.0d0) 
real*8,parameter      ::  J0      = 0.0d0
real*8,parameter      ::  J_inf   = 0.0d0
integer*8,parameter   ::  Nmax    = 200
integer*8,parameter   ::  M       = 100


!Global Variables


real*8      :: w,J1,soma,J(0:Nmax),t(0:Nmax),wmin,wmax,parte_real,parte_im,g_real,g_im
complex(8)  :: primeiro_termo,G,segundo_termo,terceiro_termo,quarto_termo,transformada,expoente,Jf
integer*8   :: i,k,contador

end module variables



program onda
use variables


    open(1,file="compliance.dat")
    do k = 0,Nmax
    read(1,*)t(k),J(k)
    end do
    close(1)

    open(2,file="saida.dat")
    open(3,file="g.dat")

wmin = 1.0d0/t(Nmax)
wmax = 1.0d0/t(0)

w = wmin
delta = (wmax/wmin)**(1.0d0/M)


do contador = 1,M
    expoente = dcmplx(0.0d0,-1.0d0*w)


    soma = 0d0
    do k = 2, Nmax
    soma = soma +   (  ( J(k) - J(k-1) ) * (  cdexp(expoente*t(k-1)   ) -   cdexp(expoente*t(k) )     )   / ( t(k) - t(k-1)  )     )
    end do


    primeiro_termo  =   dcmplx(0.0d0,w*J0) 
    segundo_termo   = ( 1.0d0   -   cdexp(expoente*t(1)) ) * ( J(1) - J(0) )/t(1)
    terceiro_termo  = J_inf*cdexp(expoente*t(Nmax))
    quarto_termo    = soma

    transformada = (-1.0d0/w**2)*(primeiro_termo + segundo_termo + terceiro_termo + quarto_termo ) 
    parte_real =   real(transformada)
    parte_im   =  aimag(transformada)
    

    g_real =  -1.0d0*parte_im  /(w*(parte_real**2 + parte_im**2))
    g_im   =  -1.0d0*parte_real/(w*(parte_real**2 + parte_im**2))


    write(2,*)w,dabs(parte_real),dabs(parte_im)
    write(3,*)w,g_real,g_im





w = w*delta
end do


end program


