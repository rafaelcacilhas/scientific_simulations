program correlacao

!Parameters
integer*8,parameter     ::   terma    =  10e5
integer*8,parameter     ::   subinter =  100
integer*8,parameter     ::   itotal   =  10e6                 !numero de linhas no arquivo de entrada
integer*8,parameter     ::   ndata    =  itotal - terma
integer*8,parameter     ::   Na       =  ndata / subinter
integer*8,parameter     ::   nmedia   =  itotal
integer*8,parameter     ::   imax     =  ndata/10000
integer*8,parameter     ::   N        =  32*32

!Global Variables

real*8		:: energia(1:ndata),magnetizacao(1:ndata),soma2,lixo1,lixo2,lixo3,xi_j_energia(Na),xi_j_mag(Na),xi_energia,xi_mag
real*8      :: soma_xk_mag,soma_xk_energia,soma_xk_mag2,soma_xk_energia2,xj_energia(Na),xj_mag(Na),xk,xk2,xj_barra,soma_xj
real*8      :: sigma2_energia,sigma2_mag,sigma2_xi_energia,sigma2_xi_mag,cv,susc
real*8      :: soma_xj_energia,soma_xj_mag,soma_xi_energia,soma_xi,mag 
integer*8   :: i,k,s,media,l,l2,j

open(1,file="1.dat")                                        ! Arquivo de entrada
open(2,file="saida.dat")                                    ! Arquivo de saída

do i = 1,terma
    read(1,*)lixo1,lixo2,lixo3
end do

xj_energia          = 0d0
xj_mag              = 0d0
soma_xj_energia     = 0d0
soma_xj_mag         = 0d0
xk                  = 0d0

sigma2_energia      =   0d0
sigma2_mag          =   0d0
sigma2_xi_energia   =   0d0
sigma2_xi_mag       =   0d0

print*,ndata,subinter
do j = 1,Na

    soma_xk_energia      = 0d0
    soma_xk_energia2     = 0d0
    soma_xk_mag          = 0d0
    soma_xk_mag2         = 0d0


    do k = 1,subinter
    if (k == j) goto 32
 
       read(1,*)lixo1,energia(k),magnetizacao(k)

        soma_xk_energia  = soma_xk_energia  + energia(k)
        soma_xk_energia2 = soma_xk_energia2 + energia(k)*energia(k)
        soma_xk_mag      = soma_xk_mag      + magnetizacao(k) 
        soma_xk_mag2     = soma_xk_mag2     + magnetizacao(k)*magnetizacao(k)

    32 continue
    end do



xj_energia(j)      = soma_xk_energia / (ndata - subinter)
xj_mag(j)          = soma_xk_mag     / (ndata - subinter)

xi_j_energia(j)    = soma_xk_energia2  - soma_xk_energia 
xi_j_mag(j)        = soma_xk_mag2      - soma_xk_mag 


end do
 
x_barra_energia = xj_energia(Na)   / Na
x_barra_mag     = xj_mag(Na)       / Na

xi_energia      = xi_j_energia(Na) / Na
xi_mag          = xi_j_mag(Na)     / Na





do j = 1, Na

sigma2_energia      = sigma2_energia    +   (   xj_energia(j)   -  x_barra_energia  )**2
sigma2_mag          = sigma2_mag        +   (   xj_mag(j)       -  x_barra_mag      )**2
sigma2_xi_energia   = sigma2_xi_energia +   (   xi_j_energia(j) -  xi_energia       )**2
sigma2_xi_mag       = sigma2_xi_mag     +   (   xi_j_mag(j)     -  xi_mag           )**2


end do

sigma2_energia      =   (   (Na-1)  /Na )*sigma2_energia
sigma2_mag          =   (   (Na-1)  /Na )*sigma2_mag
sigma2_xi_energia   =   (   (Na-1)  /Na )*sigma2_xi_energia
sigma2_xi_mag       =   (   (Na-1)  /Na )*sigma2_xi_mag

 cv = T*xi_energia
 susc = 0d0

print*,"1.dat",x_barra_energia/64.0d0,x_barra_mag/64.0d0,

end program correlacao

!============================================================================
