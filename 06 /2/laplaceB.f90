! This program solves a partial diferencial equation for the Laplacian Equation
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 21.08.2018
!
! INPUT:   
! 
! OUTPUT:   


module variables


!Parameters
real*8,parameter      ::  pi      = 4.0d0*datan(1.0d0) 
real*8,parameter      ::   L      = 1.0d0  
real*8,parameter      ::   C      = 0.4d0*L 
real*8,parameter      ::   Vc     = 0.0d0  
real*8,parameter      ::   Ve     = -100.0d0  
real*8,parameter      ::   h      = 0.005d0 
real*8,parameter      ::   epslon = 10e-9
integer*8,parameter   ::   Nmax   = L / h 



!Global Variables

real*8		:: r,dx,dy,lambda,x,y,V(0:Nmax,0:Nmax),erromax,erro
integer*8   :: linha,coluna,i,j,k

end module variables


subroutine atualiza(V,xmin,xmax,ymin,ymax,h,Vc,epslon,Nmax)

integer*8   :: linha,coluna,i,j,k,Nmax
real*8		:: r,dx,dy,lambda,x,y,V(0:Nmax,0:Nmax),erromax,erro,Vc,h,epslon

do k = 1,100000
erromax = 0d0

    y = h
    do linha = 1,Nmax-1
    x = h

        do coluna = 1,Nmax-1

            if ( x > xmin .and. x < xmax .and. y > ymin .and. y < ymax) then
                V(linha,coluna) = Vc
 
            else 
                erro = V(linha,coluna)          
                V(linha,coluna) =  0.25d0*( V(linha+1,coluna) + V(linha,coluna+1) + V(linha-1,coluna) + V(linha,coluna-1)     )
                erro = dabs(V(linha,coluna) - erro)

                if(erro > erromax) then             !Calcula o maior erro para cada k
                    erromax = erro
                end if

            end if


        
        10 continue
        x = x + h
        end do

    y = y + h

    end do

if(erromax < epslon) then                        !Se o MAIOR erro é menor do que o criterio o programa finaliza
goto 20
end if


end do

20 continue



end subroutine


subroutine saida(V,Nmax,h)
integer*8   :: i,j,Nmax
real*8      :: V(0:Nmax,0:Nmax),h


do j = 0,Nmax
    do i = 0,Nmax
            write(2,*)i*h,j*h,V(i,j)
    end do
end do


end subroutine

program laplace
use variables

open(2,file="data.dat")

x = h
y = h

V = 0.0d0
V(0,:) = Ve
V(Nmax,:) = Ve
V(:,0) = Ve
V(:,Nmax) = Ve

xmin = (L-C)/2.0d0
ymin = xmin

xmax = (L+C)/2.0d0
ymax = xmax

call atualiza(V,xmin,xmax,ymin,ymax,h,Vc,epslon,Nmax)
call saida(V,Nmax,h)

!Usar setsizesquared

close(2)
end program


