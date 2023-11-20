! This program  ...
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 2017.08.19
!
! INPUT: Data set with x, y and error (y)
! 
! OUTPUT: 



!xmgrace: 67.415 * x + 758.39

subroutine regression ( c_angular,c_linear,c_angular_2,N )

integer*4   ::  i,A,C,N
real*8      ::  y(35),x(35),B,D,c_angular,c_linear,c_angular_2,y_error(35)

A = 0
B = 0d0
C = 0
D = 0d0

open(1,file="hubble.dat")

do i=1,N
    read(1,*)x(i),y(i),y_error(i)
end do

close(1)

do i=1,N
    A = A + x(i)
    B = B + y(i)
    C = C + (   x(i)*x(i)   )
    D = D + (   x(i)*y(i)   )   

end do



c_angular = (N*D - A*B) /   (N*C - A*A)
c_linear =  (B*C - A*D) /   (N*C - A*A)

c_angular_2 = (N*D - A*A*D/C) / (N*C - A*A)     !c_angular_2 é o caso onde b = 0
print*,"Fit sem considerar erros em y:"
print*,"b diferente de 0"
print*, c_angular, " x + ",c_linear
print*,"Coeficiente angular quando b = 0: ",    c_angular_2
print*,""

end subroutine



subroutine fit(a,b,erro_a,erro_b,N)

integer*4   ::  i,N
real*8      ::  y(35),x(35),y_error(35),S,Sx,Sy,Sxx,Sxy,Delta,a,b,erro_a,erro_b,a_b_0,SS,erro_a_b_0


S = 0d0
Sx = 0d0
Sy = 0d0
Sxx = 0d0
Sxy = 0d0
Delta = 0d0
SS = 0d0

open(1,file="hubble.dat")

do i=1,N
    read(1,*)x(i),y(i),y_error(i)
end do

close(1)


do i=1,N
    
    S    =     1d0                        /   (y_error(i)**2)
    Sx   =     Sx + x(i) *   S
    Sy   =     Sy + y(i) *   S
    Sxx  =     Sxx +    x(i)*x(i)  *   S
    Sxy  =     Sxy +    x(i)*y(i)  *   S
    SS = SS + S

end do

Delta = SS*Sxx - (Sx*Sx)

a = (SS*Sxy  - Sx*Sy ) /   Delta
b = (Sxx*Sy - Sx*Sxy) /   Delta

erro_a = dsqrt(Sxx    /   Delta)
erro_b = dsqrt(SS     /   Delta)



a_b_0 = Sxy / Sxx
erro_a_b_0 = dsqrt(1d0/Sxx)

print*, "Fit consideranndo erros em y:"
print*,"b diferente de 0:"
print*, a, " x + ",b
print*, "Erro em a: ",erro_a, "; Erro em b: ",erro_b

print*,"Coeficiente angular forçando b = 0: ", a_b_0, "mais ou menos" ,erro_a_b_0




return
end subroutine fit 


program hubble
implicit none

integer*4   ::  i, ndata
real*8      ::  x(35),y(35),y_error(35),c_angular,c_linear,c_angular_2
real*8      ::  a,b,erro_a,erro_b

ndata = 35

call regression(c_angular,c_linear,c_angular_2,ndata)
call fit(a,b,erro_a,erro_b,ndata)



end program
