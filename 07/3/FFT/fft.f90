! This program calculates a FFT like in DeVries book
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 21.08.2018
!
! INPUT:   
! 
! OUTPUT:   



subroutine four1(data,nn,isign)
integer     isign,nn
real        data(2*nn)
integer     i,istep,j,m,mmax,n
real        tempi,tempr
double precision    theta,wi,wpi,wpr,wr,wtemp

n = 2*nn
j = 1

do  i = 1,n,2
    if (j .gt. i) then
        tempr = data(j)
        tempi = data(j+1)
        data(j) = data(i)
        data(j+1) = data(i+1)
        data(i) = tempr
        data(i+1) = tempi
    end if
    m = n/2

1   if ((m.ge.2) .and. (j.gt.m) ) then
        j = j-m
        m = m/2

    goto 1
    end if
    j = j+m

end do
mmax = 2

2   if(n.gt.mmax) then

        istep = 2*mmax
        theta = 6.283185/(isign*mmax)
        wpr = -2.d0*sin(0.5d0*theta)**2
        wpi = sin(theta)
        wr = 1.d0
        wi = 0d0

        do  m=1,mmax,2
            do  i=m,n,istep

                j = i+mmax
                tempr = sngl(wr)*data(j)   - sngl(wi)*data(j+1)
                tempi = sngl(wr)*data(j+1) + sngl(wi)*data(j)
                data(j)     = data(i)   - tempr
                data(j+1)   = data(i+1) - tempi
                data(i)     = data(i)   + tempr
                data(i+1)   = data(i+1) + tempi

            end do
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi

       end do
            
       mmax = istep
            
goto 2
end if

return
end






program fast

integer,parameter   ::  isign = 1
integer,parameter   ::  nn = 1024

real        data(1:2*nn),tempo(1:2*nn)
data = 0d0

open(1,file="serie2.dat")

do i = 1,nn
read(1,*)tempo(i),data(i)
end do

print*,tempo(nn),data(nn)

call four1(data,nn,isign)



open(2,file="saida2.dat")

do i = 1,nn
data(i) = data(i) **2
write(2,*)tempo(i),data(i)
end do

w0 = 2.0d0*pi / 3.255208333333333E-004
print*,w0*tempo(nn),data(nn)

end program

	
