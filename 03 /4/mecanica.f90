! This program solves a triple harmonic oscillator 
! Author: R. V. Cacilhas
! E­mail: rafael.cacilhas at ufv.br
! Last Update: 03.08.2018
!
! INPUT:    
! 
! OUTPUT:   





module variables

!Parameters
real*8      ::  w0 =    0.4d0


!Global Variables

real*8      :: H(1:3,1:3),u1(1:3),u2(1:3),u3(1:3),funcao,raiz,derivada,erro,raizes(1:3),v_inicial(1:3),v_inicial_normalizado(1:3)
real*8		:: l1,l2,l3,a,l1_exato,l2_exato,l3_exato,lm
real*8      :: soma = 0d0
integer*8   :: i,j,nraizes,k

end module variables



subroutine newton(raizes,nraizes)

real*8      ::  x0 = -0.0d0
real*8      ::  deltaX = 0.1
real*8      ::  epslon = 10.0E-5
real*8      ::  erroM = 10e-10
real*8      ::  w0 = 0.4d0
real*8      ::  xi = 0.0, xf


integer*8   ::  i, npassos,j,nraizes,contador
real*4      ::  teste1,teste2
real*8      ::  funcao,derivada,xmin,xmax,raiz,raizes(1:3),erros(1:3),seeds(1:3),alfa,temp(1:3) = 0d0

deltaX = dsqrt(2.0d0)*(w0**2)
xi = -1.0
xf = 2*(w0**2)
x0 = 0
npassos = 200
raizes = 9999                                                               !Defini este valor para que ao analilsar o array raizes o valor padrão 0 não seja confundido com uma raiz.
nraizes = 0


do i=1,npassos,1
x = x0

    do

        if(x == 0.0d0) then
            funcao =          (x**3) - 2.0d0*(w0**4)*x
            derivada =  3.0d0*(x**2) - 2.0d0*(w0**4)

            if(funcao == 0d0) then
                raiz = x
                erro = 0.0d0   
                goto 1
            else 
                x = x + deltaX
            end if
 
        end if

        funcao =    (      (x**3) - 2.0d0*(w0**4)*x)                   !A função original assume valores
        derivada =  (3.0d0*(x**2) - 2.0d0*(w0**4))

        raiz = x
        x = x - funcao/derivada


        erro = (raiz-x)/x

        if(erro < 0) then
            erro = -1.0d0*erro
        end if






1 continue
        if ( erro < erroM) then     

            do j = 1,nraizes,1                                  !Checa se a raiz encontrada ja foi encontrada anteriormente
                teste1 = raiz
                teste2 = raizes(j)

                if(teste1 == teste2) then
                    goto 32
                 endif
            end do

            if(x < xi .or. x > xf) then                         !Checa se a raiz encontrada está no intervalo de interesse
                goto 32
            end if
            
            nraizes = nraizes + 1  
            raizes(nraizes) = x
            erros(nraizes) = erro
            seeds(nraizes) = x0
!            print*,"Raiz: ",raizes(nraizes),"Erro: ", erros(nraizes), "X inicial: ", seeds(nraizes)


            32 continue
            exit
        end if


        
    end do


x0 = x0 + deltaX
end do



return
end subroutine newton


program mecanica
use variables



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definir as matrizes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

H(1,1)  =  2*(w0**2)   ;   H(1,2)  = -(w0**2)       ;   H(1,3)  = 0.0d0             ;
H(2,1)  = -(w0**2)     ;   H(2,2)  = 2*(w0**2)      ;   H(2,3)  = -(w0**2)          ;
H(3,1)  = 0.0d0        ;   H(3,2)  =  -(w0**2)      ;   H(3,3)  = 2*(w0**2)         ;


u1(1) = 1.0d0;  u1(2) = dSqrt(2.0d0);    u1(3) =  1.0d0;
u2(1) = 1.0d0;  u2(2) = 0d0;             u2(3) = -1.0d0;
u3(1) = 1.0d0;  u3(2) = -dSqrt(2.0d0);   u3(3) =  1.0d0;




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Solução do sistema  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call newton(raizes,nraizes)         !Não usaremos a equação em lambda, porque é excessivamente complicada. Usaremos 2w0² - lambda = alfa e resolveremos para alfa
raizes(3) = -1.0d0*(raizes(2))      !alfa é quadratico, então possui duas raizes com sinais diferentes. Note que raizes(3) não foi encontrada.)


l1 = 2.0d0*(w0**2) - raizes(2)  ;   l1_exato = (w0**2)*(2.0d0-dsqrt(2.0d0))             !trocamos raizes(2) por raizes(1) para ficar na mesma ordem que o resultado teorico.
l2 = 2.0d0*(w0**2) - raizes(1)  ;   l2_exato = (w0**2)*(2.0d0)
l3 = 2.0d0*(w0**2) - raizes(3)  ;   l3_exato = (w0**2)*(2.0d0+dsqrt(2.0d0))



print*,"Autovalores encontrados: ",l1,l2,l3
print*,"Autovalores exatos:      ",l1_exato,l2_exato,l3_exato
print*,


v_inicial(1) = 1.0d0;   v_inicial(2) = 1.0d0;   v_inicial(3) = 1.0d0;
v_inicial_normalizado = v_inicial / dsqrt(v_inicial(1)**2 + v_inicial(2)**2 + v_inicial(3)**2) 


do k=1,10

v_inicial = matmul(H,v_inicial)


lm = v_inicial(1)*v_inicial_normalizado(1) + v_inicial(2)*v_inicial_normalizado(2) + v_inicial(3)*v_inicial_normalizado(3) 

v_inicial_normalizado = v_inicial / dsqrt(v_inicial(1)**2 + v_inicial(2)**2 + v_inicial(3)**2) 

end do

print*,"Autovetor (normalizado) encontrado: ", v_inicial_normalizado(1),v_inicial_normalizado(2),v_inicial_normalizado(3)
print*,"Autovetor (normalizado) esperado: ", 1.0d0/2.0d0,-dsqrt(2.0d0)/2.0d0,1.0d0/2.0d0
print*,

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REPRODUÇÃO DO LIVRO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
v_inicial(1) = 2.279d0;   v_inicial(2) = 4.4710d0;
v_inicial_normalizado(1) = 0.464d0;   v_inicial_normalizado(2) = 0.886d0;

lm = v_inicial(1)*v_inicial_normalizado(1) + v_inicial(2)*v_inicial_normalizado(2)

print*,"Valor encontrado: ",lm, "Valor esperado: ", 5.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end program


