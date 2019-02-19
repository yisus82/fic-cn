program nolineal
implicit none
logical :: converged
integer :: i,j,opcion,n
real (kind=8),dimension(3,3) :: m
real (kind=8),dimension(3) :: b
real (kind=8) :: aux1,aux2,x,y,z,xold,yold,zold,error

print*,'1.- Iteracion funcional'
print*,'2.- Newton-Raphson'
read*,opcion
x=0_8
y=0_8
z=0_8
xold=0_8
yold=0_8
zold=0_8


select case(opcion)
  case(1)
       call itfunc
       call txtsalida

  case(2)
     do n=1,100
       call generasistema
       call relajacion 

       aux1 = (xold-x)**2 + (yold-y)**2 + (zold-z)**2
       aux2 =  x**2 + y**2 + z**2
       error =  sqrt(aux1) / sqrt(aux2)

       write (6,600) n,x,y,z,error
       converged = (  error < 1e-9)
       if (converged) exit
       xold=x
       yold=y
       zold=z
     end do
     call txtsalida

  case default
end select
600 format (i5,5(2x,f16.8))

contains

subroutine itfunc
implicit none
real(kind=8) :: rx,ry,rz
do n=1,100
    rx=G1(y,z)
    ry=G2(x,z)
    rz=G3(x,y)
    x=rx
    y=ry
    z=rz
    aux1=(xold-x)**2+(yold-y)**2+(zold-z)**2
    aux2=x**2+y**2+z**2
    converged = (  sqrt(aux1) / sqrt(aux2) < 1e-9)
    if (converged) exit
    xold=x
    yold=y
    zold=z
    print*, 'ITERACIONES=',n
    print *,'APROXIMACION=',x,y,z
end do
end subroutine itfunc


subroutine generasistema
implicit none
m(1,1)=DF1x(x,y,z) 
m(1,2)=DF1y(x,y,z) 
m(1,3)=DF1z(x,y,z) 
m(2,1)=DF2x(x,y,z) 
m(2,2)=DF2y(x,y,z)
m(2,3)=DF2z(x,y,z) 
m(3,1)=DF3x(x,y,z) 
m(3,2)=DF3y(x,y,z) 
m(3,3)=DF3z(x,y,z) 
b(1)=-F1(x,y,z)
b(2)=-F2(x,y,z)
b(3)=-F3(x,y,z)
end subroutine generasistema






subroutine relajacion
implicit none
integer :: iters
real (kind=8) :: sum1,sum2
logical :: converge
real(kind=8),dimension(3) :: rx,rxold
rx(1)=x
rx(2)=z
rx(3)=y
rxold(1)=x
rxold(2)=z
rxold(3)=y

do iters=1,100

  do i=1,3
    sum1=0
    do j=1,(i-1)
      sum1 = sum1 + (m(i,j) * rx(j))
    end do

    sum2=0
    do j= (i),3
      sum2 = sum2 + (m(i,j) * rx(j))
    end do
    rx(i)= rx(i) + (1/m(i,i)) * ( b(i) - sum1 - sum2 )
  end do

  aux1=0
  do j=1,3
    aux1= aux1 + (rx(j)-rxold(j))**2
  end do
  aux2=0
  do j=1,3
    aux2= aux2 + rx(j)**2
  end do

  converge = (  sqrt(aux1)/sqrt(aux2)  < 1e-9 )
  if (converge) exit

  do j=1,3
    rxold(j)=rx(j);
  end do

end do

x=x+rx(1)
y=y+rx(2)
z=z+rx(3)

end subroutine relajacion




subroutine txtsalida
open(unit=1,file='salida.txt')
100 format(A30)
200 format (3(2x,f16.8))
select case (opcion)
  case (1)
    write(1,100) 'Iteracion funcional'
  case (2)
    write(1,100) 'Newton-Raphson'
end select
write(1,100) 'Solucion:'
write(1,200) x,y,z
close(1)
end subroutine txtsalida


real function F1(x,y,z)
real(kind=8) :: x,y,z
F1 = 3_8*x+cos(y*z)-0.5_8
return
end function F1

real function G1(y,z)
real(kind=8) :: y,z
G1 = (0.5_8-cos(y*z))/3_8
return
end function G1

real function F2(x,y,z)
real(kind=8) :: x,y,z
F2 = x**2 - 81_8*(y+0.1_8)**2 + sin(z) + 1.06_8
return
end function F2

real function G2(x,z)
real(kind=8) :: x,z
G2 =(sqrt(abs((x**2+sin(z)-1.06_8)/81_8))) - 0.1_8
return
end function G2

real function F3(x,y,z)
real(kind=8) :: x,y,z,pi
pi = dacos(-1.0_8)
F3 = exp(-x*y) + 20_8*z +(5*pi-1.5_8)
return
end function F3

real function G3(x,y)
real(kind=8) :: x,y
G3 = (exp(-x*y) +(5*3.14159_8-1.5_8))/(-20_8)
return
end function G3

real function DF1x(x,y,z)
real(kind=8) :: x,y,z
DF1x = 3_8
return
end function DF1x

real function DF1y(x,y,z)
real(kind=8) :: x,y,z
DF1y = -sin(y*z)*z
return
end function DF1y

real function DF1z(x,y,z)
real(kind=8) :: x,y,z
DF1z = -sin(y*z)*y
return
end function DF1z

real function DF2x(x,y,z)
real(kind=8) :: x,y,z
DF2x = 2_8*x
return
end function DF2x

real function DF2y(x,y,z)
real(kind=8) :: x,y,z
DF2y = -162_8*y-16.2_8
return
end function DF2y

real function DF2z(x,y,z)
real(kind=8) :: x,y,z
DF2z = cos(z)
return
end function DF2z

real function DF3x(x,y,z)
real(kind=8) :: x,y,z
DF3x = (exp(-x*y))*(-y)
return
end function DF3x

real function DF3y(x,y,z)
real(kind=8) :: x,y,z
DF3y = (exp(-x*y))*(-x)
return
end function DF3y

real function DF3z(x,y,z)
real(kind=8) :: x,y,z
DF3z = 20_8
return
end function DF3z

end program

