program sel
implicit none
logical :: converged
integer :: maxit=100
integer :: i,j,k,opcion,n,iters
real (kind=8),dimension(:,:),allocatable :: m,l,u,iteraciones
real (kind=8),dimension(:),allocatable :: x,y,b,xold
real (kind=8) :: sum1,sum2,aux1,aux2,w

do
print*,'1.FACTORIZACION LU'
print*,'2.SOLUCION AL SEL USANDO RELAJACION'
print*,'3.MATRIZ TRIDIAGONAL'
print*,'4.SALIR'
read*,opcion
select case (opcion)
 
  case (1)
    call entrada
    call LU
    call txtsalida
    call cerrararray
  case(2)
    call entrada
    call relajacion
    call txtsalida  
    call cerrararray
  case(3)
    call generamatriz
    print*,'1.FACTORIZACION LU'
    print*,'2.SOLUCION AL SEL USANDO RELAJACION'
    read*,opcion
    if (opcion==1) then
      call tridiagonal
    else
      print*, 'Introduzca w:'
      read*,w
      call relajacion
    end if
    call txtsalida
    call cerrararray

  case(4)
     stop
   
  case default

end select

end do


contains

subroutine entrada
open(1,file='fichero.dat',form='formatted')
read(1,*) n
allocate(m(n,n))
allocate(l(n,n))
allocate(u(n,n))
allocate(iteraciones(maxit,n))
allocate(xold(n))
allocate(x(n))
allocate(y(n))
allocate(b(n))
do i=1,n
  read(1,*) (m(i,j),j=1,n)
end do 
do i=1,n
  read(1,*) b(i)
end do

if (opcion == 2) then
  read(1,*) w
  do i=1,n
    read(1,*) x(i)
  end do
end if

close(1)
end subroutine entrada

subroutine cerrararray 
deallocate(m)
deallocate(l)
deallocate(u)
deallocate(iteraciones)
deallocate(xold)
deallocate(x)
deallocate(y)
deallocate(b)
end subroutine cerrararray 
 
subroutine LU
call fact
call adelante
call atras
end subroutine LU

subroutine fact

do i=2,n
  do j=1,i-1
    u(i,j)=0
  end do
end do

do j=1,n
  do i=1,j
    if (i==j) then 
l(i,j)=1
    else 
l(i,j)=0;
end if
  end do
end do

do i=1,n

  do j=i,n
    aux1 = 0
    do k=1,i-1
      aux1=aux1 + l(i,k)*u(k,j)
    end do
    u(i,j)=  m(i,j)-aux1 
  end do

  do j=1,i
  if (i/=n) then
    aux2 = 0    
    do k=1,j-1
      aux2 = aux2 + l(i+1,k)*u(k,j)
    end do
    l(i+1,j)=( m(i+1,j)-aux2 )/u(j,j)
  end if
  end do

end do

print*,'Matriz U:'
do i=1,n
  print*,(u(i,j),j=1,n)
end do

print*,'Matriz L:'
do i=1,n
  print*,(l(i,j),j=1,n)
end do


end subroutine fact

subroutine adelante
  do i=1,n
    y(i)=b(i)
    do j=1,(i-1)
	y(i)=y(i)-l(i,j)*y(j)
    end do
  end do

end subroutine adelante

subroutine atras
  x(n)=y(n)/u(n,n)
  do i=n-1,1,-1
    x(i)=y(i)/u(i,i)
    do j=n,(i+1),-1
	x(i)=x(i)-( u(i,j)*x(j) )/u(i,i)
    end do
  end do
print*,'Solucion:'
print*,(x(i),i=1,n)


end subroutine atras

subroutine relajacion
implicit none

do iters=1,maxit

  do i=1,n
    iteraciones(iters,i)=x(i)
  end do
  do i=1,n
    sum1=0
    do j=1,(i-1)
      sum1 = sum1 + (m(i,j) * x(j))
    end do

    sum2=0
    do j=(i),n
      sum2 = sum2 + (m(i,j) * x(j))
    end do

    x(i)= x(i) + (w/m(i,i)) * ( b(i) - sum1 - sum2 )
  end do

  aux1=0
  do j=1,n
    aux1= aux1 + x(j)**2
  end do

  aux2=0
  do j=1,n
    aux2= aux2 + xold(j)**2
  end do

  converged = ( abs( (sqrt(aux1)-sqrt(aux2)) / sqrt(aux1) ) <= 1e-9 )
  if (converged) exit

  do j=1,n
    xold(j)=x(j);
  end do


end do

end subroutine relajacion

subroutine tridiagonal
implicit none
real (kind=8) , dimension(n) :: a,aold,b,c,b1,a1

do i=1,n
  a(i)=m(i,i)
  print*,a(i)
end do

do i=1,n-1
  c(i)=m(i,i+1)
  print*,c(i)
end do

do i=2,n
  b(i)=m(i,i-1)
  print*,b(i)
end do

a1(1)=a(1)
  do j=2,n
    b1(j) = b(j)/a(j-1)
    a1(j) = a(j) - b1(j) * c(j-1)
  end do

do i=1,n
  do j=1,n
    if (i==j) then
      l(i,j)=1
    else
      l(i,j)=0
    end if
  end do
end do

do i=1,n
  do j=1,n
    u(i,j)=0
  end do
end do

do i=2,n
  l(i,i-1)=b1(i)
end do

do i=1,n-1
  u(i,i)=a1(i)
  u(i,i+1)=c(i)
end do
u(n,n)=a1(n)

call adelante
call atras
end subroutine tridiagonal

subroutine generamatriz
real (kind=8) :: raiz
print*,'Introduzca la dimension de la matriz'
read*,n
allocate(m(n,n))
allocate(l(n,n))
allocate(u(n,n))
allocate(iteraciones(maxit,n))
allocate(xold(n))
allocate(x(n))
allocate(y(n))
allocate(b(n))
m(1,1)=1
m(1,2)=1
do j=3,n
  m(1,j)=0
end do
do j=1,n-2
  m(n,j)=0
end do
m(n,n)=1
m(n,n-1)=1

do i=2,n-1
  raiz=i
  do j=1,i-2
    m(i,j)=0
  end do
  m(i,i-1)=sqrt(raiz)
  m(i,i)=(2*i+1)
  m(i,i+1)=sqrt(raiz)
  do j=i+2,n
    m(i,j)=0
  end do
end do

b(1)=1
do i=2,n-1
  raiz=i
  b(i)=sqrt(raiz)
end do
b(n)=0

do i=1,n
  x(i)=0
end do

end subroutine generamatriz






subroutine txtsalida
implicit none
open(1,file='salida.txt')
100 format(A30)
200 format(E15.7)
write(1,100) 'Matriz de entrada:'
do i=1,n
  write(1,*) (m(i,j),j=1,n)
end do

write(1,100) 'Termino independiente:'
write(1,200) (b(i),i=1,n)

if (opcion == 2) then
  write(1,100) 'Parametro w :'
  write(1,200) w
  write(1,100) 'Aproximaciones:'
  write(1,100) 'Iteracion -----> Aproximacion'
  do i=1,(iters-1)
    write(1,*) i, (iteraciones(i,k),k=1,n,1)
  end do
end if
if (opcion == 1) then
    write(1,100) 'Solucion :'
    write(1,*) (x(i),i=1,n)
end if

end subroutine txtsalida


end program

