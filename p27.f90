PROGRAM PRACTICA23
IMPLICIT NONE
REAL(kind=8)::X,X1,X2,TOL
REAL(kind=8),DIMENSION(:),ALLOCATABLE::Y
INTEGER::MAXIT,ITERS,OPCION,I,M
LOGICAL::CONVERGED
CHARACTER(LEN=40):: TITULO

PRINT*,''
PRINT*,'*************** MENU ***************'
PRINT*,'ELIGE UNA DE LAS SIGUIENTES OPCIONES:'
PRINT*,''
PRINT*,'1.- SCHRODER'
PRINT*,'2.- NEWTON - RAPHSON'
PRINT*,'3.- CAUCHY'
PRINT*,'4.- MULLER'
PRINT*,'5.- SALIR'
PRINT*,'*************************************'

READ*,OPCION
SELECT CASE (OPCION)
CASE(1)
   PRINT*,'INTRODUCE EL PUNTO INICIAL'
   READ*,X 
   PRINT*,'INTRODUCE EL VALOR DE m'
   READ*,M
   PRINT*,'INTRODUCE NUMERO MAXIMO DE ITERACIONES'
   READ*,MAXIT
   ALLOCATE(Y(MAXIT))
   PRINT*,'INTRODUCE LA COTA DE ERROR'
   READ*,TOL
   CALL SCHRODER(X,M,MAXIT,TOL)   
   
CASE(2)
   PRINT*,'INTRODUCE EL PUNTO INICIAL'
   READ*,X
   PRINT*,'INTRODUCE NUMERO MAXIMO DE ITERACIONES'
   READ*,MAXIT
   ALLOCATE(Y(MAXIT))
   PRINT*,'INTRODUCE LA COTA DE ERROR'
   READ*,TOL
   CALL NEWTON(X,MAXIT,TOL)
CASE(3)
   PRINT*,'INTRODUCE EL PUNTO INICIAL'
   READ*,X
   PRINT*,'INTRODUCE NUMERO MAXIMO DE ITERACIONES'
   READ*,MAXIT
   ALLOCATE(Y(MAXIT))
   PRINT*,'INTRODUCE LA COTA DE ERROR'
   READ*,TOL
   CALL CAUCHY(X,MAXIT,TOL)
CASE(4)
   PRINT*,'INTRODUCE X0'
   READ*,X
   PRINT*,'INTRODUCE X1'
   READ*,X1
   PRINT*,'INTRODUCE X2'
   READ*,X2	
   PRINT*,'INTRODUCE NUMERO MAXIMO DE ITERACIONES'
   READ*,MAXIT
   ALLOCATE(Y(MAXIT))
   PRINT*,'INTRODUCE LA COTA DE ERROR'
   READ*,TOL
   CALL MULLER(X,X1,X2,MAXIT,TOL)
CASE DEFAULT 
   STOP
END SELECT

TITULO='ITERACION...APROXIMACION'
OPEN(UNIT=1,FILE='salida.txt')
REWIND(1)
IF (OPCION==4) THEN
WRITE(1,*) 'PUNTOS INICIALES= ',X,X1,X2
ELSE
WRITE(1,*) 'PUNTO INICIAL= ',X
END IF
IF (OPCION==1) WRITE(1,*) 'VALOR DE m= ',M
WRITE(1,*) 'NUMERO DE ITERACIONES= ',ITERS
WRITE(1,*) 'COTA DE ERROR= ',TOL
WRITE(1,*) ''
WRITE(1,*) TITULO
WRITE(1,100) 
100 FORMAT(A80)
DO I=1,ITERS
WRITE(1,'(I6,1x,F20.14)') I,Y(I)
END DO
CLOSE(1)
DEALLOCATE(Y)

CONTAINS

SUBROUTINE SCHRODER(X,M,MAXIT,TOL)
REAL(kind=8)::X,TOL,A
INTEGER::MAXIT,M

A=X
DO ITERS=1,MAXIT
 IF ((ABS(PP(A))<1E-10) .AND. (ABS(P(A))<1E-10)) EXIT
 Y(ITERS)=A-(M*(P(A)/PP(A)))
 CONVERGED=(ABS(Y(ITERS)-A)/ABS(Y(ITERS))<TOL)
 A=Y(ITERS)
 IF (CONVERGED) EXIT
END DO

PRINT*,'APROXIMACION=',A
PRINT*,'ITERACIONES=',ITERS

END SUBROUTINE SCHRODER


SUBROUTINE NEWTON(X,MAXIT,TOL)
REAL(kind=8)::X,TOL,A
INTEGER::MAXIT

A=X
DO ITERS=1,MAXIT
  IF ((ABS(FP(A))<1E-10) .AND. (ABS(F(A))<1E-10)) EXIT
  IF ((ABS(PP(A))<1E-10) .AND. (ABS(P(A))<1E-10)) EXIT
  IF ((ABS(PP(A))<1E-10) .AND. (ABS(P(A)*PP2(A))<1E-10)) EXIT
  Y(ITERS)=A-(F(A)/FP(A))
  CONVERGED=(ABS(Y(ITERS)-A)/ABS(Y(ITERS))<TOL)
  A=Y(ITERS)
  IF (CONVERGED) EXIT
END DO

PRINT*,'APROXIMACION=',A
PRINT*,'ITERACIONES=',ITERS

END SUBROUTINE NEWTON


SUBROUTINE CAUCHY(X,MAXIT,TOL)
REAL(kind=8)::X,TOL,A,R
INTEGER::MAXIT

A=X
DO ITERS=1,MAXIT
  R=SQRT(ABS((FP(A)**2)-(2*F(A)*FP2(A))))
  IF (((ABS(FP(A)+R)<1E-10) .OR. (ABS(FP(A)-R)<1E-10)) .AND. (ABS(F(A))<1E-10)) EXIT
  IF ((ABS(PP(A))<1E-10) .AND. (ABS(P(A))<1E-10)) EXIT
  IF ((ABS(PP(A))<1E-10) .AND. (ABS(P(A)*PP2(A))<1E-10)) EXIT
  IF ((ABS(PP(A))<1E-10) .AND. (ABS(P(X)*PP3(X))<1E-10)) EXIT
  IF (FP(A)>=0) Y(ITERS)=A-(2*F(A)/(FP(A)+R))
  IF (FP(A)<0) Y(ITERS)=A-(2*P(A)/(PP(A)-R))
  CONVERGED=(ABS(Y(ITERS)-A)/ABS(Y(ITERS))<TOL)
  A=Y(ITERS)
  IF (CONVERGED) EXIT
END DO

PRINT*,'APROXIMACION=',A
PRINT*,'ITERACIONES=',ITERS

END SUBROUTINE CAUCHY

SUBROUTINE MULLER(X,X1,X2,MAXIT,TOL)
REAL(kind=8)::X,X1,X2,TOL,E0,E1,H0,H1,A,B,C,MAXIMO,P1,P2,P3
REAL(kind=8),DIMENSION(3)::Z
INTEGER::MAXIT

P1=X
P2=X1
P3=X2
DO ITERS=1,MAXIT
    E0=F(P1)-F(P3)
    E1=F(P2)-F(P3)
    H0=P1-P3
    H1=P2-P3
    A=((E0*H1-E1*H0)/(H1*H0*H0-H0*H1*H1))
    B=((E1*H0*H0-E0*H1*H1)/(H1*H0*H0-H0*H1*H1))
    C=F(P3)
    IF (C<1E-10) EXIT
    Y(ITERS)=P3+((-2_8*C)/(B+(B/ABS(B))*SQRT(ABS(B*B-4_8*A*C))))
    CONVERGED=(ABS(Y(ITERS)-P3)/ABS(Y(ITERS))<TOL)
    Z(1)=ABS(Y(ITERS)-P1)    
    Z(2)=ABS(Y(ITERS)-P2)    
    Z(3)=ABS(Y(ITERS)-P3)    
    MAXIMO=MAXVAL(Z)
    IF (MAXIMO==Z(1)) THEN
		P2=Z(2)
		P1=Z(3)
	ELSE IF (MAXIMO==Z(2)) THEN
			P2=Z(1)
			P1=Z(3)
		ELSE IF (MAXIMO==Z(3)) THEN
				P2=Z(1)
				P1=Z(2)
    END IF
    P3=Y(ITERS) 
    IF (CONVERGED) EXIT
END DO

PRINT*,'APROXIMACION= ',P3
PRINT*,'ITERACIONES=',ITERS

END SUBROUTINE MULLER

REAL(kind=8) FUNCTION P(X)
REAL(kind=8)::X
P=X**5-6.15_8*(X**4)+13.76_8*(X**3)-12.544_8*(X**2)+2.4576_8*X+1.7384_8
RETURN
END FUNCTION P

REAL(kind=8) FUNCTION PP(X)
REAL(kind=8)::X
PP=5_8*(X**4)-24.6_8*(X**3)+41.28_8*(X**2)-25.088_8*X+2.4576_8
RETURN
END FUNCTION PP

REAL(kind=8) FUNCTION PP2(X)
REAL(kind=8)::X
PP2=20_8*(X**3)-73.8_8*(X**2)+82.56_8*X-25.088_8
RETURN
END FUNCTION PP2

REAL(kind=8) FUNCTION PP3(X)
REAL(kind=8)::X
PP3=60_8*(X**2)-147.6_8*X+82.56-8
RETURN
END FUNCTION PP3

REAL(kind=8) FUNCTION F(X)
REAL(kind=8)::X
F=P(X)/PP(X)
RETURN
END FUNCTION F

REAL(kind=8) FUNCTION FP(X)
REAL(kind=8)::X
FP=((PP(X)*PP(X))-(P(X)*PP2(X)))/(PP(X)*PP(X))
RETURN
END FUNCTION FP

REAL(kind=8) FUNCTION FP2(X)
REAL(kind=8)::X
FP2=((((2*PP(X)*PP2(X))-(PP(X)*PP2(X)+P(X)*PP3(X)))*(PP(X)*PP(X)))-(((PP(X)*PP(X))-(P(X)*PP2(X)))*(2*PP(X)*PP2(X))))/((PP(X))**4)
RETURN
END FUNCTION FP2

END PROGRAM


