      SUBROUTINE OZONE(FI,P)
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NS=3, NS1=NS+2)
      DIMENSION FI(NS1,ND),P(ND)
C
C   THIS SUBROUTINE ESTIMATES O3 MIXING RATIO AS A FUNCTION
C   OF PRESSURE.
C
      DO 1 J=1,ND
      Z = - LOG10(P(J))
      FI(4,J) = 3.5E-8 + (7.E-8 - 3.5E-8) * Z/.5
      IF(Z.GT.0.5) FI(4,J) = 7.E-8 + (5.E-7 - 7.E-8)
     2  * (Z - 0.5)/.5
      IF(Z.GT.1.) FI(4,J) = 5.E-7 + (3.3E-6 - 5.E-7)
     2  * (Z - 1.)/.5
      IF(Z.GT.1.5) FI(4,J) = 3.3E-6 + (6.E-6 - 3.3E-6)
     2  * (Z - 1.5)/.5
      IF(Z.GT.2) FI(4,J) = 6.E-6 + (5.E-6 - 6.E-6) *
     2  (Z - 2.)/.5
      IF(Z.GT.2.5) FI(4,J) = 5.E-6 + (2.E-6 - 5.E-6) *
     2  (Z - 2.5)/.5
      IF(Z.GT.3) FI(4,J) = 2.E-6 + (9.E-7 - 2.E-6) *
     2  (Z - 3.)/.5
      IF(Z.GT.3.5) FI(4,J) = 9.E-7
   1  FI(4,J) = AMAX1(FI(4,J),1.E-9)
C
C   ESTIMATE COLUMN DEPTH ABOVE THE TOP OF THE MODEL
      P1 = P(1)
      UTOP = 0.
      IF(P1.GT.6.E-4)  UTOP = .0004
      IF(P1.GT.1.3E-3) UTOP = .0013
      IF(P1.GT.2.8E-3) UTOP = .004
      IF(P1.GT.5.E-3)  UTOP = .013
      IF(P1.GT.1.2E-2) UTOP = .04
      IF(P1.GT.2.3E-2) UTOP = .1
      IF(P1.GT.3.6E-2) UTOP = .13
      IF(P1.GT.5.6E-2) UTOP = .16
      IF(P1.GT.1.1E-1) UTOP = .21
      RETURN
      END
