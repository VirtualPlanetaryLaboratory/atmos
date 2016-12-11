
      SUBROUTINE TRIDAG(A,B,C,R,U,N)
C      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NMAX=120)
      DIMENSION GAM(NMAX),A(N),B(N),C(N),R(N),U(N)
c      IF(B(1).EQ.0.)PAUSE
      BET=B(1)
      U(1)=R(1)/BET
      DO 11 J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
c        IF(BET.EQ.0.)PAUSE
        U(J)=(R(J)-A(J)*U(J-1))/BET
11    CONTINUE
      DO 12 J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
12    CONTINUE
      RETURN
      END
