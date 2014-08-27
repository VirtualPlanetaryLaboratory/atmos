      SUBROUTINE CHEMPL(D,XP,XL,K)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      character*8 REACTYPE,CHEMJ
      DIMENSION XP(NZ),XL(NZ),D(NSP2,NZ)
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'

C
C   THIS SUBROUTINE CALCULATES CHEMICAL PRODUCTION AND LOSS RATES
C   USING THE INFORMATION IN THE MATRICES JCHEM, ILOSS, AND IPROD.
C   CALLED BY SUBROUTINE DOCHEM.
C
      DO 1 I=1,NZ
      XP(I) = 0.
   1  XL(I) = 0.
C
C   LOSS FREQUENCY XL
      NL = NUML(K)        !chempl is called with given species (K)- NUML is how many reactions involve K
      DO 2 L=1,NL
      J = ILOSS(1,K,L)    !reaction number for loss process
      M = ILOSS(2,K,L)    !reactant number for loss process
      DO 2 I=1,NZ
   2  XL(I) = XL(I) + A(J,I)*D(M,I) !rate*density of reactant 1  
!XL has units of s^-1
C
C   PRODUCTION RATE XP
      NPRX = NUMP(K)
      DO 3 L=1,NPRX
      J = IPROD(K,L)   !reaction number
      M = JCHEM(1,J)   !reactant 1 for reaction number J
      N = JCHEM(2,J)   !reactant 2 for reaction number J
      DO 3 I=1,NZ      !loop over height (also for each reaction)
   3  XP(I) = XP(I) + A(J,I)*D(M,I)*D(N,I) !rate*density1*density2 
!HV has a density of 1, so rates make up for this
!XP in units of mol/cm^3/s
      RETURN
      END
