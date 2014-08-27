      SUBROUTINE AQUEOUS(X,F,I)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      DIMENSION X(NAQ),F(NAQ)
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/WBLOK.inc'
      DATA LSO2g,LH2COg,LSO2aq,LH2COaq,LHCO3_,LCO3_2,LHSO3_,
     1  LSO3_2,LH2COSO3,LOH_/1,2,3,4,5,6,7,8,9,10/
C
C     THIS SUBROUTINE DOES THE AQUEOUS PHASE CHEMISTRY DESCRIBED IN
C     SUBROUTINE RAINOUT
C
C     1)  (SO2)g + ALPHARAIN*[(SO2)aq  +  HSO3-  +  SO3=  +  CH2OHSO3-]
C                  =  (SO2)go
C     2)  (SO2)g   =  (SO2)aq
C     3)  (H2CO)g  =  CH2(OH)2
C     4)  (CO2)aq  =  HCO3-  +  H+
C     5)  (SO2)aq  =  HSO3-  +  H+
C     6)   HCO3-   =  CO3=   +  H+
C     7)   HSO3-   =  SO3=   +  H+
C     8)  CH2(OH)2  +  HSO3-  =  H2O  +  CH2OHSO3-
C     9)   H2O     =  H+  +  OH-
C    10)  (H2CO)g + ALPHARAIN*[CH2(OH)2  +  CH2OHSO3-]  =  (H2CO)go

      HPLUS = X(LHCO3_) + X(LHSO3_) + X(LH2COSO3) 
     1  + X(LOH_) + 2.*(X(LCO3_2) + X(LSO3_2) + SO4_2)
      F(1) = X(LSO2g) - SO2g0 + ALPHARAIN*( X(LSO2aq) 
     1  + X(LHSO3_) + X(LSO3_2) + X(LH2COSO3) )
      F(2) = X(LSO2aq) - HSO2*X(LSO2g)
      F(3) = X(LH2COaq) - HH2CO*X(LH2COg)
      F(4) = X(LHCO3_)*HPLUS - R4(I)*CO2aq
      F(5) = X(LHSO3_)*HPLUS - R5(I)*X(LSO2aq)
      F(6) = X(LCO3_2)*HPLUS - R6(I)*X(LHCO3_)
      F(7) = X(LSO3_2)*HPLUS - R7(I)*X(LHSO3_)
      F(8) = X(LH2COSO3) - R8(I)*X(LH2COaq)*X(LHSO3_)
      F(9) = X(LOH_)*HPLUS - R9(I)
      F(10) = X(LH2COg) - H2COg0 + ALPHARAIN*(X(LH2COaq)+X(LH2COSO3))
C
      RETURN
      END
