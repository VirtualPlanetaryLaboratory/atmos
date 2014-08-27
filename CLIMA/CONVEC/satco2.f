      SUBROUTINE SATCO2(T,PSCO2)
c This subroutine calculates the vapor saturation pressure for CO2
      
      IF (T.LT.216.56) GO TO 1
C
C   VAPOR PRESSURE OVER LIQUID
      PSL = 3.128082 - 867.2124/T + 1.865612E-2*T
     2  - 7.248820E-5*T*T + 9.3E-8*T*T*T
      PATM = 10.**PSL
      PSCO2 = 1.013*PATM 
c-rr      PSCO2 = 1.013*PATM*1.34 ! Assuming saturation occurs when e/es is 1.34 not 1
      RETURN
C
C   VAPOR PRESSURE OVER SOLID
C   (Altered to match vapor pressure over liquid at triple point)
   1  PSL = 6.760956 - 1284.07/(T - 4.718) + 1.256E-4*(T - 143.15)
      PATM = 10.**PSL
      PSCO2 = 1.013*PATM 
c-rr       PSCO2 = 1.013*PATM*1.34 ! Assume saturation occurs when e/es is 1.34 not 1 
      RETURN
      END
