      FUNCTION RELHUM(P)
         
C
C   THIS FUNCTION CALCULATES THE RELATIVE HUMIDITY (RELHUM) AT A
C   GIVEN PRESSURE P.  IT IS CURRENTLY SET UP TO YIELD EITHER A
C   STANDARD MANABE/WETHERALD RH PROFILE OR A FULLY-SATURATED
C   ATMOSPHERE, DEPENDING UPON THE VALUE OF IMW.
C  
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY

      Q = P/PG
      Q2 = AMAX1(Q-0.02,1.E-10)
      RELHUM = RSURF * Q2/0.98
      IF (IMW.EQ.0) RELHUM = 1.
      IF (IMW.EQ.3) RELHUM = 0.5
      IF (IMW.EQ.4) RELHUM = 0.1
C-KK  added for low-O2 environments, to prevent water from
C-KK  zeroing itself out. 8% rel hum is present-day atmosphere
C-KK  at approximately 15 km (cold trap level)
      IF (RELHUM .LT. 0.08) RELHUM = 0.08
      IF (IMW.EQ.5) RELHUM = 1.e-10
      END
