      SUBROUTINE LTNING
c  3-21-06   This subroutine does not work
c  it does not conserve redc  it needs to be replaced by something that isn't wrong
c-mc 4-29-06  it conserves redox now and is only as wrong as chamedies, which is probably wrong...


      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      character*8 PLANET
      real*8 mass
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/LTBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'

C
C     THIS SUBROUTINE CALCULATES LIGHTNING PRODUCTION RATES FOR O2
C     AND NO IN AN N2-O2-CO2-H2O ATMOSPHERE USING THE EQUATIONS IN
C     THESIS APPENDIX C.
c     this is incomprehensible without more hints
C
C     EQUILIBRIUM CONSTANTS AT 3500 K
      AK1 = .103                ! = pNo/sqrt(pN2*pO2)  dimensionless
      AK2 = .619                ! = pCO2/PCO/sqrt(pO2) bars
      AK3 = 5.3                 ! = pH20/pH2/sqrt(pO2) bars
      AK4 = 0.22                ! = pO/sqrt(pO2)

   
c      Pbar = P0/1.013E6   !converting from pascals to bars
      Pbar = P0/1.013   !converting from pascals to bars
      
      FCO = USOL(LCO,1) 
      PH2O = USOL(LH2O,1) * Pbar
      PH2 = USOL(LH2,1) * Pbar   ! I am having problems here if PH2 gets too big
      PCH4 = USOL(LCH4,1) * Pbar

     
     
      FH2O = PH2O /PBAR
      FN2 = 1. - FO2 - FCO2 - FCO - FAR - FH2O - FCH4   ! EWS - now includes Argon

      PN2 = FN2*Pbar
      PO2 = FO2*Pbar
      PCO2 = FCO2*Pbar
      PCO = FCO * Pbar
      
      O2T = PO2 + PCO2 + 0.5*(PH2O + PCO)
      H2T = PH2 + PH2O
      CT = PCO2 + PCO

    

      ALPHA = AK2*SQRT(O2T)
      BETA = AK3*SQRT(O2T)
      A = (AK1*SQRT(PN2) + AK4)/(2.*SQRT(O2T))    !K_T in jim's thesis
      B = 0.5*CT/O2T
      C = 0.5*H2T/O2T

C     INITIAL GUESS FOR XO2 AT 3500 K
      X = 0.1 + 0.9*PO2/O2T + 0.2*PCO2/O2T ! this initial guess seems to be causing trouble sometimes 
      X_orig = X
 
C
C     NEWTON STEP
      DO 1 N=1,20     
 400  NS = N
      XS = X                        
      X2 = SQRT(X)            
      FX = X + A*X2 - B/(1.+ALPHA*X2) - C/(1.+BETA*X2) + 2.*B + C - 1.            
       ! ^^^ C18) in JFK thesis
      FPX = 1. + (A + ALPHA*B/(1.+ALPHA*X2)**2 + BETA*C/(1.+BETA*X2)
     2  **2)/(2.*X2) 
      X = X - FX/FPX  !GNA - bad initial guess for X (if X too big) is making this occasionally go negative for me, which make PO2 negative, which wrecks havoc in the sqrt(PN2*PO2) below...
      IF(X.LT.0) THEN  
         X = X_orig * 0.95 !decrease initial guess slightly until small enough
         X_orig = X
         print *, 'fixing x in ltning.f. new x = ', x
         print *, 'old x = ', x/0.95
         GOTO 400 
      ENDIF
   
      ERR = ABS((X-XS)/X)
     
      IF(ERR.LT.1.E-5) GO TO 2
   1  CONTINUE
   2  PO2 = X*O2T
      PNO = AK1*SQRT(PN2*PO2)
      PCO = CT/(1. + AK2*SQRT(PO2))   !added 3-20-06
      PH2 = H2T/(1. + AK3*SQRT(PO2))   !added 3-20-06

c-mc 
      PO = AK4*sqrt(PO2)   !added 4-28-06
c-mc

c      print *, PO2, 0.5*(PCO+PH2-PO-PNO)

c      stop
C
C     SCALE AGAINST ESTIMATED COLUMN PRODUCTION OF NO IN THE PRESENT-
C     DAY TROPOSPHERE.  DISTRIBUTE PRODUCTION OVER LOWEST 6 KM.
C
C   COLUMN-INTEGRATED NO PRODUCTION RATE IN PRESENT ATMOSPHERE IS
C   EQUAL TO PRONO
      PRNOX = PRONO/7.942E5
      PNONOW = 3.574E-2
      PRONOP = PRONO * PNO/PNONOW

      ZAPNO = PRNOX*PNO/PNONOW
      ZAPO2 = ZAPNO*PO2/PNO

      ZAPCO = ZAPNO*PCO/PNO   ! added 3-20-06
      ZAPH2 = ZAPNO*PH2/PNO   ! added 3-20-06

      ZAPO = ZAPNO*PO/PNO     ! added 4-28-06
      ZAPCO2 = ZAPNO*PCO2/PNO
      ZAPH2O = ZAPNO*PH2O/PNO

     

c      print *, ZAPNO, ZAPO2,ZAPCO,ZAPH2,ZAPO,ZAPCO,ZAPH2O
c      stop

c  - add below
      write(14, 100) PO2,PNO,ZAPO2,ZAPNO,PRONO,PRONOP,ZAPCO,ZAPH2
 100  FORMAT(/1X,'PO2=',1PE10.3,2X,'PNO=',1PE10.3,2X,'ZAPO2=',
     2  E10.3,2X,'ZAPNO =',E10.3,2X,'PRONO =',E10.3,2X,'PRONOP =',
     3  E10.3,2X,'ZAPCO =',E10.3,2X,'ZAPH2 =',E10.3)
      write(14, 101) NS,X,ERR
 101  FORMAT(1X,'N=',I2,5X,'X=',1PE12.5,2X,'ERR=',1PE12.5)
      RETURN
      END
