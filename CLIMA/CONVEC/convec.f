      SUBROUTINE CONVEC(T1,T2,P1,P2,FH1,FH2,FC1,FC2,DZ,ITROP,
     2          cflag,Idry,IMCO2)
C   THIS SUBROUTINE FINDS THE CONVECTIVE LAPSE RATE BETWEEN GRID
C   POINTS J1 AND J. IT CONSIDERS THE CASE FOR HIGH CO2.

      INCLUDE 'CLIMA/INCLUDE/header.inc'
  
      PARAMETER(NT=76, MT=36)
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5) ! Adding parameter statement needed for FI(NS1,ND) 5/23/2011
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND),
     &  FH22
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     &  BETA1,BETA2,FVDRY,PDRY
      COMMON/FBLOK/TTAB(NT),PVAP(NT),DPVAP(NT),SVC(NT),DSV(NT),DSC(NT)
     &  ,RHOV(NT),DRHOV(NT),BETAM(70,75),TCP(75),PCP(70),DPDTL(70,75),
     &  DRDTL(70,75)
      COMMON/GBLOK/TCTAB(MT),PCVAP(MT),BETASC(MT),DPCVAP(MT),
     2  DRCVAP(MT),SVSC(MT),DSCC(MT),TKTAB(MT),TCC(25),PCC(36),
     3  BETAMC(25,36),CPC(25,36),DVDTC(25,36),DVDPC(25,36),
     4  DSVC(MT)
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,AMN, AMN2  ! AMN and AMN2 are DM and DM2 respectively (from Clima.f) 5/23/2011
      COMMON/SBLOK/P0P,T0P,R,SUBL
      COMMON/CO2BLOK/betac1,betac2,PC0,TC0,VAPCL0,SUBCL0,DLVCDT
     & ,DLSCDT,CCL,CCS 
      
C
C   THIS SUBROUTINE FINDS THE CONVECTIVE LAPSE RATE BETWEEN GRID
C   POINTS J1 AND J.
C
C JK  Idry is a parameter that has been added in order to compute the dry adiabatic lapse rate
c    
c-as This subroutine has been modified to take into account CO2 condensation (sept-2004)
c  Tag code: 1-Main loop
c            Units-Options IN the main loop   
c            1x-Moist adiabatic lapse rate
c            2x-CO2 condensation   
      itrop = itrop !EWS - no complaining
      NDIV = 10                     !number of sublevels
      PLT = DZ
      IF(P2.LT.P1) PLT = - DZ
      DLP = PLT/FLOAT(NDIV)
      AC = 0.
      CALL SATRAT(T1,PSAT)  
      IF (IMOIST.EQ.1) FH1 = PSAT/P1
c      print *, 'FH1=', FH1
      FH1 = amin1(fh1,0.99999999) ! c-rr changed from .99 to .9999 to do the runaway greenhouse 5/26/2011
      AV = 18./AMN * FH1/(1.-FH1)   !THIS VALUE OF AV IS ONLY USED BELOW 273 K
      T = T1
      P = P1
      TL = LOG(T)
      PL1 = LOG(P1)

      cflag=0.      !convection flag
c  Flag to indicate that CO2 has reach the saturation pressure
c (0= not saturated, 1=saturated)
c      imco2 = 0     
!       print 1499
! 1499  format(i2)
!       print 1500, P,T,cflag,idry,imoist,PSAT
! 1500  format(1x,'P =',1pe10.3,2x,'T =',e10.3,2x,
!     .           'CFLAG =',e10.3,2x'IDRY =',I3,2x,
!     .           'IMOIST =',I3,2x,'PSAT =',e10.3)
C-KK        Begin Sublevel Integration
      DO 1 L=1,NDIV
c Calculation of heat capacities (cal/mol/K)
c-rr          Putting new curve fit equations 3/29/11
       CPO2 = 7.47 -4.84E-3*T + 1.38E-5*T*T 
     &         -8.73E-9*T*T*T - 1.76E-9/T/T        
        CPCO2 = 5.89 + 6.06E-3*T + 2.39E-5*T*T 
     &          -3.44E-8*T*T*T  
       CPN2 = 6.76 + 6.06E-4*T + 1.3E-7*T*T

c-rr   old curve fit equations 3/29/11
c      CPCO2 = 7.7 + 5.3E-3*T - 8.3E-7*T*T
c       CPN2 = 6.76 + 6.06E-4*T + 1.3E-7*T*T

c      CPO2 = 8.27 + 2.58E-4*T - 1.877E5/T/T

      
      CPO2 = AMAX1(CPO2,CPN2)
      CPCH4 = 8.3
  !=====================================================================

c        print *, 'FCH4NC=', FCH4NC
c         if (L.eq.1)then
c          print *, 'FH1=', FH1
c         print *, 'Total F=',FI(2,J)+FN2+FO2NC+FARNC+FCH4NC
c         endif
!=====================================================================

!      SNC1  = FN2 + FO2+ FAR + FCH4 
!      print *, 'SNC1=', SNC1
!      FN2P  = FN2/SNC1
!      FO2P  = FO2/SNC1
!      FARP  = FAR/SNC1
!      FCH4P = FCH4/SNC1
!      FC1P = FC1/SNC1
       
!      CPN2 = FN2P*CPN2 + FO2P*CPO2 + FARP*4.97 + FCH4P*CPCH4 ! New normalized CPN for Co2 moist adiabat

      CPN = FN2*CPN2 + FO2*CPO2 + FAR*4.97 + FCH4*CPCH4   ! This CPN is the total CPN for the noncondensibles only. This CPN is only used for the moist CO2 adiabat (so H2O is not included because it is a condensible. In the moist adiabat, CO2 is a noncondensible).
    
!      SUMF = FN2 + FO2 + FAR + FCH4

!============================================================================
! CPNW needs to be defined before any of the adiabats are calculated

      CPNW   = FC1*CPCO2 + (1.-FC1)*CPN   ! Added CPNW
!============================================================================

      CALL SATCO2(T,PSCO2)
      FCSAT = PSCO2/P     
c Once CO2 reaches the saturation pressure it stays saturated up to JTROP
      if(FCSAT.lt.FC1) imco2=1   
c      if(FCSAT.lt.FC1/1.34) imco2=1  ! supersaturating to 1.34?
      PL = PL1 + (L-1)*DLP
      P = EXP(PL)
      TC = T - 273.15
      If (Idry.eq.1) go to 40
      if(imco2.eq.0) go to 3

c   Moist CO2 adiabat
      cflag = 3. 
      N = (TC + 130.)/5. + 1
      N = MAX0(N,1)
      IF (TC .GT. -56.595) N = (TC + 130.)/5. + 3

c-mm  Next line added to allow for surface temps above 303K, since this is
c-mm  the limit of the table.
      IF (TC.GT.29.) N = MT - 1
      N1 = N + 1
      FRC = (TC - TCTAB(N))/(TCTAB(N1) - TCTAB(N))
      BET = FRC*BETASC(N1) + (1.-FRC)*BETASC(N)
c-mm  Treat CO2 as ideal if above 303K
      IF (TC.GT.29.) BET = 1.
      AVC = 44./AMN2 * BET * FCSAT/(1. - FCSAT)  ! All AMNs in CO2 condensation section changed to AMN2. 5/23/2011
      IF (L.EQ.1) BETAC1 = BET
      TSQ = T*T
      IF (T .LT. 216.56) GO TO 21


C   Lapse rate in pure CO2 over liquid CO2
c      if (FCO2.gt.0.9) then
      DLPVLT = 2.303*T*(867.2124/TSQ + 18.65612E-3 - 2.*72.4882E-6*T
     &  + 3.*93.E-9*TSQ)
      DLADLT = 0.
c      endif
      GO TO 22


C   Lapse rate in pure CO2 over solid CO2
c21    if(FCO2.gt.0.9) then
   21     CONTINUE

       T47 = T - 4.718
       DLPVLT = 2.303*T*(1284.07/T47/T47 + 1.256E-4)
       DLADLT=0.
c      endif
C   Note that, entropy for CO2 is expressed in units of cal/g-K,
C   so that the conversion factor 4.184 disappears from
C   the expression for DLADLT
   22 DSVDLT = FRC*DSVC(N1) + (1.-FRC)*DSVC(N)
      SVCM = FRC*SVSC(N1) + (1.-FRC)*SVSC(N)
      BET = FRC*BETASC(N1) + (1.-FRC)*BETASC(N)
      DLRVLT = FRC*DRCVAP(N1) + (1.-FRC)*DRCVAP(N)
      CVN = CPN - R
      SUM1 = (R*DLRVLT - CVN)/AMN2 - AVC*DSVDLT
      SUM2 = AVC*SVCM + R/AMN2
      DLADLT = SUM1/SUM2
      DLPDLT = DLPVLT- DLADLT/(1. + AVC*AMN2/44./BET)
C
      GO TO 2
c   WATER MOIST ADIABAT




   3   CALL SATRAT(T,PSAT)
       cflag =1.
       
C-KK        Calculate a dry standard adiabat above tropopause. 
c--JFK 7/14/08 Dont allow these dry adiabats. They seem to be
c      occurring within the mesosphere, which is unphysical.
c      IF (ITROP .EQ. 0) THEN 
       IF (IMW.eq.5) THEN   ! This does the dry adiabat only
        DLPDLT = CPN/R
c        print *, 'CPN=',CPN, 'R=', R
        cflag = 2.
        GO TO 2
      END IF

      IF (IMOIST.EQ.1) GO TO 11

      IF (L.EQ.1) BET = BETA2
      FP = 1./(1. + BET*(1.-FVDRY)/FVDRY)
      PH2O = FP * P
      IF (PH2O.LT.PSAT) GO TO 12


  11  IMOIST = 1
      IF(T.LT.273.16) GO TO 13

C   INGERSOLLS FORMULATION BETWEEN 0 AND 374 C
      N = TC/5. + 1
      N1 = N + 1
      FR = (TC - TTAB(N))/(TTAB(N1) - TTAB(N))
      DLPVLT = FR*DPVAP(N1) + (1.-FR)*DPVAP(N)
      DSVDLT = FR*DSV(N1) + (1.-FR)*DSV(N)
      SVCM = FR*SVC(N1) + (1.-FR)*SVC(N)
      PV = FR*PVAP(N1) + (1.-FR)*PVAP(N)
      RV = FR*RHOV(N1) + (1.-FR)*RHOV(N)
      DLRVLT = FR*DRHOV(N1) + (1.-FR)*DRHOV(N)
      BET = 41.84*RV*R*T/(18.*PV)


!    New statements to normalize mixing ratios in water moist adiabat 8/31/2011.. Adds FCO2


!      SNCW  = 1 - FH1
!      FN2PW  = FN2/SNCW
!      FO2PW  = FO2/SNCW
!      FARPW  = FAR/SNCW
!      FCH4PW = FCH4/SNCW
!      FC1PW = FC1/SNCW

!      CPNW2 = FC1PW*CPCO2 + (1.- FC1PW)*CPN ! Added CPNW

! --------------------------------------------------------------------

      CVN = CPNW - R
      IF (L.NE.1) GO TO 14
      AV = 18./AMN * BET * FH1/(1.-FH1)
      BETA1 = BET
  14  CONTINUE

      SUM1 = 4.184*(R*DLRVLT - CVN)/AMN - AV*DSVDLT
      SUM2 = AV*SVCM + 4.184*R/AMN
      DLADLT = SUM1/SUM2
      DADLT = AV*DLADLT
      SUM3 = 1. + DLRVLT - DLADLT
      FAC = BET*18./AV/AMN
      DLPDLT = PV/P * (DLPVLT + FAC*SUM3)
      GO TO 2

C   DRY ADIABAT ABOVE 374 C
  12  DADLT = 0.
      DLADLT = 0.
      AV = 18./AMN * FVDRY/(1.-FVDRY)
      TC = T - 273.15
      IMOIST = 0
      PDRY = P
      N = TC/10. + 2
      IF (TC.GT.600.) N = (TC-600.)/100. + 62
      PV = FH1 * P
      M = PV/5. + 1
      M1 = M - 1
      IF (M.EQ.1) M1 = M

      FP = 0.
      IF (M.GT.1) FP = (PV - PCP(M))/(PCP(M) - PCP(M1))
   
      FT = (TC - TCP(N))/(TCP(N+1) - TCP(N))
      
      
      DLPVLTsave = DLPVLT

      DLPVLT = (1.+FP-FT)*DPDTL(M,N) - FP*DPDTL(M1,N) + FT*
     &  DPDTL(M,N+1)

 
!      DLPVLT = (1.+FP-FT)*DPDTL(M,N) - FP*(1.-FT)*DPDTL(M1,N) + FT*
!     &  (1.+FP)*DPDTL(M,N+1) - FP*FT*DPDTL(M1,N+1)




!===========================================================================
! The steam tables do not have DPDTL values below 1 bar. But below 1 bar
! water vapor behaves as an ideal gas at high surface temperatures (>=1400 K)
! so we can just use the idea gas law equation for DLPVLT

        if( P < 1.0) DLPVLT = CPNW/R
!===========================================================================




        BETAM(M,N) = AMIN1(4., BETAM(M,N))! Statements to correct BETAMN
!      print 16333, BETAM(M,N),T
!16333 format (1p2e12.5)
      BET = (1.+FP-FT)*BETAM(M,N) - FP*BETAM(M1,N) + FT*BETAM(M,N+1)


      IF (L.EQ.1) BETA1 = BET
      DLPDLT = FH1*DLPVLT + (1.-FH1)*CPNW/R
!      print 6969, TC,TCP(N),DLPDLT,CPNW,R,M,N
c 6969  format(1p5e14.5,0p,2x,i3,2x,i3) !EWS - not used
      GO TO 2

C   CLAUSIUS-CLAPEYRON EQUATION BELOW 0 C
  13  HL = SUBL
      HLT = HL/T
      DDLT = 0.
      CC = 0.5


!     New statements to normalize mixing ratios in water moist adiabat 8/31/2011


!      SNCW  = 1 - FH1
!      FN2PW  = FN2/SNCW
!      FO2PW  = FO2/SNCW
!      FARPW  = FAR/SNCW
!      FCH4PW = FCH4/SNCW
!      FC1PW = FC1/SNCW

!      print *, 'SUM=', FN2PW + FO2PW +FARPW + FCH4PW + FC1PW
      

      CPNW   = FC1*CPCO2 + (1.-FC1)*CPN  ! Correct CP noncondensible without CO2.
!      CPNW2 = FC1PW*CPCO2 + (1.- FC1PW)*CPN ! Added CPNW
!      print *, 'CPNW=', CPNW, 'CPNW2=', CPNW2
!-----------------------------------------------------

      SUM1 = (HLT*18. - CPNW)/AMN - AC*CC
      SUM2 = DDLT + CC - HLT
      SUM3 = HLT + R/AV/AMN
      DADLT = (SUM1 - AV*SUM2)/SUM3
      DLPDLT = 18.*HLT/R - DADLT/AV/(1. + AV*AMN/18.)
!      print 6969,DLPDLT,HLT,R,DADLT,AV,T
      DLPVLT = 18.*HLT/R
      BET = 1.
      GO TO 2


C Dry adiabat
  40  DLPDLT = CPNW*AMN/(4.184*R)


   2  DLT = DLP/DLPDLT     

      TL = TL + DLT
      T = EXP(TL)

      if(imco2.eq.0)then
      DAV = DADLT*DLT
      AV = AV + DAV
      endif

      if(imco2.eq.1)then
      DAV = DADLT*DLT
      AVC = AVC + DAV
      endif

   1  CONTINUE       !End Sublevel integration


      T2 = T
       
      
      FH2 = 1./(1. + BET*18./AMN/AV)

      

      CALL SATRAT(T2,PSAT)
      IF (PSAT.GT.P2) GOTO 30
      FH2 = RELHUM(P2) * PSAT/P2
!      print *,FH2,PSAT,P2,RELHUM(P2)
     
         
  30  IF(IMCO2.EQ.1) THEN
       BETA2 = 1.
       BETAC2 = BET
      ELSE
       BETA2=BET
      ENDIF  
   
      FC2 = FCO2
      CALL SATCO2(T2,PSCO2)
      FCSAT = PSCO2/P2
      IF (IMCO2.EQ.1) FC2 = FCSAT
!      imco2 = 0
      
!      print 1501,P2,T2
c 1501 format(1x,'P2 =',1pe10.3,2x,'T2 =',e10.3) !EWS - not used
      RETURN
      END
