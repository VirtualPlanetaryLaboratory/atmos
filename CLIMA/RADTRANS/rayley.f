       SUBROUTINE RAYLEY(SIGR,AL2,Fwater,FNCR)  ! 6/7/2012 added FNCR to RAYLEY argument
C
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NSOL=38, NGS=8) !gna ngas changed to 8
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5) ! Adding parameter statement needed for FI(NS1,ND) 5/23/2011          
c      REAL kmatrix_sol, weights ! EWS - not used
      COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1)

C
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND),
     &   FH22
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM, DM2  ! Passing DM2 which is AMN2 in CONVEC. 5/3/2011
      DIMENSION A(4),B(4),DEL(6),SIG(6)  ! Changed dimension statement so that A,B,and SIG have 4 entries. SIG,DEL has 5,6 for water methane and water, respectively. SIGR is just a scalar. 5/3/2011
C

      REAL SIGR, AIRSCAT, NM,RM, WN, KCF, SIGH
      DATA A/29.06,26.63,43.9,27.92/
      DATA B/7.7,5.07,6.4,5.6/
      DATA DEL/.0305,.054,.0805,0.032,0.0002,0.17/
   
C

!      print *, 'FNCR=', FNCR, Fwater, FCO2
!      pause

      !DO 1135 I=1,NSOL                 Sending SIGRs one wavelength at a time so do not need these statements 5/3/2011
      !   AL2=ALAMBDA(I)*ALAMBDA(I)
  
         AL4=AL2*AL2
         AL = AL4**.25      ! AL is wavelength (in microns)
         WN = 1E4*(1/AL)   ! WN is Wavenumber (1/cm)

cRKK- THE EXPRESSION FOR SIG(J) IS FROM VARDAVAS & CARVER, 1984, 
cRKK- PLANETARY SPACE SCIENCE, VOL 32, P-1307.
C

         DO 1140 J=1,4
            PAREN=(1.E-5*A(J)*(1.+1.E-3*B(J)/AL2))**2
            SMDEL=(6.+3.*DEL(J))/(6.-7.*DEL(J))
            SIG(J)=4.577E-21*SMDEL*PAREN/AL4
 1140    CONTINUE
C           c-rr 5/29/2011
          ! SIG(1) is Nitrogen. SIG(2) is oxygen . SIG(3) is carbon dioxide. SIG(4) IS ARGON.SIG(5) is methane. N2,O2,AR, and CH4 make up the non-condensible (NC)
          !rayleigh scattering components for air. SIG(6) is the rayleigh scattering component due to water. There is no known A term for water but it is assumed that A ~30. 
         ! Then, and since the term (A(J)*1E-5)^2 is of order ~1E-7, it is multiplied by 4.557E-21, resulting in the 4.55E-28 term.
          
!          AIRSCAT =  FO2*SIG(2) + FN2*SIG(1) + FAR*SIG(4) 
        ! methane rayleigh scattering (Follows procedure from Sneep and Ubachs, 2004)
         RM = 1 + .00046662 + 4.02E-14*WN**2
         NM = ((RM**2-1)/(RM**2+2))**2
         KCF = (6.+3.*DEL(5))/(6.-7.*DEL(5))
         SIG(5) = 9.75036E-37*NM*KCF*WN**4                 ! = .85*NM*KCF*((24*3.14159**3*WN**4)/(25.47E18**2))

         
        ! Hydrogen rayleigh scattering (From Dalgarno and Williams 1962)
         SIGH = 8.14E-29/(AL**4) + 1.28E-30/(AL**6) + 1.61E-32/(AL**8)

 
         AIRSCAT =  FO2*SIG(2) + FN2*SIG(1) + FAR*SIG(4)+FCH4*SIG(5) 
     &    + FH22*SIGH

     
        ! water rayleigh scattering
        C1 = (.05791817/(238.0185-(1/AL)**2))   ! C1 and C2 expression
        C2 = (.00167909/(57.362-(1/AL)**2))
        RDA = 1 +C1 + C2
        R = .85*(RDA-1)

   
        SIG(6) = 4.577E-21*((6 + 3.*DEL(6))/(6.-7.*DEL(6)))*(R**2/AL**4)
!        SIG(6) = (4.577E-28*(6 + 3.*DEL(6))/(6.-7.*DEL(6))*(R/AL)**4)
        

         SIGR = Fwater*SIG(6) + FCO2*SIG(3) + FNCR*AIRSCAT   ! Water is less than that of air. Current cross-section. REMOVED FNCR from FNCR*AIRSCAT 8/31/2011
!        SIGR = Fwater*SIG(1) + FCO2*SIG(3) + FNCR*AIRSCAT ! Old rayleigh scattering coefficient, assuming water has the same scattering cross-section as N2. 
!        Valid for inner edge comparison with above SIGR only.

!         print *, 'SUMMIX=', Fwater + FCO2 + FNC
c 1135 CONTINUE !EWS - not used

      RETURN
      END  
