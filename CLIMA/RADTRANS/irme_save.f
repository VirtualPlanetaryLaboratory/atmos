

      SUBROUTINE IRME(T,PF,P,FNC,CGAS)
c
c jfk 6/27/08 The call sequence now contains both PF and P. This was
c     being done incorrectly before, because PF was being converted
c     into P during the call.
c
c  This subroutine calculates the infrared flux

c  This subroutine contains CH4 and C2H6 (gna merging in changes from Fung's version that has ethane)
   
 
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NF=55,NGS=8, IK=8)
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5)
C new common block, von Paris, 21/04/2006
      COMMON/IRDATA/WEIGHTCH4(6),xkappa(3,12,NF,8), ! weightch4 is for methane 3/20/2012
     & CIA(7,NF), CPRW(ND,NF)
c-rr !3/23/11 put CIA matrix in IRDATA

      COMMON/VARIR/kappa_irh2o(NF,8,8,IK), kappa_irco2(NF,8,8,IK)! Added kappa matrix in IR for kpsectrum Co2 and H2O coefficients 8/26/2012  
      COMMON/weightsIR/weightco2_h2oIR(IK)
      COMMON/HYDN2CIA/ H2N2CIA(5,55), H2N2FIN(ND,NF) ! c-rr 5/23/12 H2 CIA coefficients from Borysow et al., (2002).
      COMMON/HYDCIA/ H2H2CIA(6,55), H2H2FIN(ND,55) ! c-rr 7/02/2012 H2-H2 CIA coefficients from Borysow et al. 
      COMMON/OXYCIA/O2O2CIA(15,55), O2O2FIN(ND,55) ! c-rr 6/17/2012 O2-O2 CIA coefficients from HITRAN CIA database
      COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
     & ASYAIR(NF,ND-1),IO3,QEXTIR(NF,ND-1)
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND),
     &  FH22
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     2 ALT(ND)
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM, DM2
      COMMON/WAVE/AV(NF),LAM(NF),W(NF)
	  COMMON/BPS_IR/s_abir(NF), f_abir(NF), TDir(NF), Bsir(NF), 
     &    Bfir(NF)  ! Added COMMON BLOCK FOR BPS CONTINUUM FOR IR  8/30/2012

!      REAL KAPPALAYERC2H6(NF,6) !gna
!      COMMON/TFBLOK12/KAPPALAYERC2H6 !gna - what is TFBLOK12? in Fung's version but in our version I think it gets defined below
      
c  
c jfk 6/25/08 Temporary (PF should have been passed in the call
c     statement. Actually, this statement does remain once the code is
c     corrected.)
      DIMENSION PF(ND), MS1(ND),MS(ND), FX(ND),KMS(ND),KMS1(ND),
     &          kmatrix_irh2o(NF,IK),kmatrix_irco2(NF,IK),
     &          MSH1(ND), 
     &          MSH(ND), FXH(ND), FNC(ND),MSO(ND),FXO(ND), 
     &          MSO1(ND), MSHH(ND), MSHH1(ND), FXHH(ND)
 

c     DIMENSION AV(NF)
c      DIMENSION P(ND),BPLANCK(ND),CPR(NF),TPR(NF),TPRIND(ND), 
c     2 TAUGCO2(ND),TAUGH2O(ND),TAUGCH4(ND),TAUGIR(ND),FUP(ND),FDN(ND),
c     3 FUPA(ND),FDNA(ND),T(ND),CGAS(ND,NGS),TAUCONTIN(ND),WEIGHT(8,3),
c     & xkappa(8,12,55,8,3)
      DIMENSION P(ND),BPLANCK(ND),CPR(NF),TPR(NF),TPRIND(ND), 
     2 TAUGCO2(ND),TAUGH2O(ND),TAUGCH4(ND),TAUGIR(ND),FUP(ND),FDN(ND),
     3 FUPA(ND),FDNA(ND),T(ND),CGAS(ND,NGS),TAUCONTIN(ND), TAUH2N2(ND)  ! Added TAUH2 5/29/2012 c-rr
     4 ,TAUO2O2(ND), TAUH2H2(ND),SELF_ABSIR(NF,ND),FORN_ABSIR(NF,ND)  ! Added self and foreign broadening continuum matrices 8/30/2012

      REAL KAPPALAYER(NF,IK,3,ND),kappa(55,6),KAPPALAYEROZ(8,ND), ! redimensioned KAPPALAYER and kappa for CH4 3/21/2012
     2 WEIGHTOZC(8),KAPPAOZC(8),TAUGOZ(ND),TAUCONTINT(NF),CNUT,
     3 TAUTOTAL(NF),TRANSLAYER(ND),TWGHTT(8),OMG0IR(ND-1),ASYIR(ND-1),
     4 TAULAMIR(ND-1),TAUAEXTIR(ND-1),TAUASIR(ND-1),TAUSIR(ND-1)
      REAL KAPPALAYERC2H6(NF,6) !c2h6
      REAL DPLYR(ND), PLAYR(ND), CRAY(ND), TAUR(ND)
      DATA HP,SIGMA/6.63E-27, 5.67E-5/
      DIMENSION WEIGHTC2H6(6),CGASC2H6(ND) !c2h6
      DIMENSION TAUGC2H6(ND)
      REAL KAPPAC2H6(NF), kmatrix_irh2o, kmatrix_irco2, kappa_ir
      REAL np
c-rr Added CIAMS1L, CIAMSL, and CPRL and CIA Common block items 3/24/2011
      REAL TAUSUM(6), CIAMS1L, CIAMSL, CPRL, CIA, CPRW
      

C   PRESSURE-INDUCED CO2 ABSORPTION (FROM JIM POLLACK)
      DATA CPR/4.3E-5, 3.8E-5, 1.2E-5, 2.8E-6, 7.6E-7, 4.5E-7, 2.3E-7,
     2  5.4E-7, 1.6E-6, 10*0., 4.E-7, 4.E-6, 1.4E-5, 1.0E-5,  ! 4e-7 value updates 7.5e-7 in Kasting et al. (1984)
     3  1.2E-6, 2.0E-7, 5.0E-8, 3.0E-8, 28*0./



c-jdh No 20-um band
c      DATA CPR/0, 0, 0, 0, 0, 0, 0,
c     2  0, 0, 10*0., 4.E-7, 4.E-6, 1.4E-5, 1.0E-5,
c     3  1.2E-6, 2.0E-7, 5.0E-8, 3.0E-8, 28*0./
c-jdh No 7.7-um band
c      DATA CPR/4.3E-5, 3.8E-5, 1.2E-5, 2.8E-6, 7.6E-7, 4.5E-7, 2.3E-7,
c     2  5.4E-7, 1.6E-6, 10*0., 0, 0, 0, 0,
c     3  0, 0, 0, 0, 28*0./

      DATA TPR/3.4, 2.2, 1.9, 52*1.7/

C   FREQUENCIES AT ENDS OF SPECTRAL INTERVALS (1/CM)
c-jdh commented out, declared common in Clima.f, converted to units of (1/S)
c      DATA AV/40., 100., 160., 220., 280., 330., 380., 440., 495.,
c     2  545., 617., 667., 720., 800., 875., 940., 1000., 1065.,
c     3  1108., 1200., 1275., 1350., 1450., 1550., 1650., 1750., 1850.,
c     4  1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760.,
c     5  4030., 4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315.,
c     6  8850., 9350., 9650., 10400., 11220., 11870., 12790., 13300.,
c     7  14470., 15000./
C
c-jdh Ethane data (for polynomial fit approximation)
      DATA KAPPAC2H6/0,0,0,0,0,0,0,0,0,0,
     &  0,0,0,8.627161E-21,1.036469E-20,1.098662E-21,0,0,0,0,
     &  0,2.377869E-22,9.640977E-21,1.726593E-20,1.426406E-21,0,
     &  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &  0,0,0,0,0,0,0,0,0,0/

c-jdh Ethane k-coefficient weights
      DATA WEIGHTC2H6/0.08566,0.18038,0.23396,0.23396,0.18038,0.08566/ !c2h6

       INTEGER COUNTERIR

       
       COUNTERIR = 0
       SRFALBIR = 0.
       NLAYERS = ND - 1    
c      BCON = 2.*HP/C/C
       HK = HP/BK
C       np = 1.E+1

       
       C2 = 1.4388 ! constant for BPS continuum (in cmK) 8/31/2012 , c-rr
        
c-jdh Set C2H6 absorption coefficients
      DO NI=1,NF
        DO NK=1,6
          KAPPALAYERC2H6(NI,NK) = 0.0       
         END DO
      END DO
!so these are like 10-20 orders of magnitude higher than the kappas I get for co2 and ch4...
!what are these values and what are their units 
!gna - the 2.68E+19 is something I'm testing out.  I don't know how kappalayerc2h6 was calculated,
!but it's resulting in tau values that are obviously much too large (and larger than those of all the other species
!even when the amount of c2h6 in the atmosphere is miniscule).  In gascon.f, cgas in the "new" version (from 2012?)
!is divided by 2.68E+19 to convert to molecules/cm**2, which is not something gasconold.f did.  Possibly this kappalayerc2h6
!is assuming units we used pre-2012?
!how were the methane k coeffs calculated?
      KAPPALAYERC2H6(14,1) = 1.0140E-23/2.687E+19
      KAPPALAYERC2H6(14,2) = 6.6404E-23/2.687E+19
      KAPPALAYERC2H6(14,3) = 4.4811E-22/2.687E+19
      KAPPALAYERC2H6(14,4) = 2.2630E-21/2.687E+19
      KAPPALAYERC2H6(14,5) = 5.1880E-21/2.687E+19
      KAPPALAYERC2H6(14,6) = 7.8162E-21/2.687E+19
      KAPPALAYERC2H6(15,1) = 4.1940E-21/2.687E+19
      KAPPALAYERC2H6(15,2) = 5.3313E-21/2.687E+19
      KAPPALAYERC2H6(15,3) = 7.0534E-21/2.687E+19
      KAPPALAYERC2H6(15,4) = 9.2389E-21/2.687E+19
      KAPPALAYERC2H6(15,5) = 1.3508E-20/2.687E+19
      KAPPALAYERC2H6(15,6) = 2.7746E-20/2.687E+19
      KAPPALAYERC2H6(16,1) = 2.3159E-23/2.687E+19
      KAPPALAYERC2H6(16,2) = 8.6746E-23/2.687E+19
      KAPPALAYERC2H6(16,3) = 2.6336E-22/2.687E+19
      KAPPALAYERC2H6(16,4) = 9.1975E-22/2.687E+19
      KAPPALAYERC2H6(16,5) = 2.3556E-21/2.687E+19
      KAPPALAYERC2H6(16,6) = 3.9631E-21/2.687E+19
      KAPPALAYERC2H6(22,1) = 3.0909E-24/2.687E+19
      KAPPALAYERC2H6(22,2) = 1.5773E-23/2.687E+19
      KAPPALAYERC2H6(22,3) = 4.9806E-23/2.687E+19
      KAPPALAYERC2H6(22,4) = 1.3673E-22/2.687E+19
      KAPPALAYERC2H6(22,5) = 3.3965E-22/2.687E+19
      KAPPALAYERC2H6(22,6) = 8.2701E-22/2.687E+19
      KAPPALAYERC2H6(23,1) = 1.0923E-21/2.687E+19
      KAPPALAYERC2H6(23,2) = 2.7999E-21/2.687E+19
      KAPPALAYERC2H6(23,3) = 7.0277E-21/2.687E+19
      KAPPALAYERC2H6(23,4) = 1.1359E-20/2.687E+19
      KAPPALAYERC2H6(23,5) = 1.5960E-20/2.687E+19
      KAPPALAYERC2H6(23,6) = 2.1332E-20/2.687E+19
      KAPPALAYERC2H6(24,1) = 6.0680E-21/2.687E+19
      KAPPALAYERC2H6(24,2) = 9.5224E-21/2.687E+19
      KAPPALAYERC2H6(24,3) = 1.3251E-20/2.687E+19
      KAPPALAYERC2H6(24,4) = 1.6683E-20/2.687E+19
      KAPPALAYERC2H6(24,5) = 2.1862E-20/2.687E+19
      KAPPALAYERC2H6(24,6) = 4.6329E-20/2.687E+19
      KAPPALAYERC2H6(25,1) = 3.1607E-23/2.687E+19
      KAPPALAYERC2H6(25,2) = 7.4713E-23/2.687E+19
      KAPPALAYERC2H6(25,3) = 2.2452E-22/2.687E+19
      KAPPALAYERC2H6(25,4) = 9.2931E-22/2.687E+19
      KAPPALAYERC2H6(25,5) = 3.0558E-21/2.687E+19
      KAPPALAYERC2H6(25,6) = 5.7388E-21/2.687E+19
      DO I=1,NLAYERS
!gna: added ethane to cgas in gascon.f so these lines no longer needed
c        CGASC2H6(I) = 1.3386e-3*BK*273.16/(SM*G)
c     &                * (P(I+1) - P(I))*FC2H6/DM
c        CGASC2H6(I) = 1.0e-5*BK*273.16/(SM*G)
c     &                * (P(I+1) - P(I))*FC2H6/DM
c jfk 6/27/08 P was changed to PF in the 2 lines below.
!       CGASC2H6(I) = 1.0e-5*BK*273.16/(SM*G) !see no factor of 2.68e19
!     &                * (PF(I+1) - PF(I))*FC2H6/DM !c2h6

      END DO


! Intializing Kappalayer, temperature, and FX indices    4/25/2012 c-rr
        DO IL = 1, NLAYERS
          DO J = 1,IK
            DO I = 1,NF   !frequencies?
              DO K = 1,3 ! Gases are H2O(gas 1), CO2(gas 2), and CH4(gas 3)
        KAPPALAYER(I,J,K,IL) = 1.e-60
              ENDDO
            ENDDO
          ENDDO
          MS(IL) =0.0
          MS1(IL)= 0.0
          FX(IL) = 0.0
          MSH(IL)=0.0
          MSH1(IL)=0.0
          FXH(IL)=0.0
          MSO(IL)=0.0
          MSO1(IL)=0.0
          FXO(IL)=0.0
          MSHH(IL)= 0.0
          MSHH1(IL)=0.0
          FXHH(IL)=0.0
        ENDDO

! Initializing kmatrix

       DO I = 1, IK
           DO J = 1,NF
           kmatrix_irh2o(J,I) = 0.0d0
           kmatrix_irco2(J,I) = 0.0d0
          ENDDO
      ENDDO


	DO JS = 1, ND
	PLAYR(JS) = PF(JS)*1.E6      !PF IN BARS; PLAYR IN DYNES/CM^2
	ENDDO

C Read the IR exponential sums
c-jdh Moved to Clima.f 
c      CALL IREXPSUMS(WEIGHT,xkappa)

      DO 7 IL = 1,NLAYERS  !
       
        TIL = AMIN1(T(IL),600.)
        TIL = AMAX1(TIL,100.)
        

       CALL INTERPIR(TIL,P(IL), IL,kappa,
     &  kmatrix_irco2, kmatrix_irh2o) ! Outputting kmatrix_irco2 and kmatrix_irh2o 8/27/2012
!      print *, 'calling interpir'         
 
       CALL INTERPCO2CIA(T(IL),IL,MS(IL),MS1(IL),FX(IL))
!         KMS(IL)  = MS(IL)
!         KMS1(IL) = MS1(IL)
       
        
       CALL INTERPH2N2CIA(T(IL),IL,MSH(IL),MSH1(IL),FXH(IL))

       CALL INTERPO2CIA(T(IL),IL,MSO(IL),MSO1(IL),FXO(IL))

       CALL INTERPH2H2CIA(T(IL),IL,MSHH(IL),MSHH1(IL),FXHH(IL))

          
 !      print *, 'MS=', MS(IL), IL
             DO I =1, 55
                  
              DO J = 1, IK
         KAPPALAYER(I,J,1,IL) = kmatrix_irh2o(I,J) ! Gas 1 is H2O
         KAPPALAYER(I,J,2,IL) = kmatrix_irco2(I,J) ! Gas 2 is CO2

 

! 3383 format(1pe14.5,3x,0p,a5,3x, i0, 3x, a6, i0, 4x, a6,4x,
!     &  1pe12.5, 0p,4x,a5,4x,f8.2)  
           ENDDO

           ENDDO
            
!            DO J = 1,16
!                KAPPALAYER(33,J,1,IL) = KAPPALAYER(32,J,1,IL)   ! Address numerical instability in interval 33 with negative fluxes
!            ENDDO
           
        
        DO 8 I=1, 55  ! Methane loop
        DO 9 J=1, 6   ! number of sums, coefficients
        KAPPALAYER(I,J,3,IL) = kappa(I,J) ! Gas 3 is CH4
         
        
 9      CONTINUE
 8      CONTINUE
 7      CONTINUE  ! We are precalculating the INTERPOLATION for both INTERPIR AND INTERPCO2CIA for all layers 3/21/2012

       
     


!      DO IL = 1, NLAYERS
c     Added new KAPPALAYER loop for kmatrix_IR (which is CO2/H2O) 3/20/2012
!      DO I =1, 55
!        DO J = 1, 16
!         KAPPALAYER(I,J,1,IL) = kmatrix_ir(I,IL,J) ! Gas 1 is CO2/H2O
         

!        ENDDO
!      ENDDO       
!      pause
!       DO 8 I=1, 55  !K species
!         DO 9 J=1, 8   ! number of sums, coefficients
!         KAPPALAYER(I,J,2,IL) = kappa(I,J) ! Gas 2 is CH4
 
!  9   CONTINUE
!  8   CONTINUE
      
         
!      ENDDO
3333  format(1p1e14.5,2x,i3,2x,i3,2x,i3)
c-------------------------------
     

      DO IL = 1,NLAYERS

       CALL INTERPOZONE(P(IL),WEIGHTOZC,KAPPAOZC)
        DO J = 1, 8
        KAPPALAYEROZ(J,IL) = KAPPAOZC(J)
        ENDDO 
        
      ENDDO
c       PRINT*, T,P
       DO 99 I=1, ND
        FUPIR(I) = 0.
        FDNIR(I) = 0.
 99   CONTINUE
       DO I=1,55
         TAUTOTAL(I) = 0.0
       ENDDO


        write(5552,2333)
2333   format(3x, 'TAUGH2OCO101', 3x, 'PATH_L', 9x, 'Temp@top',
     &    8x, 'KAPPA', 9x, 'FH2O', 9x,'INT', 4x, 'NST')

****** Loop over frequency      
       DO 1 I=1, NF
                   AL2   =  (1E4*(C/AV(I)))**2   ! wavenumber 
                   VAC = AV(I)
       DO 20 J=1,ND 
C-jdh  VAC = C*AV(I)


!        if (J.eq.1) then
!        print *, AV(I)/C
!        endif
       BPLANCK(J) = PLANCK(VAC,T(J),HP,C,HK)
  20   CONTINUE
       DO 15 J=1,ND
         FUPA(J) = 0.0
         FDNA(J) = 0.0
         FUP(J) = 0.0
         FDN(J) = 0.0
         TRANSLAYER(J) = 0.0 
  15   CONTINUE
****** AEROSOLS

       DO IL=1,NLAYERS
        r = RAER(IL)
        np = PARTICLES(IL)
        TAUAEXTIR(IL) = QEXTIR(I,IL)*PI*r*r*DALT(IL)*np*1.E+5
        TAUASIR(IL) = OMG0AIR(I,IL)*TAUAEXTIR(IL)
         
       ENDDO
       IF (I.EQ.15) THEN
       TAUEXTIRTOTAL = 0.
       TAUASIRTOTAL = 0.
       DO IL=1,NLAYERS
       TAUEXTIRTOTAL = TAUEXTIRTOTAL + TAUAEXTIR(IL)
       TAUASIRTOTAL = TAUASIRTOTAL + TAUASIR(IL)
       ENDDO
c       PRINT *, '******************************'
c       PRINT *, 'TAUEXTIRTOTAL'
c      PRINT *, TAUEXTIRTOTAL
c       PRINT *, 'TAUSIRTOTAL'
c       PRINT *, TAUASIRTOTAL
       TAUAABSIRTOTAL = TAUEXTIRTOTAL - TAUASIRTOTAL
c       PRINT *, 'TAUAABSIRTOTAL'
c       PRINT *, TAUAABSIRTOTAL
c       PRINT *, '*******************************'
       ENDIF 

!----------------BPS WATER CONTINUUM 8/30/2012 c-rr
        ! PF(IL) are pressures at the layer boundaries. P(IL) are pressures in the middle of the layers


 !          IF (I.eq.15)sumcont = 0.0d0
 
  !     IF(((I.ge.15).and.(I.le.20)).or.((I.ge.27).and.(I.le.33)))THEN		
  !      IF((I.ge.15).and.(I.le.20))THEN		
!	IF(I.eq.15)THEN

           DO IL = 1, NLAYERS

                AAV  = AV(I)/C ! convert from frequency to wavenumber (taken at interval midpoints).
                PH2O = P(IL)*FI(1,IL)  ! water partial pressure (in bars)
                PDRY = P(IL)*(1.-FI(1,IL)) !partial pressure dry air (in bars)

                
                RHOW = (PH2O/1.)*(296./T(IL))  ! RHO_water. Different from MT_CKD? Why does BPS use Ph2o instead of P? 
                RHOF = (PDRY/1.)*(296./T(IL))  ! RHO_foreign. Different from MT_CKD? Why does BPS use Pdry instead of P?

               
                ANUM = (exp((C2*AAV)/296.)-1.)*(exp((C2*AAV)/T(IL))+1.)  
                DEN = (exp((C2*AAV)/296.)+1.)*(exp((C2*AAV)/T(IL))-1.)
                RADFLD = ANUM/DEN ! radiation field hyperbolic expression written in terms of exponents 

!         Bsir(i) = 0.0d0
!		 Bfir(i) = 0.0d0
 
             if(AAV.ge.500)RADFLD = 1.  ! radiation field term is negligible above 500cm-1.
          self_absir(I,IL) = RHOW*(s_abir(i)*exp(TDir(i)*(296.-T(IL)))
     &        + RADFLD*Bsir(i)) ! Self broadening coefficient with radiation field and temperature dependence included (cm^2/molecule)
               forn_absir(I,IL) = RHOF*(f_abir(i) + Bfir(i))*RADFLD  ! foreign broadening coefficient with radiation field and temperature dependence included(cm^2/molecule)

                ABSCONT = self_absir(I,IL) + forn_absir(I,IL)
                TAUCONTIN(IL) = ABSCONT*CGAS(IL,6) !DP*FI(1,IL)/(DM*SM*G) !optical depth
!                TAUCONTIN(IL) = 0.0d0

!            if(T(ND)>2200.0d0)then
!              print 4444,self_absir(I,IL),forn_absir(I,IL),
!     &                    FI(1,IL),IL,I
!            endif
4444        format(1p3e14.5,2(2x,i3))
!        print *, TAUCONTIN(IL), DP, FI(1,IL), DM, SM, G, 
!     &   PF(IL), PF(IL+1),P(IL)   ! DIFFERENCE BETWEEN PF and P?
!                 print *, ABSCONT, DP*FI(1,IL)/(DM*SM*G), CGAS(IL,6)
!		pause
         
    !    IF ((I.eq.15).and.(IL.ge.97))sumcont = sumcont + TAUCONTIN(IL)   ! only add approximately 1 km. tau ~ 0.1 or so for Earth
                  ENDDO  ! ENDS LAYER LOOP IN CONTINUUM
    
!	    ENDIF
   !      IF (I.eq.15)print *, sumcont  
!---------------------------------------------------------	   






         sumtwght = 0.0
         SUM_TPRIND=0.   ! summing TPRIND experiment. Sets SUMTPRIND when start a new frequency (I)
         SUM_TAU=0.
         SUM_CO2H2O = 0.

!        print *, VAC, AV(I), I

       DO 2 K1 = 1,IK ! H2O loop
       Do 3 K2 = 1,IK ! CO2 loop
       DO 4 K3 = 1,6 ! methane loop  ! Activated when methane is on (IMET = 1)
       DO 5 K0 = 1,6 ! ethane loop
  
!        TWGHT = weightco2_h2oIR(K1)*weightco2_h2oIR(K2) ! no methane (CO2 and H2O weights are the same) 8/27/2012
 
!          print *, weightco2_h2oIR(K1)
!          print *, weightco2_h2oIR(K2)
!          print *, weightch4(k2)
!          print *, weightch4(k3)
!          pause
       	
!       TWGHT = weightco2_h2oIR(K1)*weightco2_h2oIR(K2)*weightch4(K3) ! with methane 8/27/2012
        TWGHT = weightco2_h2oIR(K1)*weightco2_h2oIR(K2)*weightch4(K3) 
     &           *weightc2h6(K0) 
!        TWGHT = TWGHT * weightc2h6(K4) !with ethane


        DO 11 IL = 1,NLAYERS
         PPE = (1. + 0.3*FI(2,IL))*P(IL)     !CO2
         TPE = (300./T(IL))**TPR(I)
         
         TAUGH2O(IL) = KAPPALAYER(I,K1,1,IL)*CGAS(IL,6)
!         print *, cgas(IL, 6)
!         call sleep(1)
!   1145011266883586.5     
!   1842019049386489.0     
!   0.0000000000000000     
!   90771608043177040.     
!   15.106465570752329     
!   3951422779566604.5     
!   0.0000000000000000  

         TAUGCO2(IL) = KAPPALAYER(I,K2,2,IL)*CGAS(IL,5)
        ! if(taugco2(il) .gt. 0) print *, taugco2(il)
        ! print *, 'CGAS(IL, 5) CO2', CGAS(IL, 5)

         TAUGCH4(IL) = KAPPALAYER(I,K3,3,IL)*CGAS(IL,2)!! Activated when methane is on (IMET = 1)
          !if(taugch4(il) .gt. 0) print *, taugch4(il)
      !   if (kappalayer(i,k3,3,il) .gt. 0) then 
      !      print *, kappalayer(I,k3,3,il)  !kappalayer for ch4 much smaller than t
      !      endif
        ! print *, 'CGAS(IL, 2) CH4 ', CGAS(IL, 2)

             
!4242     format(1p3e14.5,2x,i3,2x,i3,2x,i3)
 4242     format(1p7e14.5, 4x, i3)

!         TAUGC2H6(IL) = 0. !gna - commented out
!         TAUGC2H6(IL) = KAPPAC2H6(I)*CGAS(IL,5)*2.687E24*(FC2H6/FCH4)  !gna - is THIS correct?
!         TAUGC2H6(IL) = KAPPALAYERC2H6(I,K4)*CGASC2H6(IL)*2.687E24 !gna - from Fung's version (is this right??)
          TAUGC2H6(IL) = KAPPALAYERC2H6(I,K0)*CGAS(IL, 8) !gna - gas 8 is ethane
         !  if(taugc2h6(il) .gt. 0) print *, taugc2h6(il)
         ! if(taugc2h6(il) .gt. 0)print *, 'taugc2h6', TAUGC2H6(IL)
!          if(KAPPALAYERC2H6(I,K0) .gt. 0) print *, 
!     &       'KAPPALAYERC2H6(I,K0)', KAPPALAYERC2H6(I,K0)
         ! print *, 'CGAS(IL, 8) C2H6', CGAS(IL, 8)
        !  call sleep(1)
          
         ! if (kappalayerc2h6(i,k4) .gt. 0) print *, kappalayerc2h6(I,k4)
          

   !       print *, cgas(IL, 8)
   !       print *, 'hi'
    !      call sleep(1)
          !try this:
!          TAUGC2H6(IL) = KAPPAC2H6(I)*CGAS(IL,8)*2.687E24  !need to multiply by 2.68E24?
!          print *, cgas(IL, 8)
!         print *, H2CIA(MS1(IL),I)
!         pause
c-rr	This the H2-N2 CIA loop calculation 5/29/2012
        H2N2CIAMS1L = log(H2N2CIA(MSH1(IL),I))
        H2N2CIAMSL  = log(H2N2CIA(MSH(IL),I))
        H2N2TOT = H2N2CIAMS1L + FXH(IL)*(H2N2CIAMSL - H2N2CIAMS1L)
        H2N2FIN(IL,I) = exp(H2N2TOT)


c-rr	This the H2-H2 CIA loop calculation 7/02/2012
        H2H2CIAMS1L = log(H2H2CIA(MSHH1(IL),I))
        H2H2CIAMSL  = log(H2H2CIA(MSHH(IL),I))
        H2H2TOT = H2H2CIAMS1L + FXHH(IL)*(H2H2CIAMSL - H2H2CIAMS1L)
        H2H2FIN(IL,I) = exp(H2H2TOT)


          
c-rr	This is the O2 CIA loop calculation 6/17/2012
        O2O2CIAMS1L = log(O2O2CIA(MSO1(IL),I))
        O2O2CIAMSL  = log(O2O2CIA(MSO(IL),I))
        O2O2TOT = O2O2CIAMS1L + FXO(IL)*(O2O2CIAMSL - O2O2CIAMS1L)
        O2O2FIN(IL,I) = exp(O2O2TOT)


c-rr	This is the CO2 CIA loop that takes the saved values for FX, MS, and MS1 at a given height and calculates
c-rr	the correct interpolated CIA values. These are from the model of Wordsworth et al. (2010) that uses the GBB
c-rr	parametrization scheme. 3/24/11
c        print *, 'is this working?'
        
!        print *,MS1(IL),MS(IL),IL,I
!         pause
!        print *,CIA(MS1(IL),I)
	CIAMS1L = log(CIA(MS1(IL),I))
	CIAMSL = log(CIA(MS(IL),I))
	CPRL = CIAMS1L + FX(IL)*(CIAMSL-CIAMS1L)
	CPRW(IL,I) = exp(CPRL)


       IF (CPRW(IL,I).lt.1.E-45)CPRW(IL,I)= 0. !c-rr Ensures bands where CIA is supposed to be zero, are zero. 5/28/2011
       
       

       IF(I .EQ. 38)THEN
          CPRW(IL,I) = 3.5E-8  !Rewriting CIA at intervals 37 and 40 with appropriate CIA at 2.3 and 1.73 microns, respectively(Tsang et al. 2008).
         
       ELSEIF (I.EQ.41)THEN
          CPRW(IL,I) = .57*6.0E-9  !Because the 1.73 micron window is only 57% of this bin
       ELSEIF ((I.EQ.45).OR.(I.EQ.46).OR.(I.EQ.47))THEN  !Rewriting CIA at intervals 45-47 with 1.2 micron complex
	      CPRW(IL,I) = 1.5E-9
               
       ENDIF
       
   

****** Pressure induced absorption by CO2 
         CGAS1 = CGAS(IL,5)/2.687E19  ! Gas 5 in CGAS is CO2  3/20/2012. Converts mol/cm^2 into atm-cm. 3/30/2012
c-rr 3/25/11 commenting out Kasting et al. (1984) parametrization

c-rr 3/24/11 	 PUT NEW TPRIND parametrization here! 
          PCGS = P(IL)* 1.e6
!       TPRIND(IL) = CPRW(IL,I)*(PCGS/(BK*T(IL)*2.687E19))*CGAS1
!    &   *((1/1.3)+(FI(2,IL)/1.3))

!  c-rr Corrected bug in CO2 CIA. Removed extraneous pressure broadening factor of (1/1.3). 2/7/2012

        TPRIND(IL) = CPRW(IL,I)*(PCGS/(BK*T(IL)*2.687E19))*CGAS1
     &   *(1.+0.3*FI(2,IL))
!        print *, BK,CGAS1,FI(2,IL), IL
!        pause

!        IF(IL.eq.100) THEN
!        print *, 'TPRIND=', TPRIND(IL), 'CGAS1=', CGAS1
!        ENDIF
!         TPRIND(IL) = 0.  ! zeroes out CIA

	SUM_TPRIND = SUM_TPRIND + TPRIND(IL) !summing TPRIND experiment
        
         
c          print *, TPRIND(IL),IL
c        TPRIND(IL) = 0
*******
       
!-------------------Pressure-induced absorption by N2-H2 and then add 1 +0.3FCO2 to simulate CO2-H2

         CGASH2 = CGAS(IL,7)/2.687E19 ! Converts molec./cm^2 into atm-cm.

         TAUH2N2(IL) = H2N2FIN(IL,I)*(PCGS/(BK*T(IL)*2.687E19))
     &   *CGASH2*(1.-(FH22*FNC(IL)))*(1. + 0.3*FI(2,IL))
!----------------------------------------------------------------
!         if ((K1.eq.1))then
!         print *, TAUH2N2(IL), FH22*FNC(IL)
!         pause
!         endif


!-------------- Pressure-induced absorption by O2-O2 
        CGASO2 = CGAS(IL,3)/2.687E19 ! Converts molec./cm^2 to atm-cm

       
        TAUO2O2(IL) = O2O2FIN(IL,I)*(PCGS/(BK*T(IL)*2.687E19))
     &  *CGASO2*FO2*FNC(IL)
!               TAUO2O2(IL)=1.e-60
      

!-------------- Pressure-induced absorption by H2-H2 
 
        TAUH2H2(IL) = H2H2FIN(IL,I)*(PCGS/(BK*T(IL)*2.687E19))
     &  *CGASH2*FH22*FNC(IL)

         
!                  if (K1.eq.1)then
!            print *, H2H2FIN(IL,I), TAUH2H2(IL), FH22, FNC(IL), CGASH2
!     &      , CGAS1, TAUGIR(IL), T(IL), IL
!            pause
!                 endif       

!----------------------------------------------------------

           TAUH2H2(IL) = AMAX1(TAUH2H2(IL),1.E-60)
           TAUO2O2(IL) = AMAX1(TAUO2O2(IL),1.E-60)
           TAUH2N2(IL) = AMAX1(TAUH2N2(IL),1.E-60)


         TAUGIR(IL) = TAUGH2O(IL)+ TAUGCO2(IL) + TAUGCH4(IL)
     &                 + TPRIND(IL)+
     &                TAUH2N2(IL)+TAUCONTIN(IL)+TAUGC2H6(IL)+ 
     &                TAUO2O2(IL)+TAUH2H2(IL)

!          print 2222,TAUGH2O(IL),TAUGCO2(IL),IL,K1,K2
!           if ((K1.eq.1).and.(I.eq.1))then
!            print *, FNC(IL)
!     &      , CGAS1, TPRIND(IL), TAUGIR(IL), T(IL), 
!     &        TAUGH2O(IL),TAUGCO2(IL),IL
!            pause
!          endif         

 
!            if((K1.eq.16).and.(IL.eq.1))then
!             print 444, TAUGCH4(IL), TAUGH2O_CO2(IL), TPRIND(IL), 
!     &              TAUCONTIN(IL), TAUGC2H6(IL), I
!            pause
!            endif




!        if(NST==6)then
!        print 3321,TAUGH2O_CO2(IL),real(IL),real(I)
!         print *,TPRIND(IL),
!           IF (IL.eq.100)THEN
!          print 2222,TAUGH2O_CO2(IL),TAUGCH4(IL),TPRIND(IL),
!     &    TAUCONTIN(IL),TAUGC2H6(IL), IL
 2222     format(1p2e14.5, 4(2x,i3))
!          ENDIF

             


  11    CONTINUE
444      format(1p5e14.5,2x,i3)	
        
       
        
!        GOTO 55555
!	print *, AV(I), SUM_TPRIND,I
!        read(*,*)
                        
******* Ozone absorption
C-TF  THIS SPECIAL TREATMENT FOR OZONE IS DONE TO SAVE COMPUTER
C-TF  TIME BECAUSE ONLY ONE WAVELENGTH IS AFFECTED BY OZONE.
C-TF  HERE A NEW WEIGHT FUNCTION TWGHTT(K00) IS CALCULATED
C-TF  AND IS USED TO COMPUTE FUPA AND FDNA.
        IF (I.EQ.18) THEN
!		 sumoz = 0.  
      DO K4=1,8
        TWGHTT(K4) = TWGHT*WEIGHTOZC(K4)
      DO IL=1,NLAYERS
!        CGAS(IL,4) = CGAS(IL,4)*.97d0
         
!		  sumoz = TWGHTT(K4) + sumoz
        TAUGOZ(IL) =KAPPALAYEROZ(K4,IL)*CGAS(IL,4) ! Rederived correct units for TAUGOZ rr and rv
!		TAUGOZ(IL) = 0.0d0

!        TAUGOZ(IL) = 0.
!        print *, TAUGOZ(IL),KAPPALAYEROZ(K4,IL),CGAS(IL,4),IL,I
!        pause
c       TAUGOZ(IL) =KAPPALAYEROZ(K4,IL)*CGAS(IL,4)*(4.46E-5*DM)/(48.*SM)
c -PJK There was a bug in the line of code below. The IR optical depth was
c      being compounded improperly in the ozone band. Its done right below.
c      TAUGIR(IL) =TAUGIR(IL) + TAUGOZ(IL)
      TAUGIR(IL) = TAUGCH4(IL)+TAUGH2O(IL)+ TAUGCO2(IL) + TPRIND(IL)  
     & +TAUCONTIN(IL) + TAUGOZ(IL) + TAUGC2H6(IL)+ TAUH2N2(IL)
     & + TAUO2O2(IL)+TAUH2H2(IL)
!gna - taugir is computed twice??  I guess it just overwrites taugir that's computed above...

	! print *, TAUGOZ(IL)
	! pause
	       ENDDO
      

      DO IL=1, NLAYERS
           Fwater = FI(1,IL) ! Needed for rayley
	   FCO2 = FI(2,IL)  ! Needed for rayley
           FNCR = FNC(IL) ! Needed for rayley
           DPLYR(IL)=PLAYR(IL+1)-PLAYR(IL)
           AM = DM
           CONS0=6.0255E23/GNEW(IL)   ! Put constant to use for rayleigh scattering calculation in IR
           CONS=CONS0/AM
           CRAY(IL)=CONS*DPLYR(IL)
           CALL RAYLEY(SIGR, AL2,Fwater,FNCR) ! Call rayley to output SIGR at a given altitude and wavelength, inputting AL2. 5/8/2011
           TAUR(IL)= SIGR*CRAY(IL)  ! SIGR is a scalar now 5/8/2011
           TAUSIR(IL) = TAUASIR(IL) + TAUR(IL)


          TAULAMIR(IL) = TAUAEXTIR(IL) + TAUGIR(IL) + TAUSIR(IL)
          ASYIR(IL) = ASYAIR(I,IL) 
          OMG0IR(IL) = TAUSIR(IL)/TAULAMIR(IL)
          OMG0IR(IL) = AMIN1(OMG0IR(IL),0.99999)
          OMG0IR(IL) = AMAX1(OMG0IR(IL),1.E-5)  ! changed lower limit of single-scattering albedo to 1e-12 (from 1.e-5) c-rr 4/30/2012
      ENDDO

!       print *,'k1.....',K1

c jfk  6/25/08  Include TAUTOP in the call sequence to do the upper BC.
c      This requires the use of PF, not P.
      TAUTOP = TAULAMIR(1)*PF(1)/(PF(2)-PF(1))
      TAUTOP = AMIN1(TAUTOP,1.)


!               if((I.eq.19).and.(K4.eq.8))then
!               print *, 'going into DELTATWOSTRIR', TAULAMIR(1)
!               pause
!		endif
      
      CALL DELTA2STRIR(SRFALBIR,ASYIR,TAULAMIR,OMG0IR,
     & FUP,FDN,BPLANCK,TAUTOP, I ,K1, IL)
C

                  DO J = 1,ND
                  FUPA(J)=FUPA(J)+TWGHTT(K4)*FUP(J)
                  FDNA(J)=FDNA(J)+TWGHTT(K4)*FDN(J)
                  ENDDO
        

        ENDDO
		 !  IF (I.eq.18) print *, sumoz
        GOTO 4 ! 3/20/2012 changed from GOTO 5 to GOTO 4
        ENDIF
***************

          
!55555      CONTINUE


      DO IL=1, NLAYERS
           Fwater = FI(1,IL) ! Needed for rayley
	   FCO2 = FI(2,IL)  ! Needed for rayley
           FNCR = FNC(IL) ! Needed for rayley
           DPLYR(IL)=PLAYR(IL+1)-PLAYR(IL)
           AM = DM
           CONS0=6.0255E23/GNEW(IL)   ! Put constant to use for rayleigh scattering calculation in IR
           CONS=CONS0/AM
           CRAY(IL)=CONS*DPLYR(IL)
          CALL RAYLEY(SIGR, AL2,Fwater,FNCR) ! Call rayley to output SIGR at a given altitude and wavelength, inputting AL2. 5/8/2011
           TAUR(IL)= SIGR*CRAY(IL)  ! SIGR is a scalar now 5/8/2011
           TAUSIR(IL) = TAUASIR(IL) + TAUR(IL)
          TAULAMIR(IL) = TAUAEXTIR(IL) + TAUGIR(IL) + TAUSIR(IL)

             
!                       if ((I.eq.20).and.(IL.eq.1).and.(K1.eq.16))then
!            print *, 'TAUOZ=',TAUGOZ(IL),'TAUGIR=',
!     &              TAUGIR(IL), 'TAUL=', TAULAMIR(IL)
!           endif


          ASYIR(IL) = ASYAIR(I,IL) 
          OMG0IR(IL) = TAUSIR(IL)/TAULAMIR(IL)
!          if (I.eq.49) then
!          print 23656, TAUSIR(IL), TAULAMIR(IL), OMG0IR(IL), T(IL), IL
!            pause
!          endif
!         
          
          OMG0IR(IL) = AMIN1(OMG0IR(IL),0.99999)
          OMG0IR(IL) = AMAX1(OMG0IR(IL),1.E-5)! changed lower limit of single-scattering albedo to 1e-12 (default 1e-5) c-rr 4/30/2012
!          if(NST==6)then
!          print 3111,OMG0IR(IL),ASYIR(IL),real(IL),real(I)
!          endif
!           print *,OMG0IR(IL),IL,K1,K2,I


!                if(T(ND) >2200.0d0)then
!                print 23656,TAUGIR(IL),TAUGCH4(IL),TAUGOZ(il),
!     &                 TAUH2N2(IL),TAUO2O2(IL),TAUH2H2(IL),
!     &                 TAUCONTIN(IL),TPRIND(IL),I,K1,IL
!                endif
23656	  format(1p8e14.5,0p,3(2x,i3))

      ENDDO



!               if((I.eq.10).and.(K1.eq.16))then
!              print *, 'going into DELTATWOSTRIR', TAULAMIR(1),
!     &           TAUGIR(1), TAUSIR(1)
!               pause
!		endif




!            if((K1.eq.16).and.(IL.eq.1))then
!             print 444, TAUGCH4(IL), TAUGH2O_CO2(IL), TPRIND(IL), 
!     &              TAUCONTIN(IL), TAUGOZ(IL), TAUGC2H6(IL), I
!            pause
!            endif


!          if (I.eq.49)then
!          print *, I
!          endif
!          pause
c
c jfk  6/25/08  Include TAUTOP in the call sequence to do the upper BC
      TAUTOP = TAULAMIR(1)*PF(1)/(PF(2)-PF(1))
      TAUTOP = AMIN1(TAUTOP,1.)
      
      CALL DELTA2STRIR(SRFALBIR,ASYIR,TAULAMIR,OMG0IR,
     & FUP,FDN,BPLANCK,TAUTOP,I, K1,IL)
     
C
C
 
!                if((I.eq.9).and.(K1.eq.16))then
!         print *,'OUT OF DELTA2STR', FUP(1), FDN(1)
!         pause
!       endif


                   DO J = 1, ND
                  FUPA(J)=FUPA(J)+TWGHT*FUP(J)
                  FDNA(J)=FDNA(J)+TWGHT*FDN(J)
!                  if( (NST==6) .and. (J<30) .and. (K1==1))then
!                  PRINT 3321,TWGHT,FUP(J),FDN(J),real(J)
!                  endif
                   ENDDO
!               pause
               sumtwght = sumtwght + TWGHT

               
  5     CONTINUE  ! ETHANE               
  4     CONTINUE ! METHANE           
  3     CONTINUE ! Co2
  2     CONTINUE ! H2O       
          write(5552,*) ! puts an extra space between each set of 16 coefficients

!         print *,'sum of TWGHT is.....',sumtwght, I
         
!         write(89,*) AV(I)/C, SUM_TPRIND   !summing tprind experiment
          write(89,3111) AV(I)/C, SUM_TPRIND
3111     FORMAT(1P5E14.5)

         DO 14 J=1,ND
             FDNIR(J)=FDNIR(J)+W(I)*FDNA(J)
             FUPIR(J)=FUPIR(J)+W(I)*FUPA(J)
!             PRINT 3321,W(I),FDNA(J),FUPA(J),real(J),real(I)
    
                 if (j.eq.1) then
              write(6969, *) FUPA(j), FDNA(j), I,NST
                  endif

3131          format(1p2e14.5,2(2x,i3))

!               if(NST==6)then
!                 print 3111,FUPIR(J),FUPA(J),real(J),real(I)
!               endif

             IF (J.eq.1) THEN
             write(90, *)AV(I)/C,ABS(FDNIR(J) - FUPIR(J)) ! Outgoing IR at top of atmosphere   
             ENDIF
c-jdh uncomment to get band-by-band fluxes            
c             IF (J.EQ.1) THEN
c  17            FORMAT(1PE13.5,3x,1PE22.10)
c               PRINT 17, AV(I), W(I)*FUPA(J)
c             END IF

c-jdh        FDNIR(J)=FDNIR(J)+VAC1*FDNA(J)
c-jdh        FUPIR(J)=FUPIR(J)+VAC1*FUPA(J)
             BPLANCK(J) = 0.    
  14    CONTINUE

!         pause

C        PRINT*, 'TRANSLAYER'
C        DO J=1,NLAYERS
C        PRINT 19, TRANSLAYER(J)
C         TAUTOTAL(I) = TAUTOTAL(I) - log(TRANSLAYER(J))
C        ENDDO
   1     CONTINUE                   !**END LOOP over frequency**

          write(6969,*)
19     FORMAT(/1X,1PE10.4)
100    FORMAT(1X,1P10E12.5) 
 !     PRINT*,'TAUTOTAL'
 !     PRINT 100, TAUTOTAL
 !     pause
      RETURN
      END
