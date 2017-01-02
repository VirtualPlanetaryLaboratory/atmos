
      SUBROUTINE PROFILE(TSTRAT,P,T,DZ,FH2O,FCO2V,BETA,JCOLD,
     &            IDRY,FLAGCONVEC)

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5)
      DIMENSION P(ND),T(ND),DZ(ND),FH2O(ND),BETA(ND),
     & FLAGCONVEC(ND),FCO2V(ND) ! YO(ND) - not used - EWS
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND),
     &  FH22
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     &  BETA1,BETA2,FVDRY,PDRY
      COMMON/CO2BLOK/betac1,betac2,PC0,TC0,VAPCL0,SUBCL0,DLVCDT
     & ,DLSCDT,CCL,CCS 


!      DATA YO/5.204400e-08,5.204400e-08,5.204400e-08,5.204400e-08,
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 ,
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 
!     &    5.204400e-08 , 5.204400e-08 , 5.204400e-08 , 5.583400e-08 , 
!     &    5.983700e-08 , 6.406100e-08 , 6.851200e-08 , 1.234600e-07 , 
!     &    2.171800e-07 , 3.733500e-07 , 6.278100e-07 , 1.033600e-06 , 
!     &    1.667600e-06 , 2.638700e-06 , 4.098000e-06 , 6.251500e-06 , 
!     &    9.373500e-06 , 1.382300e-05 , 2.006100e-05 , 2.866700e-05 ,
!     &    4.035300e-05 , 5.598300e-05 , 7.657500e-05 , 1.033000e-04 ,
!     &    1.374900e-04 , 1.805900e-04 , 2.341600e-04 , 2.997800e-04 , 
!     &    3.790400e-04 , 4.734400e-04 , 5.843300e-04 , 7.128400e-04 , 
!     &    8.598100e-04 , 1.025700e-03 , 1.210700e-03 , 1.414600e-03 , 
!     &    1.524300e-03/
C
C   THIS SUBROUTINE CALCULATES TROPOSPHERIC TEMPERATURES AND
C   HUMIDITIES BY CALLING CONVEC.
C
      BETA1 = 1.
      BETA2 = 1.
      ITROP = 1

c Settinbg all the temperatures at the stratospheric temperature
      do i=1, ND
       FLAGCONVEC(i)=0.
       T(i)=TSTRAT
       FCO2V(i)=FCO2
      enddo

C Calculating tropospheric  temperatures and water
C   SOLVE FROM THE GROUND UP
      PG = P(ND)
      CALL SATRAT(TG,PSAT)
      T(ND) = TG
      FH2O(ND) = RELHUM(PG) * PSAT/PG
       
      IF(PSAT.GT.PG) FH2O(ND) = POCEAN/PG
      IMCO2=0
    
      DO 2 J1=ND,JCOLD+1,-1
       T1 = T(J1)
       F1 = FH2O(J1)
        
       if (F1.ge.0.99) F1=.98 ! test
       P1 = P(J1)
       P2 = P(J1-1)
       DZP = DZ(J1)
       FC1 = FCO2
       
       CALL CONVEC(T1,T2,P1,P2,F1,FH2,FC1,FC2,DZP,ITROP,cflag,
     & IDRY,IMCO2)
       
       
       FLAGCONVEC(J1)=cflag
       BETA(J1) = BETA1
       BETA(J1-1) = BETA2
       JCOLD = J1
       IF(T2.LT.T(J1-1)) GOTO 6
       T(J1-1) = T2
       FCO2V(J1-1)= FC2
   2   FH2O(J1-1) = FH2
    
c Calculating stratospheric water contents

   6    IF(IMW.NE.2) GO TO 4

      DO 3 J=1,ND
       TJ = T(J)
       PJ = P(J)
       CALL SATRAT(TJ,PSAT)
       FH2 = RELHUM(PJ) * PSAT/PJ
   3   FH2O(J) = AMAx1(FH2,4.e-6) 
       
      RETURN
C
   4  CALL SATRAT(T(JCOLD),PSAT)
      FSAT= PSAT/P(JCOLD)
      DO 5 J=(JCOLD-1),1,-1
      CALL SATRAT(T(J),PSAT)
      FSATUR = (PSAT/P(J))
      FH2O(J) = FSATUR*RELHUM(P(J))
   5  FH2O(J) = AMIN1(FSAT,FH2O(J+1))
      

C-JFK 5/122011  Set stratospheric CO2 to the value at the CO2 cold trap
      ND1 = ND-1
      DO J=ND1,1,-1
      FCO2V(J) = AMIN1(FCO2V(J),FCO2V(J+1))
      END DO
	  
c to zero out water  
      IF (IMW.eq.5)THEN
   
      DO I = 1, ND
       FH2O(I)= 1.E-40
      ENDDO

      ENDIF

      RETURN
      END
C ----------------------------------------------------------------
