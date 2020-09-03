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

      RELHUMMW = Q2/0.98

      OMGA = 1  ! For normal Manabe-Wetherald
!--------------------------------------------Using Modified Cess et al. 1976 expression to calculation FH2O for all layers
       IF (IMW.eq.1)THEN
       CALL SATRAT(TG,PSAT)

       FSAT = PSAT/PG
       FP = 1.66e-02


      IF(FSAT.lt.FP)THEN
! If the saturation mixing ratio is less than saturation mixing ratio for present atmosphere (288 K)
! use normal Manabe-Wetherald.
       OMGA = 1.

      ELSEIF((FSAT.ge.FP).and.(FSAT.le.0.1))THEN ! If saturation mixing ratio is between
!present saturation mixing ratio and 0.1 use Kasting '86 sigma parameterization

       OMGA=1.-(FSAT-FP)/(0.1-FP)

      ELSEIF(FSAT.gt.0.1)THEN ! If the saturation mixing ratio is gt 0.1 then set RH = 80%
! throughout entire troposphere
        OMGA = 0.

      ENDIF
!S      print*, OMGA, '= this is OMEGA!'

       ENDIF   ! Only use this logic when Manabe-Wetherald is used. c-rr 10/23/2012
!---------------------------------------------------------------------------------------


      RELHUM = RSURF*(RELHUMMW**OMGA)

!      print 55,RELHUM, OMGA, P,PG,
!     &   FSAT

55    format(1x,'RELHUM=',1pe10.3,2x,'OMGA=',e10.3,2x,'P=',e10.3,2x,
     &  'PG=',e10.3, 3x, 'FSAT=', e10.3)


      IF (IMW.EQ.0) RELHUM = 1.
      IF (IMW.EQ.3) RELHUM = 0.5
      IF (IMW.EQ.4) RELHUM = 0.1
C-KK	added for low-O2 environments, to prevent water from
C-KK	zeroing itself out. 8% rel hum is present-day atmosphere
C-KK	at approximately 15 km (cold trap level)
      IF (RELHUM .LT. 0.08) RELHUM = 0.08

      RETURN
      END
