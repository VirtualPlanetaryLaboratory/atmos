C     path:      $Source: /storm/rc1/cvsroot/rc/rrtm_lw/src/rtr.f,v $
C     author:    $Author: jdelamer $
C     revision:  $Revision: 3.1 $
C     created:   $Date: 2002/08/15 18:33:27 $
      SUBROUTINE RTR

C *** This program calculates the upward fluxes, downward fluxes,
C     and heating rates for an arbitrary atmosphere.  The input to
C     this program is the atmospheric profile and all Planck function
C     information.  Only one angle is used from standard 
C     Gaussian quadrature.

      PARAMETER (MG=16)
      PARAMETER (MXLAY=203)
      PARAMETER (MXANG = 4)
      PARAMETER (NBANDS = 16)
      PARAMETER (NTBL = 10000,TBLINT=10000.0)

      IMPLICIT DOUBLE PRECISION (V)                                     

      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND
      COMMON /PROFI/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SURFACE/   TBOUND,IREFLECT,SEMISS(NBANDS)
      COMMON /PLNKDAT/   PLANKLAY(MXLAY,NBANDS),
     &                   PLANKLEV(0:MXLAY,NBANDS),PLANKBND(NBANDS)
      COMMON /PLANKG/    FRACS(MXLAY,MG)
      COMMON /TAUGCOM/   TAUG(MXLAY,MG)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /RTTBL/     BPADE,
     &                   TAUTBL(0:NTBL),TRANS(0:NTBL),TF(0:NTBL)

      COMMON /CVRRTR/    HVRRTR

      CHARACTER*16       HVRRTR

      DIMENSION BBU1(MXLAY)
      DIMENSION ATRANS1(MXLAY)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY)
      DIMENSION DRAD1(0:MXLAY),URAD1(0:MXLAY)

      HVRRTR = '$Revision: 3.1 $'

      DATA SECDIFF /1.66/
      DATA WTDIFF /0.5/
      DATA REC_6 /0.166667/

      DO 200 LAY = 0, NLAYERS
         URAD1(LAY) = 0.0
         DRAD1(LAY) = 0.0
         TOTUFLUX(LAY) = 0.0
         TOTDFLUX(LAY) = 0.0
 200  CONTINUE

C *** Loop over frequency bands.
      DO 6000 IBAND = ISTART, IEND
        IF (IBAND .EQ. 1) THEN
            CALL TAUGB1
         ELSEIF (IBAND .EQ. 2) THEN
            CALL TAUGB2
         ELSEIF (IBAND .EQ. 3) THEN
            CALL TAUGB3
         ELSEIF (IBAND .EQ. 4) THEN
            CALL TAUGB4
         ELSEIF (IBAND .EQ. 5) THEN
            CALL TAUGB5
         ELSEIF (IBAND .EQ. 6) THEN
            CALL TAUGB6
         ELSEIF (IBAND .EQ. 7) THEN
            CALL TAUGB7
         ELSEIF (IBAND .EQ. 8) THEN
            CALL TAUGB8
         ELSEIF (IBAND .EQ. 9) THEN
            CALL TAUGB9
         ELSEIF (IBAND .EQ. 10) THEN
            CALL TAUGB10
         ELSEIF (IBAND .EQ. 11) THEN
            CALL TAUGB11
         ELSEIF (IBAND .EQ. 12) THEN
            CALL TAUGB12
         ELSEIF (IBAND .EQ. 13) THEN
            CALL TAUGB13
         ELSEIF (IBAND .EQ. 14) THEN
            CALL TAUGB14
         ELSEIF (IBAND .EQ. 15) THEN
            CALL TAUGB15
         ELSEIF (IBAND .EQ. 16) THEN
            CALL TAUGB16
         ENDIF
C ***    Loop over g-channels.

         IG = 1
 1000    CONTINUE
C ***    Radiative transfer starts here.

         RADLD1 = 0.
C ***    Downward radiative transfer.  Due to the simple form taken by
C        certain equations when the optical depth is small,
C        these conditions are tested for.  
         DO 2500 LEV = NLAYERS, 1, -1
            BLAY = PLANKLAY(LEV,IBAND)
            PLFRAC = FRACS(LEV,IG)
            DPLANKUP = PLANKLEV(LEV,IBAND) - BLAY
            DPLANKDN = PLANKLEV(LEV-1,IBAND) - BLAY
            ODEPTH1 = SECDIFF*TAUG(LEV,IG)
            IF (ODEPTH1 .LE. 0.06) THEN
               IF (ODEPTH1 .LT. 0.0) ODEPTH1 = 0.0
               ATRANS1(LEV) = ODEPTH1-0.5*ODEPTH1*ODEPTH1
               ODEPTH1 = REC_6*ODEPTH1
               BBD1 =  PLFRAC*(BLAY+DPLANKDN*ODEPTH1)
               RADLD1 = RADLD1 + (BBD1-RADLD1)*ATRANS1(LEV)
               BBU1(LEV) =  PLFRAC*(BLAY+DPLANKUP*ODEPTH1)
            ELSE
               TBLIND = ODEPTH1/(BPADE+ODEPTH1)
               ITR1 = TBLINT*TBLIND+0.5
               TRANS1 = TRANS(ITR1)
               ATRANS1(LEV) = 1. - TRANS1
               TAUSFAC1 = TF(ITR1)
               BBD1 = PLFRAC * (BLAY + TAUSFAC1 * DPLANKDN)
               RADLD1 = RADLD1 + (BBD1-RADLD1)*ATRANS1(LEV)
               BBU1(LEV) = PLFRAC * (BLAY + TAUSFAC1 * DPLANKUP)
            ENDIF   
            DRAD1(LEV-1) = DRAD1(LEV-1) + RADLD1
 2500    CONTINUE

C ***    Upward radiative transfer.
         RAD0 = FRACS(1,IG) * PLANKBND(IBAND)
C        Add in reflection of surface downward radiance.
         REFLECT = 1. - SEMISS(IBAND)
         IF (IREFLECT .EQ. 1) THEN
C           Specular reflection.
            RADLU1 = RAD0 + REFLECT * RADLD1
         ELSE
C           Lambertian reflection.
            RAD = 2. * (RADLD1*WTDIFF)
            RADLU1 = RAD0 + REFLECT * RAD
         ENDIF
         URAD1(0) = URAD1(0) + RADLU1
         DO 2600 LEV = 1, NLAYERS
            RADLU1 = RADLU1 + (BBU1(LEV)-RADLU1)*ATRANS1(LEV)
            URAD1(LEV) = URAD1(LEV) + RADLU1
 2600    CONTINUE
c 4000    CONTINUE !EWS - label not used
         IG = IG + 1
         IF (IG .LE. NG(IBAND)) GO TO 1000
         
C ***    Process longwave output from band.
C ***    Calculate upward, downward, and net flux.
         DO 5000 LEV = NLAYERS, 0, -1
            UFLUX(LEV) = URAD1(LEV)*WTDIFF
            DFLUX(LEV) = DRAD1(LEV)*WTDIFF
            URAD1(LEV) = 0.0
            DRAD1(LEV) = 0.0
            TOTUFLUX(LEV) = TOTUFLUX(LEV) + UFLUX(LEV) * DELWAVE(IBAND)
            TOTDFLUX(LEV) = TOTDFLUX(LEV) + DFLUX(LEV) * DELWAVE(IBAND)
 5000    CONTINUE
 6000 CONTINUE
      TOTUFLUX(0) = TOTUFLUX(0) * FLUXFAC
      TOTDFLUX(0) = TOTDFLUX(0) * FLUXFAC
      FNET(0) = TOTUFLUX(0) - TOTDFLUX(0)
      DO 7000 LEV = 1, NLAYERS
         TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
         TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
         FNET(LEV) = TOTUFLUX(LEV) - TOTDFLUX(LEV)
         L = LEV - 1
C        Calculate Heating Rates.
         HTR(L)=HEATFAC*(FNET(L)-FNET(LEV))/(PZ(L)-PZ(LEV)) 
 7000 CONTINUE
      HTR(NLAYERS) = 0.0

c 9000 CONTINUE !EWS - label not used

      RETURN

      END   

