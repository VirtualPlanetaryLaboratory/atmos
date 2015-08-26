C     path:      $Source: /storm/rc1/cvsroot/rc/rrtm_lw/src/rtregcld.f,v $
C     author:    $Author: jdelamer $
C     revision:  $Revision: 3.2 $
C     created:   $Date: 2002/08/15 18:33:27 $
      SUBROUTINE RTREGCLD

C *** This program calculates the upward fluxes, downward fluxes, and
C     heating rates for an arbitrary cloudy atmosphere.  The input
C     to this program is the atmospheric profile, including cloud
C     properties,  and all needed Planck function information.  First  
C     order standard Gaussian quadrature is used for the angle 
C     integration.

C     Clouds are treated with random overlap scheme.

      PARAMETER (MG=16)
      PARAMETER (MXLAY=203)
      PARAMETER (MXANG = 4)
      PARAMETER (NBANDS = 16)
      PARAMETER (NTBL = 10000,TBLINT = 10000.0)

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
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)
      COMMON /PLNKDAT/   PLANKLAY(MXLAY,NBANDS),
     &                   PLANKLEV(0:MXLAY,NBANDS),PLANKBND(NBANDS)
      COMMON /PLANKG/    FRACS(MXLAY,MG)
      COMMON /TAUGCOM/   TAUG(MXLAY,MG)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /RTTBL/     BPADE,
     &                   TAUTBL(0:NTBL),TRANS(0:NTBL),TF(0:NTBL)

      COMMON /CVRRGC/    HVRRGC

      CHARACTER*16       HVRRGC

      DIMENSION ATRANS(MXLAY,MXANG),BBUGAS(MXLAY,MXANG)
      DIMENSION ATOT(MXLAY,MXANG),ODCLD(MXLAY,NBANDS,MXANG)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY),BBUTOT(MXLAY,MXANG)
      DIMENSION DRAD(0:MXLAY-1,MXANG),URAD(0:MXLAY,MXANG)
      DIMENSION SECANG(MXANG),ANGWEIGH(MXANG),RAD(MXANG)
      DIMENSION SECREG(MXANG,MXANG),WTREG(MXANG,MXANG)
      DIMENSION EFCLFRAC(MXLAY,NBANDS,MXANG)
      DIMENSION ABSCLD(MXLAY,NBANDS,MXANG)
      DIMENSION IPAT(16,0:2)

C Dimensions for cloud 
      DIMENSION ICLDLYR(MXLAY)

      DATA IPAT /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1,
     &           1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5,
     &           1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

C *** When standard first-order Gaussian quadrature is chosen as
C     the method to approximate the integral over angles that yields
C     flux from radiances, then SECREG(I,J) is the secant of the Ith  
C     (out of a total of J angles) and WTREG(I,J) is the corresponding
C     weight.
      DATA SECREG(1,1) / 1.5/
      DATA SECREG(2,2) / 2.81649655/, SECREG(1,2) / 1.18350343/
      DATA SECREG(3,3) / 4.70941630/, SECREG(2,3) / 1.69338507/
      DATA SECREG(1,3) / 1.09719858/
      DATA SECREG(4,4) / 7.15513024/, SECREG(3,4) / 2.40148179/
      DATA SECREG(2,4) / 1.38282560/, SECREG(1,4) / 1.06056257/
      DATA WTREG(1,1) / 0.50000000/
      DATA WTREG(2,2) /0.1819586183/, WTREG(1,2) /0.3180413817/
      DATA WTREG(3,3) /0.0698269799/, WTREG(2,3) /0.2292411064/
      DATA WTREG(1,3) /0.2009319137/
      DATA WTREG(4,4) /0.0311809710/, WTREG(3,4) /0.1298475476/
      DATA WTREG(2,4) /0.2034645680/, WTREG(1,4) /0.1355069134/
      DATA REC_6 /0.166667/

      HVRRGC = '$Revision: 3.2 $'

      RADSUM = 0.
      NUMANG = ABS(NUMANGS)
C *** Load angle data in arrays depending on angular quadrature scheme.
      DO 100 IANG = 1, NUMANG
         SECANG(IANG) = SECREG(IANG,NUMANG)
         ANGWEIGH(IANG) = WTREG(IANG,NUMANG)
 100  CONTINUE

      
      TOTUFLUX(0) = 0.0
      TOTDFLUX(0) = 0.0
      DO 120 IANG = 1, NUMANG
         URAD(0,IANG) = 0.0
         DRAD(0,IANG) = 0.0
 120  CONTINUE
      DO 200 LAY = 1, NLAYERS
         TOTUFLUX(LAY) = 0.0
         TOTDFLUX(LAY) = 0.0
         DO 180 IANG = 1, NUMANG
            URAD(LAY,IANG) = 0.
            DRAD(LAY,IANG) = 0.
            DO 150 IB = 1, NCBANDS
               IF (CLDFRAC(LAY) .GE. 1.E-6) THEN
                  ODCLD(LAY,IB,IANG) = SECANG(IANG) * TAUCLOUD(LAY,IB)
                  TRANSCLD = EXP(-ODCLD(LAY,IB,IANG))
                  ABSCLD(LAY,IB,IANG) = 1. - TRANSCLD
                  EFCLFRAC(LAY,IB,IANG) = ABSCLD(LAY,IB,IANG) * 
     &                 CLDFRAC(LAY)
                  ICLDLYR(LAY) = 1
               ELSE
                  ODCLD(LAY,IB,IANG) = 0.0
                  ABSCLD(LAY,IB,IANG) = 0.0
                  EFCLFRAC(LAY,IB,IANG) = 0.0
                  ICLDLYR(LAY) = 0
               ENDIF
 150        CONTINUE
 180     CONTINUE
 200  CONTINUE

C *** Loop over frequency bands.
      DO 6000 IBAND = ISTART, IEND
         IF (NCBANDS .EQ. 1) THEN
            IB = IPAT(IBAND,0)
         ELSEIF (NCBANDS .EQ.  5) THEN
            IB = IPAT(IBAND,1)
         ELSEIF (NCBANDS .EQ. 16) THEN
            IB = IPAT(IBAND,2)
         ENDIF
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
         ELSE
            CALL TAUGB16
         ENDIF

         
C ***    Loop over g-channels.
         IG = 1
 1000    CONTINUE
         
C ***    Loop over each angle for which the radiance is to be computed.
         DO 3000 IANG = 1, NUMANG
C ***       Radiative transfer starts here.
            RADLD = 0.

C ***       Downward radiative transfer.  
            DO 2500 LEV = NLAYERS, 1, -1
               PLFRAC = FRACS(LEV,IG)
               BLAY = PLANKLAY(LEV,IBAND)
               DPLANKUP = PLANKLEV(LEV,IBAND) - BLAY
               DPLANKDN = PLANKLEV(LEV-1,IBAND) - BLAY
               ODEPTH = SECANG(IANG) * TAUG(LEV,IG)
               IF (ODEPTH .LT. 0.0) ODEPTH = 0.0
               IF (ICLDLYR(LEV).EQ.1) THEN
                  ODTOT = ODEPTH + ODCLD(LEV,IB,IANG)
                  IF (ODTOT .LT. 0.06) THEN
                     ATRANS(LEV,IANG) = ODEPTH - 
     &                    0.5*ODEPTH*ODEPTH
                     ODEPTH_REC = REC_6*ODEPTH
                     GASSRC = PLFRAC*(BLAY+DPLANKDN*ODEPTH_REC)
     &                    *ATRANS(LEV,IANG)

                     ATOT(LEV,IANG) =  ODTOT - 0.5*ODTOT*ODTOT
                     ODTOT_REC = REC_6*ODTOT
                     BBDTOT =  PLFRAC * (BLAY+DPLANKDN*ODTOT_REC)
                     RADLD = RADLD - RADLD * (ATRANS(LEV,IANG) +
     &                    EFCLFRAC(LEV,IB,IANG) * 
     &                    (1. - ATRANS(LEV,IANG))) +
     &                    GASSRC + CLDFRAC(LEV) * 
     &                    (BBDTOT * ATOT(LEV,IANG) - GASSRC)
                     DRAD(LEV-1,IANG) = DRAD(LEV-1,IANG) + RADLD
                  
                     BBUGAS(LEV,IANG) =  PLFRAC *
     &                    (BLAY+DPLANKUP*ODEPTH_REC)
                     BBUTOT(LEV,IANG) =  PLFRAC * 
     &                    (BLAY+DPLANKUP*ODTOT_REC)
                  ELSEIF (ODEPTH .LE. 0.06) THEN
                     ATRANS(LEV,IANG) = ODEPTH -
     &                    0.5*ODEPTH*ODEPTH
                     ODEPTH_REC = REC_6*ODEPTH
                     GASSRC = PLFRAC*(BLAY+DPLANKDN*ODEPTH_REC)
     &                    *ATRANS(LEV,IANG)

                     ODTOT = ODEPTH + ODCLD(LEV,IB,IANG)
                     TBLIND = ODTOT/(BPADE+ODTOT)
                     ITTOT = TBLINT*TBLIND + 0.5
                     TFACTOT = TF(ITTOT)
                     BBDTOT = PLFRAC * (BLAY + TFACTOT*DPLANKDN)
                     ATOT(LEV,IANG) = 1. - TRANS(ITTOT)

                     RADLD = RADLD - RADLD * (ATRANS(LEV,IANG) +
     &                    EFCLFRAC(LEV,IB,IANG) * 
     &                    (1. - ATRANS(LEV,IANG))) +
     &                    GASSRC + CLDFRAC(LEV) * 
     &                    (BBDTOT * ATOT(LEV,IANG) - GASSRC)
                     DRAD(LEV-1,IANG) = DRAD(LEV-1,IANG) + RADLD

                     BBUGAS(LEV,IANG) = PLFRAC * 
     &                    (BLAY + DPLANKUP*ODEPTH_REC)
                     BBUTOT(LEV,IANG) = PLFRAC * 
     &                    (BLAY + TFACTOT * DPLANKUP)
                  ELSE
                     TBLIND = ODEPTH/(BPADE+ODEPTH)
                     ITGAS = TBLINT*TBLIND+0.5
                     ODEPTH = TAUTBL(ITGAS)
                     ATRANS(LEV,IANG) = 1. - TRANS(ITGAS)
                     TFACGAS = TF(ITGAS)
                     GASSRC = ATRANS(LEV,IANG) * PLFRAC * 
     &                    (BLAY + TFACGAS*DPLANKDN)

                     ODTOT = ODEPTH + ODCLD(LEV,IB,IANG)
                     TBLIND = ODTOT/(BPADE+ODTOT)
                     ITTOT = TBLINT*TBLIND + 0.5
                     TFACTOT = TF(ITTOT)
                     BBDTOT = PLFRAC * (BLAY + TFACTOT*DPLANKDN)
                     ATOT(LEV,IANG) = 1. - TRANS(ITTOT)

                     RADLD = RADLD - RADLD * (ATRANS(LEV,IANG) +
     &                    EFCLFRAC(LEV,IB,IANG) * 
     &                    (1. - ATRANS(LEV,IANG))) +
     &                    GASSRC + CLDFRAC(LEV) * 
     &                    (BBDTOT * ATOT(LEV,IANG) - GASSRC)
                     DRAD(LEV-1,IANG) = DRAD(LEV-1,IANG) + RADLD
                     BBUGAS(LEV,IANG) = PLFRAC * 
     &                    (BLAY + TFACGAS * DPLANKUP)
                     BBUTOT(LEV,IANG) = PLFRAC * 
     &                    (BLAY + TFACTOT * DPLANKUP)
                  ENDIF
               ELSE
                  IF (ODEPTH .LE. 0.06) THEN
                     ATRANS(LEV,IANG) = ODEPTH-0.5*ODEPTH*ODEPTH
                     ODEPTH = REC_6*ODEPTH
                     BBD = PLFRAC*(BLAY+DPLANKDN*ODEPTH)
                     BBUGAS(LEV,IANG) = PLFRAC*
     &                    (BLAY+DPLANKUP*ODEPTH)
                  ELSE
                     TBLIND = ODEPTH/(BPADE+ODEPTH)
                     ITR = TBLINT*TBLIND+0.5
                     TRANSC = TRANS(ITR)
                     ATRANS(LEV,IANG) = 1.-TRANSC
                     TAUSFAC = TF(ITR)
                     BBD = PLFRAC*(BLAY+TAUSFAC*DPLANKDN)
                     BBUGAS(LEV,IANG) = PLFRAC * 
     &                    (BLAY + TAUSFAC * DPLANKUP)
                  ENDIF   
                  RADLD = RADLD + (BBD-RADLD)*ATRANS(LEV,IANG)
                  DRAD(LEV-1,IANG) = DRAD(LEV-1,IANG) + RADLD
               ENDIF
 2500       CONTINUE
            RAD(IANG) = RADLD
            RADSUM = RADSUM + ANGWEIGH(IANG) * RADLD
 3000    CONTINUE 

         RAD0 = FRACS(1,IG) * PLANKBND(IBAND)
         REFLECT = 1. - SEMISS(IBAND)
         DO 4000 IANG = 1, NUMANG
C           Add in reflection of surface downward radiance.
            IF (IREFLECT .EQ. 1) THEN
C              Specular reflection.
               RADLU = RAD0 + REFLECT * RAD(IANG)
            ELSE
C              Lambertian reflection.
               RADLU = RAD0 + 2. * REFLECT * RADSUM
            ENDIF
            URAD(0,IANG) = URAD(0,IANG) + RADLU
C ***       Upward radiative transfer.
            DO 2600 LEV = 1, NLAYERS
               IF (ICLDLYR(LEV) .EQ. 1) THEN
                  GASSRC = BBUGAS(LEV,IANG) * ATRANS(LEV,IANG)
                  RADLU = RADLU - RADLU * (ATRANS(LEV,IANG) +
     &                 EFCLFRAC(LEV,IB,IANG) * 
     &                 (1. - ATRANS(LEV,IANG))) +
     &                 GASSRC + CLDFRAC(LEV) * 
     &                 (BBUTOT(LEV,IANG) * ATOT(LEV,IANG) - GASSRC)
                  URAD(LEV,IANG) = URAD(LEV,IANG) + RADLU
               ELSE
                  RADLU = RADLU + (BBUGAS(LEV,IANG)-RADLU)
     &                 *ATRANS(LEV,IANG)
                  URAD(LEV,IANG) = URAD(LEV,IANG) + RADLU
               ENDIF
 2600       CONTINUE
 4000    CONTINUE 
         RADSUM = 0.

         IG = IG + 1
         IF (IG .LE. 16) GO TO 1000

C ***    Calculate upward, downward, and net flux.
         DO 5000 LEV = NLAYERS, 0, -1
            UFLUX(LEV) = 0.0
            DFLUX(LEV) = 0.0
            DO 4500 IANG = 1, NUMANG
               UFLUX(LEV) = UFLUX(LEV) + URAD(LEV,IANG)*ANGWEIGH(IANG)
               DFLUX(LEV) = DFLUX(LEV) + DRAD(LEV,IANG)*ANGWEIGH(IANG)
               URAD(LEV,IANG) = 0.0
               DRAD(LEV,IANG) = 0.0
 4500       CONTINUE
            TOTUFLUX(LEV) = TOTUFLUX(LEV) + UFLUX(LEV) * DELWAVE(IBAND)
            TOTDFLUX(LEV) = TOTDFLUX(LEV) + DFLUX(LEV) * DELWAVE(IBAND)
 5000    CONTINUE
 6000 CONTINUE 

      DO 7000 LEV = 0, NLAYERS
         TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
         TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
         FNET(LEV) = TOTUFLUX(LEV) - TOTDFLUX(LEV)
 7000 CONTINUE

C *** Calculate Heating Rates.
      HTR(NLAYERS) = 0.0
      DO 8000 LEV  = NLAYERS-1, 0, -1
         LAY = LEV + 1
         HTR(LEV)=HEATFAC*(FNET(LEV)-FNET(LAY))/(PZ(LEV)-PZ(LAY)) 
 8000 CONTINUE

      RETURN
      END   


