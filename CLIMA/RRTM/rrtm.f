C     path:      $Source: /storm/rc1/cvsroot/rc/rrtm_lw/src/rrtm.f,v $
C     author:    $Author: jdelamer $
C     revision:  $Revision: 3.2 $
C     created:   $Date: 2002/08/15 18:33:27 $
****************************************************************************
*                                                                          *
*                               RRTM                                       *
*                                                                          *
*                                                                          *
*                                                                          *
*                   A RAPID RADIATIVE TRANSFER MODEL                       *
*                       FOR THE LONGWAVE REGION                            * 
*                                                                          *
*                                                                          *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                  *
*                        131 HARTWELL AVENUE                               *
*                        LEXINGTON, MA 02421                               *
*                                                                          *
*                                                                          *
*                         ELI J. MLAWER                                    *
*                         JENNIFER S. DELAMERE                             *
*                         STEVEN J. TAUBMAN~                               *
*                         SHEPARD A. CLOUGH                                *
*                                                                          *
*                                                                          *
*                         ~currently at GFDL                               *
*                                                                          *
*                                                                          *
*                                                                          *
*                       email:  mlawer@aer.com                             *
*                       email:  jdelamer@aer.com                           *
*                                                                          *
*        The authors wish to acknowledge the contributions of the          *
*        following people:  Karen Cady-Pereira, Patrick D. Brown,          *
*        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom. *
*                                                                          *
****************************************************************************

C THIS PROGRAM HAS BEEN MODIFIED TO RUN AS A SUBROUTINE OF THE CLIMATE MODEL
C ALL CHANGES ARE MARKED WITH THE WORD ojo (MINDS "EYE" IN SPANISH)
C CHANGES WERE ORIGINALLY MADE BY KARA KRELOVE FOR RRTM V2.3 
C THIS PROGRAM WAS MODIFIED BY ANTIGONA SEGURA (AUGUST/2003)

c ojo PROGRAM was changed for SUBROUTINE
       SUBROUTINE RRTM
                    
C *** This program is the driver for RRTM, the AER rapid model.  
C     For each atmosphere the user wishes to analyze, this routine
C     a) calls READPROF to read in the atmospheric profile
C     b) calls SETCOEF to calculate various quantities needed for 
C        the radiative transfer algorithm
C     c) calls RTR or RTREG (depending on angular quadrature
C         method) to do the radiative transfer calculation for clear sky 
C         calculstions OR calls RTRCLD or RTREGCLD for calculations 
C         with cloudy skies and a random cloud overlap scheme OR
C         calls RTRCLDMR or RTREGCLDMR for calculations with cloud skies
C         and a maximum/random cloud overlap scheme
C     d) writes out the upward, downward, and net flux for each
C        level and the heating rate for each layer

      PARAMETER (MXLAY=203)
      PARAMETER (MG = 16)
      PARAMETER (NBANDS = 16)
      PARAMETER (NTBL = 10000,TBLINT=10000.0)

      CHARACTER :: DIRINOUT*2,DIRDATA*4
c-ojo Common file to direct the output file to the right subdirectory
      COMMON/DIR/DIRINOUT,DIRDATA

      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/   ONEMINUS
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /IFIL/      IRD,IPR,IPU,IDUM(15)
      COMMON /PROFI/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /RTTBL/     BPADE,
     &                   TAUTBL(0:NTBL),TRANS(0:NTBL), TF(0:NTBL)

      COMMON /CVRRTM/    HVRRTM
      COMMON /CVRREG/    HVRREG
      COMMON /CVRRTR/    HVRRTR
      COMMON /CVRATM/    HVRATM
      COMMON /CVRSET/    HVRSET
      COMMON /CVRTAU/    HVRTAU
      COMMON /CVRRGC/    HVRRGC
      COMMON /CVRRTC/    HVRRTC
      COMMON /CVRCLD/    HVRCLD
      COMMON /CVRUTL/    HVRUTL
      COMMON /CVREXT/    HVREXT
      COMMON /CVRRTX/    HVRRTX
      COMMON /CVRRGX/    HVRRGX

      COMMON /HVRSN1/    HVRKG1
      COMMON /HVRSN2/    HVRKG2
      COMMON /HVRSN3/    HVRKG3
      COMMON /HVRSN4/    HVRKG4
      COMMON /HVRSN5/    HVRKG5
      COMMON /HVRSN6/    HVRKG6
      COMMON /HVRSN7/    HVRKG7
      COMMON /HVRSN8/    HVRKG8
      COMMON /HVRSN9/    HVRKG9
      COMMON /HVRSN10/   HVRKG10
      COMMON /HVRSN11/   HVRKG11
      COMMON /HVRSN12/   HVRKG12
      COMMON /HVRSN13/   HVRKG13
      COMMON /HVRSN14/   HVRKG14
      COMMON /HVRSN15/   HVRKG15
      COMMON /HVRSN16/   HVRKG16

C ojo COMMON BLOCK added
      COMMON /ALTI/  ALTZ(0:MXLAY)


      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBRODL(MXLAY),
     &                  COLMOL(MXLAY),NMOL


      CHARACTER*16 HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVRRGC,HVRRTC,HVRCLD,HVRUTL,HVREXT,
     *            HVRRTX,HVRRGX

      CHARACTER*16 HVRKG1,HVRKG2,HVRKG3,HVRKG4,HVRKG5
      CHARACTER*16 HVRKG6,HVRKG7,HVRKG8,HVRKG9,HVRKG10,HVRKG11
      CHARACTER*16 HVRKG12,HVRKG13,HVRKG14,HVRKG15,HVRKG16
      CHARACTER PAGE

C      DATA WAVENUM1(1) /10./, WAVENUM2(1) /350./, DELWAVE(1) /340./
C      DATA WAVENUM1(2) /350./, WAVENUM2(2) /500./, DELWAVE(2) /150./
C      DATA WAVENUM1(3) /500./, WAVENUM2(3) /630./, DELWAVE(3) /130./
C      DATA WAVENUM1(4) /630./, WAVENUM2(4) /700./, DELWAVE(4) /70./
C      DATA WAVENUM1(5) /700./, WAVENUM2(5) /820./, DELWAVE(5) /120./
C      DATA WAVENUM1(6) /820./, WAVENUM2(6) /980./, DELWAVE(6) /160./
C      DATA WAVENUM1(7) /980./, WAVENUM2(7) /1080./, DELWAVE(7) /100./
C      DATA WAVENUM1(8) /1080./, WAVENUM2(8) /1180./, DELWAVE(8) /100./
C      DATA WAVENUM1(9) /1180./, WAVENUM2(9) /1390./, DELWAVE(9) /210./
C      DATA WAVENUM1(10) /1390./,WAVENUM2(10) /1480./,DELWAVE(10) /90./
C      DATA WAVENUM1(11) /1480./,WAVENUM2(11) /1800./,DELWAVE(11) /320./
C      DATA WAVENUM1(12) /1800./,WAVENUM2(12) /2080./,DELWAVE(12) /280./
C      DATA WAVENUM1(13) /2080./,WAVENUM2(13) /2250./,DELWAVE(13) /170./
C      DATA WAVENUM1(14) /2250./,WAVENUM2(14) /2380./,DELWAVE(14) /130./
C      DATA WAVENUM1(15) /2380./,WAVENUM2(15) /2600./,DELWAVE(15) /220./
C      DATA WAVENUM1(16) /2600./,WAVENUM2(16) /3250./,DELWAVE(16) /650./

C      DATA NG  /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
C      DATA NSPA /1, 1,10, 9, 9, 1, 9, 1,11, 1, 1, 9, 9, 1, 9, 9/
C      DATA NSPB /1, 1, 5, 6, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
C      DATA HEATFAC /8.4391/

      ONEMINUS = 1. - 1.E-6
      FLUXFAC = PI * 2.D4  

      IWR = 10
      PAGE = CHAR(12)

      HVRRTM = '$Revision: 3.2 $'

c uncomment to get band-by-band output
c     IOUT = 99


C  Compute lookup tables for transmittance, tau transition function,
C  and clear sky tau (for the cloudy sky radiative transfer).  Tau is 
C  computed as a function of the tau transition function, transmittance 
C  is calculated as a function of tau, and the tau transition function 
C  is calculated using the linear in tau formulation at values of tau 
C  above 0.01.  TF is approximated as tau/6 for tau < 0.01.  All tables 
C  are computed at intervals of 0.001.  The inverse of the constant used
C  in the Pade approximation to the tau transition function is set to b.

      TAUTBL(0) = 0.0
      TAUTBL(NTBL) = 1.E10
      TRANS(0) = 1.0
      TRANS(NTBL) = 0.0
      TF(0) = 0.0
      TF(NTBL) = 1.0
      PADE  = 0.278
      BPADE = 1.0/PADE
      DO 500 ITR = 1,NTBL-1
         TFN = ITR/FLOAT(NTBL)
         TAUTBL(ITR) = BPADE*TFN/(1.-TFN)
         TRANS(ITR) = EXP(-TAUTBL(ITR))
         IF (TAUTBL(ITR) .LT. 0.06) THEN
            TF(ITR) = TAUTBL(ITR)/6.
         ELSE
            TF(ITR) = 1.-
     &           2.*((1./TAUTBL(ITR))-(TRANS(ITR)/(1.-TRANS(ITR))))
         ENDIF
 500     CONTINUE

C     Open the INPUT set of atmospheres
c ojo INPUT FILE commented
c      IRD = 9
c      OPEN (IRD,FILE='INPUT_RRTM',FORM='FORMATTED')

c Multiple atmosphere option implemented      
c ojo NUMATMOS must be set = 1
      NUMATMOS = 1
      DO 4000 IATMOS = 1, NUMATMOS
C ***    Input atmospheric profile from INPUT_RRTM.

         CALL READPROF

         ISTART = 1
         IEND = 16
         IFLAG = IOUT
     
 1000    CONTINUE
         IF (IFLAG .GT. 0 .AND. IFLAG .LE. 40) THEN
            ISTART = IFLAG
            IEND = IFLAG
         ENDIF


C ***    Calculate information needed by the radiative transfer routine
C        that is specific to this atmosphere, especially some of the 
C        coefficients and indices needed to compute the optical depths
C        by interpolating data from stored reference atmospheres. 
         ICLDATM = 0
         IF (ICLD .GE. 1) CALL CLDPROP(ICLDATM)

         CALL SETCOEF

c      DO J=0,NLAYERS
c   17   FORMAT (I2,3x,F7.3,3x,F8.3,3x,PE11.5,3x,PE11.5,3x,PE11.5)  
c        PRINT 17, J,TZ(J),PZ(J),WKL(1,J+1),WKL(2,J+1),WKL(6,J+1)
c        PRINT 17, J,TAVEL(J),PAVEL(J),WBRODL(J)
c      END DO
         
C ***    Call the radiative transfer routine.
         IF (NUMANGS .EQ. 0 .AND. ICLDATM .EQ. 0) THEN
            CALL RTR
         ELSEIF (NUMANGS .EQ. 0 .AND. ICLDATM .EQ. 1) THEN
            IF (ICLD .EQ. 2) THEN
               CALL RTRCLDMR 
            ELSE 
               CALL RTRCLD
            ENDIF
         ELSEIF (ICLDATM .EQ. 1) THEN
            IF (ICLD .EQ. 2) THEN
               CALL RTREGCLDMR 
            ELSE 
               CALL RTREGCLD
            ENDIF
         ELSE
            CALL RTREG
         ENDIF
         IF (IOUT .LT. 0) GO TO 4000

C ***    Process output for this atmosphere.
         OPEN (IWR,FILE=DIRINOUT//'/OUTPUT_RRTM',FORM='FORMATTED')
         WRITE(IWR,9899)WAVENUM1(ISTART),WAVENUM2(IEND),IATMOS
         WRITE(IWR,9900)
         WRITE(IWR,9901)
C

         PRINT *, "NLAYERS = ", NLAYERS

         DO 3000 I = NLAYERS, 0, -1
            IF (PZ(I) .LT. 1.E-2) THEN
               WRITE(IWR,9952) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 1.E-1) THEN
               WRITE(IWR,9953) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 1.) THEN
               WRITE(IWR,9954) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 10.) THEN
               WRITE(IWR,9955) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 100.) THEN
               WRITE(IWR,9956) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 1000.) THEN
               WRITE(IWR,9957) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSE
               WRITE(IWR,9958) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ENDIF
 3000    CONTINUE
         WRITE(IWR,9903)PAGE
         
         IF (IOUT .LE. 40 .OR. IFLAG .EQ. 16) GO TO 3500
         IF (IFLAG .EQ. 99) THEN
            IFLAG = 1
         ELSEIF (IOUT .EQ. 99) THEN
            IFLAG = IFLAG + 1
         ENDIF
         GO TO 1000

 3500    CONTINUE
C
C ***    Output module version numbers
C
         WRITE(IWR,9910) HVRRTM,HVRATM,HVRRTR,HVRRTC,HVRREG,HVRRGC,
     *     HVRRTX,HVRRGX,HVRSET,HVRCLD,HVRUTL,HVRTAU,HVRKG1,HVRKG2,
     *     HVRKG3,HVRKG4,HVRKG5,HVRKG6,HVRKG7,HVRKG8,HVRKG9,HVRKG10,
     *     HVRKG11,HVRKG12,HVRKG13,HVRKG14,HVRKG15,HVRKG16

 4000 CONTINUE

c ojo CLOSE output file commented
c         CLOSE(IRD)
         CLOSE(IWR)

 9952 FORMAT(1X,I3,9X,F7.6,3X,F8.4,6X,F8.4,6X,F9.4,10X,F9.5)
 9953 FORMAT(1X,I3,9X,F6.5,4X,F8.4,6X,F8.4,6X,F9.4,10X,F9.5)
 9954 FORMAT(1X,I3,9X,F5.4,5X,F8.4,6X,F8.4,6X,F9.4,10X,F9.5)
 9955 FORMAT(1X,I3,8X,F5.3,6X,F8.4,6X,F8.4,6X,F9.4,10X,F9.5)
 9956 FORMAT(1X,I3,7X,F5.2,7X,F8.4,6X,F8.4,6X,F9.4,10X,F9.5)
 9957 FORMAT(1X,I3,6X,F5.1,8X,F8.4,6X,F8.4,6X,F9.4,10X,F9.5)
 9958 FORMAT(1X,I3,5X,F5.0,9X,F8.4,6X,F8.4,6X,F9.4,10X,F9.5)
 9899 FORMAT(1X,'Wavenumbers: ',F6.1,' - ',F6.1,' cm-1, ATM ',i6)
 9900 FORMAT(1X,'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET
     &FLUX       HEATING RATE')
 9901 FORMAT(1X,'            mb          W/m2          W/m2           W/
     &m2          degree/day')
c 9902 FORMAT(1X,I3,3X,F11.6,4X,1P,2(G12.6,2X),G13.6,3X,G16.9,0P) !EWS - not used
 9903 FORMAT(A)
 9910 FORMAT('  Modules and versions used in this calculation:',/,/,5X,
     *        '    rrtm.f: ',6X,A15,10X, 'rrtatm.f: ',6X,A15,/,5X,
     *        '     rtr.f: ',6X,A15,10X, 'rtrcld.f: ',6X,A15,/,5X,
     *        '   rtreg.f: ',6X,A15,8X, 'rtregcld.f: ',6X,A15,/,5X, 
     *        'rtrcldmr.f: ',6x,A15,8x, 'rtregcldmr.f:',5x,A15,/,5x,
     *        ' setcoef.f: ',6X,A15,9X, 'cldprop.f: ',6X,A15,/,5X,
     *        'util_xxx.f: ',6X,A15,10X, 'taumol.f: ',6X,A15,/,5X,
     *        '  k_gB01.f: ',6X,A15,10X, 'k_gB02.f: ',6X,A15,/,5X,
     *        '  k_gB03.f: ',6X,A15,10X, 'k_gB04.f: ',6X,A15,/,5X,
     *        '  k_gB05.f: ',6X,A15,10X, 'k_gB06.f: ',6X,A15,/,5X,
     *        '  k_gB07.f: ',6X,A15,10X, 'k_gB08.f: ',6X,A15,/,5X,
     *        '  k_gB09.f: ',6X,A15,10X, 'k_gB10.f: ',6X,A15,/,5X,
     *        '  k_gB11.f: ',6X,A15,10X, 'k_gB12.f: ',6X,A15,/,5X,
     *        '  k_gB13.f: ',6X,A15,10X, 'k_gB14.f: ',6X,A15,/,5X,
     *        '  k_gB15.f: ',6X,A15,10X, 'k_gB16.f: ',6X,A15,/)

c ojo New SUBROUTINE added, STOP commented and RETURN added
       CALL TranslateSO
c      STOP
       RETURN
       
      END

C************************  SUBROUTINE READPROF  *****************************C

      SUBROUTINE READPROF
                    
C     Read in atmospheric profile.

      IMPLICIT DOUBLE PRECISION (V)
                    
      PARAMETER (MXLAY=203)
      PARAMETER (NBANDS = 16)
      PARAMETER (MAXINPX=35)
      PARAMETER (MAXXSEC=4)
      PARAMETER (MAXPROD = MXLAY*MAXXSEC)

c ojo  ALTZ eliminated. Original statement commented
      DIMENSION IXTRANS(14),SEMIS(NBANDS)
C     DIMENSION ALTZ(0:MXLAY), IXTRANS(14),SEMIS(NBANDS)

      COMMON /CONTROL/  NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SURFACE/  TBOUND,IREFLECT,SEMISS(NBANDS)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBRODL(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /IFIL/     IRD,IPR,IPU,IDUM(15)
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /PATHX/    IXMAX,NXMOL0,IXINDX0(MAXINPX),WX0(MAXINPX,MXLAY)
      COMMON /XRRTATM/  IXSECT

c ojo COMMON BLOCK added      
      COMMON /ALTI/ ALTZ(0:MXLAY)
      
      CHARACTER*80 FORM1(0:1),FORM2(0:1),FORM3(0:1)
      CHARACTER*1 CDOLLAR, CPRCNT ! variables CDUM, CTEST not used

      DATA CDOLLAR /'$'/
      DATA CPRCNT /'%'/
      DATA IXTRANS /0,0,0,1,2,3,0,0,0,0,0,4,0,0/
c      DATA WX /MAXPROD*0.0/

      FORM1(0) = '(3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))'
      FORM2(0) = '(3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))'
      FORM3(0) = '(8E10.3)'
      FORM1(1) = '(G15.7,G10.4,G10.4,A3,I2,1X,2(G7.2,G8.3,G7.2))'
      FORM2(1) = '(G15.7,G10.4,G10.4,A3,I2,23X,(G7.2,G8.3,G7.2))'
      FORM3(1) = '(8G15.7)'

      IXMAX = MAXINPX

c 1000 CONTINUE
 
c ojo ALL READ lines commented
c      READ (IRD,9010,END=8800) CTEST
c     IF (CTEST .EQ. CPRCNT) GO TO 8900 
c     IF (CTEST .NE. CDOLLAR) GO TO 1000

c      READ (IRD,9011) IATM, IXSECT, NUMANGS, IOUT, ICLD
c     If numangs set to -1, reset to default rt code for
c     backwards compatibility with original rrtm
c      IF (NUMANGS .EQ. -1) NUMANGS = 0

C     If clouds are present, read in appropriate input file, IN_CLD_RRTM.
c      IF (ICLD .GE. 1) CALL READCLD

C     Read in surface information.
c      READ (IRD,9012) TBOUND,IEMISS,IREFLECT,(SEMIS(I),I=1,16)

c ojo Definition of INPUT parameters
      NUMANGS = 0
      IEMISS = 0
      IPTHAK = 3
      IXSECT = 0
      do j = 1,NBANDS
        SEMIS(j) = 1.
      end do
      IREFLECT = 0 
c ojo  END definitions

      DO 1500 IBAND = 1, NBANDS
         SEMISS(IBAND) = 1.0
         IF (IEMISS .EQ. 1 .AND. SEMIS(1) .NE. 0.) THEN
            SEMISS(IBAND) = SEMIS(1)
         ELSEIF (IEMISS .EQ. 2) THEN
            IF (SEMIS(IBAND) .NE. 0.) THEN
               SEMISS(IBAND) = SEMIS(IBAND)
            ENDIF
         ENDIF
 1500 CONTINUE

c ojo Lines commented
c      IF (IATM .EQ. 0) THEN
c         READ (IRD,9013) IFORM,NLAYERS,NMOL
c         IF (NMOL.EQ.0) NMOL = 7                                    
c         READ (IRD,FORM1(IFORM)) PAVEL(1),TAVEL(1),SECNTK,CINP,
c     &        IPTHAK,ALTZ(0),PZ(0),TZ(0),ALTZ(1),PZ(1),TZ(1)
c         READ (IRD,FORM3(IFORM)) (WKL(M,1),M=1,7), WBRODL(1)
c         IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,1),M=8,NMOL)
c
c         DO 2000 L = 2, NLAYERS
c            READ (IRD,FORM2(IFORM)) PAVEL(L),TAVEL(L),SECNTK,CINP,
c     &           IPTHRK,ALTZ(L),PZ(L),TZ(L)
c            READ (IRD,FORM3(IFORM)) (WKL(M,L),M=1,7), WBRODL(L)
c           IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,L),M=8,NMOL)
c 2000    CONTINUE       
c           
c         IF (IXSECT .EQ. 1) THEN                                 
c            READ (IRD,9300) NXMOL0
c            NXMOL = NXMOL0
c            CALL XSIDENT(IRD)
c            READ (IRD,9301) IFORMX
C     
c            DO 3000 L = 1, NLAYERS       
c               READ (IRD,9010) CDUM
c               READ (IRD, FORM3(IFORMX)) (WX0(M,L),M=1,7),WBRODX    
c               IF (NXMOL0 .GT. 7) READ (IRD,FORM3(IFORMX)) 
c     &              (WX0(M,L),M=8,NXMOL0)
c 3000       CONTINUE
c         ENDIF
c      ELSE
c         IPU = 7
c         IPR = 66
c         OPEN(UNIT=IPR,FILE='TAPE6',STATUS='UNKNOWN')
c         CALL RRTATM
c         IF (IXSECT .EQ. 1) THEN
c            DO 3300 MX = 1, NXMOL0
c               IXINDX(MX) = IXTRANS(IXINDX0(MX))
c 3300       CONTINUE
c         ENDIF
c      ENDIF
c      IF (TBOUND .LT. 0) TBOUND = TZ(0)

c ojo CALL added

      CALL TRANSLATESI

C     Test for mixing ratio input.
      IMIX = 1
      DO 3500 M = 1, NMOL
         IF (WKL(M,1) .GT. 1.0) THEN
            IMIX = 0
            GO TO 3600
         ENDIF
 3500 CONTINUE
 3600 CONTINUE

      IF (IXSECT .EQ. 1) THEN
         IMIXX = 0
         IF (WX0(1,1) .LE. 1.0) IMIXX = 1
      ENDIF
      DO 5000 L = 1, NLAYERS
         SUMMOL = 0.0
         DO 4100 IMOL = 2, NMOL
            SUMMOL = SUMMOL + WKL(IMOL,L)
 4100    CONTINUE
         IF (IMIX .EQ. 1) THEN
            COLDRY(L) = WBRODL(L) / (1. - SUMMOL)
            DO 4200 IMOL = 1, NMOL
               WKL(IMOL,L) = COLDRY(L) * WKL(IMOL,L)
 4200       CONTINUE
         ELSE
            COLDRY(L) = WBRODL(L) + SUMMOL
         ENDIF

c ojo Next lines commented
c         IF (IXSECT .EQ. 1) THEN
c            DO 4400 IX = 1, NXMOL0
c               IF (IXINDX(IX) .NE. 0) THEN
c                  IF (IMIXX .EQ. 1) THEN
c                     WX(IXINDX(IX),L) = COLDRY(L) * WX0(IX,L) * 1.E-20
c                  ELSE
c                     WX(IXINDX(IX),L) = WX0(IX,L) * 1.E-20
c                  ENDIF
c               ENDIF
               
c 4400       CONTINUE
c         ENDIF
 5000 CONTINUE

      GO TO 9000

c 8800 CONTINUE
c 8900 IF (CTEST.EQ.'%') STOP 'END OF INPUT FILE'
 9000 CONTINUE

c 9010 FORMAT (A1)
c 9011 FORMAT (49X,I1,19X,I1,13X,I2,2X,I3,4X,I1)
c 9012 FORMAT (E10.3,1X,I1,2X,I1,16E5.3)
c 9013 FORMAT (1X,I1,I3,I5)                                     
c 9300 FORMAT (I5)
c 9301 FORMAT (1X,I1)

      RETURN
      END 

C************************  SUBROUTINE READCLD  *****************************C

      SUBROUTINE READCLD

C     Purpose:  To read in IN_CLD_RRTM, the file that contains input 
C               cloud properties.

      PARAMETER (MXLAY=203)
      PARAMETER (NBANDS = 16)

      COMMON /PROFI/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)

      CHARACTER*1 CTEST, CPERCENT

      DATA CPERCENT /'%'/
      IRDCLD = 11

      OPEN(IRDCLD,FILE='IN_CLD_RRTM',FORM='FORMATTED')

C     Read in cloud input option.  
      READ(IRDCLD,9050) INFLAG, ICEFLAG, LIQFLAG
      DO 500 LAY = 1, NLAYERS
         CLDFRAC(LAY) = 0.
 500  CONTINUE

 1000 CONTINUE
C     For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
C     DAT1 are pertinent.  If CTEST = '%', then there are no more 
C     cloudy layers to process.
      READ (IRDCLD,9100,END=9000) CTEST,LAY,FRAC,DAT1,DAT2,DAT3,DAT4
      IF (CTEST .EQ. CPERCENT) GO TO 9000
      CLDFRAC(LAY) = FRAC
      CLDDAT1(LAY) = DAT1
      CLDDAT2(LAY) = DAT2
      CLDDAT3(LAY) = DAT3
      CLDDAT4(LAY) = DAT4
      GO TO 1000

 9000 CONTINUE
      CLOSE(IRDCLD)

 9050 FORMAT (3X,I2,4X,I1,4X,I1)
 9100 FORMAT (A1,1X,I3,5E10.5)

      RETURN
      END

C************************  SUBROUTINE XSIDENT  *****************************C

      SUBROUTINE XSIDENT(IRD)
C                    
C     This subroutine identifies which cross-sections are to be used.

      PARAMETER (MAXINPX=35)
      PARAMETER (MAXXSEC=4)

      IMPLICIT DOUBLE PRECISION (V)  ! 
C                    
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
C                    
C     NXMOL     - number of cross-sections input by user
C     IXINDX(I) - index of cross-section molecule corresponding to Ith
C                 cross-section specified by user
C                 = 0 -- not allowed in RRTM
C                 = 1 -- CCL4
C                 = 2 -- CFC11
C                 = 3 -- CFC12
C                 = 4 -- CFC22
C                    
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          
C                    
      CHARACTER*10 XSNAME(MAXINPX),ALIAS(MAXXSEC,4),BLANK               
C                    
      DATA (ALIAS(1,I),I=1,4)/                                           
     *    'CCL4      ', 'CCL3F     ', 'CCL2F2    ', 'CHCLF2    '/
      DATA (ALIAS(2,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFCL3     ', 'CF2CL2    ', 'CHF2CL    '/ 
      DATA (ALIAS(3,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFC11     ', 'CFC12     ', 'CFC22     '/
      DATA (ALIAS(4,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'F11       ', 'F12       ', 'F22       '/

      DATA BLANK / '          '/
C                    
      DO 10 I = 1, NXMOL
         XSNAME(I) = BLANK
   10 CONTINUE       
C                    
C     READ IN THE NAMES OF THE MOLECULES                                  
C                    
      IF (NXMOL.GT.7) THEN
         READ (IRD,'(7A10)') (XSNAME(I),I=1,7)
         READ (IRD,'(8A10)') (XSNAME(I),I=8,NXMOL)
      ELSE           
         READ (IRD,'(7A10)') (XSNAME(I),I=1,NXMOL)
      ENDIF          
C                    
C     MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS           
C     AND DETERMINE THE INDEX VALUE.  
      IXMAX = 4     
      DO 40 I = 1, NXMOL
C        Left-justify all inputed names.                                      
         CALL CLJUST (XSNAME(I),10)
         IXINDX(I) = 0
         DO 20 J = 1, IXMAX
            IF ((XSNAME(I).EQ.ALIAS(1,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(2,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(3,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(4,J))) THEN
               IXINDX(I) = J
            ENDIF    
   20    CONTINUE
   40 CONTINUE       

      RETURN
      END

      BLOCK DATA

      PARAMETER (MG = 16)
      PARAMETER (NBANDS = 16)
      PARAMETER (MXLAY=203)
      PARAMETER (MAXXSEC=4)
      PARAMETER (MAXPROD = MXLAY*MAXXSEC)

      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(NBANDS),NSPA(MG),NSPB(MG)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)


      COMMON /CVRRTM/    HVRRTM
      COMMON /CVRREG/    HVRREG
      COMMON /CVRRTR/    HVRRTR
      COMMON /CVRATM/    HVRATM
      COMMON /CVRSET/    HVRSET
      COMMON /CVRTAU/    HVRTAU
      COMMON /CVRRGC/    HVRRGC
      COMMON /CVRRTC/    HVRRTC
      COMMON /CVRCLD/    HVRCLD
      COMMON /CVRUTL/    HVRUTL
      COMMON /CVREXT/    HVREXT
      COMMON /CVRRTX/    HVRRTX
      COMMON /CVRRGX/    HVRRGX

      COMMON /HVRSN1/    HVRKG1
      COMMON /HVRSN2/    HVRKG2
      COMMON /HVRSN3/    HVRKG3
      COMMON /HVRSN4/    HVRKG4
      COMMON /HVRSN5/    HVRKG5
      COMMON /HVRSN6/    HVRKG6
      COMMON /HVRSN7/    HVRKG7
      COMMON /HVRSN8/    HVRKG8
      COMMON /HVRSN9/    HVRKG9
      COMMON /HVRSN10/   HVRKG10
      COMMON /HVRSN11/   HVRKG11
      COMMON /HVRSN12/   HVRKG12
      COMMON /HVRSN13/   HVRKG13
      COMMON /HVRSN14/   HVRKG14
      COMMON /HVRSN15/   HVRKG15
      COMMON /HVRSN16/   HVRKG16

      CHARACTER*15 HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVRRGC,HVRRTC,HVRCLD,HVRUTL,HVREXT,
     *            HVRRTX,HVRRGX

      CHARACTER*15 HVRKG1,HVRKG2,HVRKG3,HVRKG4,HVRKG5
      CHARACTER*15 HVRKG6,HVRKG7,HVRKG8,HVRKG9,HVRKG10,HVRKG11
      CHARACTER*15 HVRKG12,HVRKG13,HVRKG14,HVRKG15,HVRKG16

      DATA WAVENUM1(1) /10./, WAVENUM2(1) /350./, DELWAVE(1) /340./
      DATA WAVENUM1(2) /350./, WAVENUM2(2) /500./, DELWAVE(2) /150./
      DATA WAVENUM1(3) /500./, WAVENUM2(3) /630./, DELWAVE(3) /130./
      DATA WAVENUM1(4) /630./, WAVENUM2(4) /700./, DELWAVE(4) /70./
      DATA WAVENUM1(5) /700./, WAVENUM2(5) /820./, DELWAVE(5) /120./
      DATA WAVENUM1(6) /820./, WAVENUM2(6) /980./, DELWAVE(6) /160./
      DATA WAVENUM1(7) /980./, WAVENUM2(7) /1080./, DELWAVE(7) /100./
      DATA WAVENUM1(8) /1080./, WAVENUM2(8) /1180./, DELWAVE(8) /100./
      DATA WAVENUM1(9) /1180./, WAVENUM2(9) /1390./, DELWAVE(9) /210./
      DATA WAVENUM1(10) /1390./,WAVENUM2(10) /1480./,DELWAVE(10) /90./
      DATA WAVENUM1(11) /1480./,WAVENUM2(11) /1800./,DELWAVE(11) /320./
      DATA WAVENUM1(12) /1800./,WAVENUM2(12) /2080./,DELWAVE(12) /280./
      DATA WAVENUM1(13) /2080./,WAVENUM2(13) /2250./,DELWAVE(13) /170./
      DATA WAVENUM1(14) /2250./,WAVENUM2(14) /2380./,DELWAVE(14) /130./
      DATA WAVENUM1(15) /2380./,WAVENUM2(15) /2600./,DELWAVE(15) /220./
      DATA WAVENUM1(16) /2600./,WAVENUM2(16) /3250./,DELWAVE(16) /650./

      DATA NG /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      DATA NSPA /1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/
      DATA NSPB /1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
      DATA HEATFAC /8.4391/

      DATA WX /MAXPROD*0.0/

      DATA HVRRTM / 'NOT USED' /,   HVRREG / 'NOT USED' /,
     *     HVRRTR / 'NOT USED' /,   HVRATM / 'NOT USED' /,
     *     HVRSET / 'NOT USED' /,   HVRTAU / 'NOT USED' /,
     *     HVRUTL / 'NOT USED' /, HVREXT / 'NOT USED' /,
     *     HVRRTC /'NOT USED'/, HVRRTX /'NOT USED'/,
     *     HVRRGC /'NOT USED'/, HVRRGX /'NOT USED'/,
     *     HVRCLD /'NOT USED'/

c     DATA HVRSN1 / 'NOT USED' /, HVRSN2 / 'NOT USED' /,
c          HVRSN3 / 'NOT USED' /, HVRSN4 / 'NOT USED' /,
c          HVRSN5 / 'NOT USED' /, HVRSN6 / 'NOT USED' /,
c          HVRSN7 / 'NOT USED' /, HVRSN8 / 'NOT USED' /,
c          HVRSN9 / 'NOT USED' /, HVRSN10/ 'NOT USED' /,
c          HVRSN11/ 'NOT USED' /, HVRSN12/ 'NOT USED' /,
c          HVRSN13/ 'NOT USED' /, HVRSN14/ 'NOT USED' /,
c          HVRSN15/ 'NOT USED' /, HVRSN16/ 'NOT USED' /


      END
c**********************************************************************
      Block Data phys_consts
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
c
      DATA PI / 3.1415927410125732 /
c
c    Constants from NIST 01/11/2002
c
      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     *     CLIGHT / 2.99792458E+10 /, 
     *     AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     *     GASCON / 8.314472  E+07 /
     *     RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c
c     Pi was obtained from   PI = 2.*ASIN(1.)                             A03980
c
c     units are genrally cgs
c
c     The first and second radiation constants are taken from NIST.
c     They were previously obtained from the relations:
c                            RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      A03990
c                            RADCN2 = PLANCK*CLIGHT/BOLTZ                 A04000
      end
c



c ojo SUBROUTINES TRANSLATESI and TRANSLATESO added
c (made by Kara Krelove)
C-------------------------------------------------------
      SUBROUTINE TRANSLATESI
C  This subroutine is designed to set some basic parameters for 
C   Mlawer's RRTM code, in lieu of actually reading them in in the
C   original subroutine, since that part of the code is out. 

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NL=ND-1, MXLAY=203, NBANDS=16)
c     PARAMETER (ND = 52, NL=51, MXLAY=203, NBANDS=16)

      COMMON/ MLAWERI/  layers, numspec, newalt(ND), TempT(0:NL), 
     &                         Pres(0:NL), gasses(7, 0:NL), COLDEP(ND)
      COMMON /PROFI/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY), TZ(0:MXLAY)
      COMMON /SURFACE/   TBOUND,IREFLECT,SEMISS(NBANDS)
      COMMON /SPECIES/   COLDRY(MXLAY),WKL(35,MXLAY),WBRODL(MXLAY),
     &                   COLMOL(MXLAY),NMOL
      COMMON /ALTI/  altz(0:MXLAY)


        integer x, z
        real newalt

        do x = 0, layers
          TZ(x) = TempT(x)
          PZ(x) = Pres(x)
          WBRODL(x) = COLDEP(x)
          altz(x) = newalt(x)
          do z = 1, 7
             wkl(z, x+1) = gasses(z,x)
          end do
        end do

        altz(0) = 0.0

        NMOL = 7
        NLAYERS = layers

C  calculate average Ts and Ps

        do x = 1, NLAYERS
          TAVEL(x) = (TZ(x-1) + TZ(x))/2.
          PAVEL(x) = (PZ(x-1) + PZ(x))/2.
        end do
        
        TBOUND = TAVEL(1)

      return

      END
C--------------------------------------------------------
      SUBROUTINE TRANSLATESO
C  This subroutine is designed to translate Mlawer's results
C   back into a format SurfTem can use, including flipping the 
C   ordering of the layers. 

      INCLUDE 'CLIMA/INCLUDE/header.inc'
c     PARAMETER (ND=52, MXLAY=203)
      PARAMETER (MXLAY=203)

      COMMON /IRBLK/     FUPIR(ND),FDNIR(ND)
      COMMON /PROFI/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY), TZ(0:MXLAY)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY)


C  align Temp-grid points with SurfTem's flux grid points by taking average;
C  puts them halfway between T-grid points. 

C  This is temporarily commented out due to the problem of not having
C  an extra layer with which to properly calculate a midpoint for the 
C  topmost layer. 

C      do n = 1, NLAYERS
C        tempUflux(n) = (TOTUFLUX(n-1) + TOTUFLUX(n))/2.
C        tempDflux(n) = (TOTDFLUX(n-1) + TOTDFLUX(n))/2.
C      end do

C        tempUflux(0)=TOTUFLUX(0)
C        tempDflux(0)=TOTDFLUX(0)

C   Previous line of code not quite correct. TOT and F*IR have same # of
C   points, and therefor I have to substitute at the bottom. 
 
C  Now reverse the order of the layers.

      do n = 0, NLAYERS
        FUPIR(n+1) = TOTUFLUX(NLAYERS-n)*1000
        FDNIR(n+1) = TOTDFLUX(NLAYERS-n)*1000
      end do

      return

      END
C--------------------------------------------------------
