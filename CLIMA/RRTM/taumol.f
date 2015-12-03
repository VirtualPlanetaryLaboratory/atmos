C     path:      $Source: /storm/rc1/cvsroot/rc/rrtm_lw/src/taumol.f,v $
C     author:    $Author: jdelamer $
C     revision:  $Revision: 3.3 $
C     created:   $Date: 2002/08/15 18:33:27 $
*******************************************************************************
*                                                                             *
*                  Optical depths developed for the                           *
*                                                                             *
*                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
*                                                                             *
*                                                                             *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
*                        131 HARTWELL AVENUE                                  *
*                        LEXINGTON, MA 02421                                  *
*                                                                             *
*                                                                             *
*                           ELI J. MLAWER                                     * 
*                         JENNIFER DELAMERE                                   * 
*                         STEVEN J. TAUBMAN                                   *
*                         SHEPARD A. CLOUGH                                   *
*                                                                             *
*                                                                             *
*                                                                             *
*                                                                             *
*                       email:  mlawer@aer.com                                *
*                       email:  jdelamer@aer.com                              *
*                                                                             *
*        The authors wish to acknowledge the contributions of the             *
*        following people:  Karen Cady-Pereira, Patrick D. Brown,             *  
*        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
*                                                                             *
*******************************************************************************
*     TAUMOL                                                                  *
*                                                                             *
*     This file contains the subroutines TAUGBn (where n goes from            *
*     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
*     per g-value and layer for band n.                                       *
*                                                                             *
*  Output:  optical depths (unitless)                                         *
*           fractions needed to compute Planck functions at every layer       *
*               and g-value                                                   *
*                                                                             *
*     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
*     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
*                                                                             *
*  Input                                                                      *
*                                                                             *
*     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
*     COMMON /PRECISE/  ONEMINUS                                              *
*     COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
*     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
*     COMMON /PROFDATA/ LAYTROP,                                              *
*    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
*    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
*    &                  COLO2(MXLAY)
*     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
*    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
*     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
*     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
*                                                                             *
*     Description:                                                            *
*     NG(IBAND) - number of g-values in band IBAND                            *
*     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
*                   atmospheres that are stored for band IBAND per            *
*                   pressure level and temperature.  Each of these            *
*                   atmospheres has different relative amounts of the         *
*                   key species for the band (i.e. different binary           *
*                   species parameters).                                      *
*     NSPB(IBAND) - same for upper atmosphere                                 *
*     ONEMINUS - since problems are caused in some cases by interpolation     *
*                parameters equal to or greater than 1, for these cases       *
*                these parameters are set to this value, slightly < 1.        *
*     PAVEL - layer pressures (mb)                                            *
*     TAVEL - layer temperatures (degrees K)                                  *
*     PZ - level pressures (mb)                                               *
*     TZ - level temperatures (degrees K)                                     *
*     LAYTROP - layer at which switch is made from one combination of         *
*               key species to another                                        *
*     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
*               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
*               respectively (molecules/cm**2)                                *

*     FACij(LAY) - for layer LAY, these are factors that are needed to        *
*                  compute the interpolation factors that multiply the        *
*                  appropriate reference k-values.  A value of 0 (1) for      *
*                  i,j indicates that the corresponding factor multiplies     *
*                  reference k-value for the lower (higher) of the two        *
*                  appropriate temperatures, and altitudes, respectively.     *
*     JP - the index of the lower (in altitude) of the two appropriate        *
*          reference pressure levels needed for interpolation                 *
*     JT, JT1 - the indices of the lower of the two appropriate reference     *
*               temperatures needed for interpolation (for pressure           *
*               levels JP and JP+1, respectively)                             *
*     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
*               (water vapor density)/(atmospheric density at 296K and        *
*               1013 mb)                                                      *
*     SELFFRAC - factor needed for temperature interpolation of reference     *
*                water vapor self-continuum data                              *
*     INDSELF - index of the lower of the two appropriate reference           *
*               temperatures needed for the self-continuum interpolation      *
*     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
*     FORFRAC - factor needed for temperature interpolation of reference      *
*                water vapor foreign-continuum data                           *
*     INDFOR  - index of the lower of the two appropriate reference           *
*               temperatures needed for the foreign-continuum interpolation   *
*                                                                             *
*  Data input                                                                 *
*     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
*                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
*        (note:  n is the band number,'MGAS' is the species name of the minor *
*         gas)                                                                *
*                                                                             *
*     Description:                                                            *
*     KA - k-values for low reference atmospheres (key-species only)          *
*          (units: cm**2/molecule)                                            *
*     KB - k-values for high reference atmospheres (key-species only)         *
*          (units: cm**2/molecule)                                            *
*     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
*          (units: cm**2/molecule)                                            *
*     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
*          (units: cm**2/molecule)                                            *
*     SELFREF - k-values for water vapor self-continuum for reference         *
*               atmospheres (used below LAYTROP)                              *
*               (units: cm**2/molecule)                                       *
*     FORREF  - k-values for water vapor foreign-continuum for reference      *
*               atmospheres (used below/above LAYTROP)                        *
*               (units: cm**2/molecule)                                       *
*                                                                             *
*     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
*     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
*                                                                             *
*******************************************************************************


      SUBROUTINE TAUGB1

C     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

C     BAND 1:  10-350 cm-1 (low key - H2O; low minor - N2)
C                          (high key - H2O; high minor - N2)

C     NOTE: Previous versions of RRTM BAND 1: 
C           10-250 cm-1 (low - H2O; high - H2O)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K1/       KA(5,13,MG), KB(5,13:59,MG),FORREF(4,MG), 
     &                  SELFREF(10,MG),KA_MN2(19,MG),KB_MN2(19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)


C Planck fraction mapping level: P = 212.7250 mbar, T = 223.06 K
      DATA FRACREFA /
     &2.1227E-01,1.8897E-01,1.3934E-01,1.1557E-01,9.5282E-02,8.3359E-02,
     &6.5333E-02,5.2016E-02,3.4272E-02,4.0257E-03,3.1857E-03,2.6014E-03,
     &1.9141E-03,1.2612E-03,5.3169E-04,7.6476E-05/

C Planck fraction mapping level: P = 212.7250 mbar, T = 223.06 K
C These planck refractions were calculated using lower 
C atmosphere parameters.
      DATA FRACREFB /
     &2.1227E-01,1.8897E-01,1.3934E-01,1.1557E-01,9.5282E-02,8.3359E-02,
     &6.5333E-02,5.2016E-02,3.4272E-02,4.0257E-03,3.1857E-03,2.6014E-03,
     &1.9141E-03,1.2612E-03,5.3169E-04,7.6476E-05/

C Minor gas mapping levels:
C     LOWER - N2, P = 142.5490 mbar, T = 215.70 K
C     UPPER - N2, P = 142.5490 mbar, T = 215.70 K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA, KB, KA_MN2, KB_MN2, MINORFRAC

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.

      HVRTAU = '$Revision: 3.3 $'
      DO 2500 LAY = 1, LAYTROP

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(1) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(1) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         PP = PAVEL(LAY)
         CORRADJ =  1.
         IF (PP .LT. 250.) THEN
            CORRADJ = 1. - 0.15 * (250.-PP) / 154.4
         ENDIF

         SCALEN2 = COLBRD(LAY) * SCALEMINORN2(LAY)
         DO 2000 IG = 1, NG(1)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUN2 = SCALEN2*(KA_MN2(INDM,IG) +  
     &           MINORFRAC(LAY) *
     &           (KA_MN2(INDM+1,IG) - KA_MN2(INDM,IG)))
            TAUG(LAY,IG) = CORRADJ * (COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG))  
     &           + TAUSELF + TAUFOR
     &           + TAUN2)
             FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(1) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(1) + 1
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         PP = PAVEL(LAY)
         CORRADJ =  1. - 0.15 * (PP / 95.6)

         SCALEN2 = COLBRD(LAY) * SCALEMINORN2(LAY)
         DO 3000 IG = 1, NG(1)
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUN2 = SCALEN2*(KB_MN2(INDM,IG) +  
     &           MINORFRAC(LAY) *
     &           (KB_MN2(INDM+1,IG) - KB_MN2(INDM,IG)))
            TAUG(LAY,IG) = CORRADJ * (COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))   
     &           + TAUFOR
     &           + TAUN2)
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------


C----------------------------------------------------------------------------

      SUBROUTINE TAUGB2

C     BAND 2:  350-500 cm-1 (low key - H2O; high key - H2O)

C     NOTE: Previous version of RRTM BAND 2: 
C           250 - 500 cm-1 (low - H2O; high - H2O)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K2/       KA(5,13,MG), KB(5,13:59,MG) , FORREF(4,MG), 
     &                  SELFREF(10,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

C Planck fraction mapping level: P = 1053.630 mbar, T = 294.2 K

      DATA FRACREFA /
     & 1.6388E-01, 1.5241E-01, 1.4290E-01, 1.2864E-01, 
     & 1.1615E-01, 1.0047E-01, 8.0013E-02, 6.0445E-02, 
     & 4.0530E-02, 4.3879E-03, 3.5726E-03, 2.7669E-03,
     & 2.0078E-03, 1.2864E-03, 4.7630E-04, 6.9109E-05/

C Planck fraction mapping level: P = 3.206e-2 mb, T = 197.92 K
      DATA FRACREFB /
     & 1.4697E-01, 1.4826E-01, 1.4278E-01, 1.3320E-01, 
     & 1.1965E-01, 1.0297E-01, 8.4170E-02, 6.3282E-02, 
     & 4.2868E-02, 4.6644E-03, 3.8619E-03, 3.0533E-03,
     & 2.2359E-03, 1.4226E-03, 5.3642E-04, 7.6316E-05/


      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.
      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(2) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(2) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         PP = PAVEL(LAY)
         CORRADJ = 1. - .05 * (PP - 100.) / 900.
         DO 2000 IG = 1, NG(2)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUG(LAY,IG) = CORRADJ * (COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR)
            FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(2) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(2) + 1
         INDF = INDFOR(LAY)
         DO 3000 IG = 1, NG(2)
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
     &           + TAUFOR
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB3

C     BAND 3:  500-630 cm-1 (low key - H2O,CO2; low minor - n2o)
C                           (high key - H2O,CO2; high minor - n2o)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K3/       KA(9,5,13,MG), KB(5,5,13:59,MG), FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MN2O(9,19,MG), 
     &                  KB_MN2O(5,19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      REAL KA,KB
      REAL KA_MN2O, KB_MN2O, MINORFRAC
      REAL N2OM1,N2OM2
      DIMENSION ABSA(585,MG),ABSB(1175,MG)
      DIMENSION FRACREFA(MG,9), FRACREFB(MG,5)

C Planck fraction mapping level: P=212.7250 mbar, T = 223.06 K
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.6251E-01,1.5572E-01,1.4557E-01,1.3208E-01,1.1582E-01,9.6895E-02,
     &7.8720E-02,5.8462E-02,3.9631E-02,4.3001E-03,3.5555E-03,2.8101E-03,
     &2.0547E-03,1.3109E-03,4.9403E-04,6.9515E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.6006E-01,1.5576E-01,1.4609E-01,1.3276E-01,1.1594E-01,9.7336E-02,
     &7.9035E-02,5.8696E-02,3.9723E-02,4.3001E-03,3.5555E-03,2.8101E-03,
     &2.0547E-03,1.3109E-03,4.9403E-04,6.9515E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.5952E-01,1.5566E-01,1.4590E-01,1.3294E-01,1.1599E-01,9.7511E-02,
     &7.9127E-02,5.8888E-02,3.9874E-02,4.3001E-03,3.5555E-03,2.8102E-03,
     &2.0547E-03,1.3109E-03,4.9403E-04,6.9515E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.5907E-01,1.5541E-01,1.4585E-01,1.3316E-01,1.1596E-01,9.7647E-02,
     &7.9243E-02,5.9024E-02,4.0028E-02,4.3112E-03,3.5555E-03,2.8102E-03,
     &2.0547E-03,1.3109E-03,4.9403E-04,6.9515E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.5862E-01,1.5517E-01,1.4588E-01,1.3328E-01,1.1585E-01,9.7840E-02,
     &7.9364E-02,5.9174E-02,4.0160E-02,4.3403E-03,3.5900E-03,2.8102E-03,
     &2.0547E-03,1.3109E-03,4.9403E-04,6.9515E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.5830E-01,1.5490E-01,1.4582E-01,1.3331E-01,1.1567E-01,9.8079E-02,
     &7.9510E-02,5.9369E-02,4.0326E-02,4.3343E-03,3.5908E-03,2.8527E-03,
     &2.0655E-03,1.3109E-03,4.9403E-04,6.9515E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.5789E-01,1.5435E-01,1.4595E-01,1.3304E-01,1.1566E-01,9.8426E-02,
     &7.9704E-02,5.9618E-02,4.0520E-02,4.3812E-03,3.6147E-03,2.8395E-03,
     &2.1301E-03,1.3145E-03,4.9403E-04,6.9515E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.5704E-01,1.5398E-01,1.4564E-01,1.3222E-01,1.1586E-01,9.9230E-02,
     &8.0011E-02,6.0149E-02,4.0790E-02,4.4253E-03,3.6534E-03,2.9191E-03,
     &2.1373E-03,1.3558E-03,5.1631E-04,7.8794E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.5270E-01,1.5126E-01,1.4264E-01,1.3106E-01,1.1740E-01,1.0137E-01,
     &8.3057E-02,6.2282E-02,4.2301E-02,4.6486E-03,3.8159E-03,3.0472E-03,
     &2.2870E-03,1.4818E-03,5.6773E-04,7.8794E-05/
C Planck fraction mapping level: P = 95.8 mbar, T = 215.7 K
      DATA (FRACREFB(IG, 1),IG=1,16) /
     &1.6413E-01,1.5665E-01,1.4606E-01,1.3184E-01,1.1517E-01,9.6243E-02,
     &7.7982E-02,5.8165E-02,3.9311E-02,4.2586E-03,3.5189E-03,2.7793E-03,
     &2.0376E-03,1.2938E-03,4.8853E-04,6.8745E-05/
      DATA (FRACREFB(IG, 2),IG=1,16) /
     &1.6254E-01,1.5674E-01,1.4652E-01,1.3221E-01,1.1535E-01,9.6439E-02,
     &7.8155E-02,5.8254E-02,3.9343E-02,4.2586E-03,3.5189E-03,2.7793E-03,
     &2.0376E-03,1.2938E-03,4.8853E-04,6.8745E-05/
      DATA (FRACREFB(IG, 3),IG=1,16) /
     &1.6177E-01,1.5664E-01,1.4669E-01,1.3242E-01,1.1541E-01,9.6536E-02,
     &7.8257E-02,5.8387E-02,3.9431E-02,4.2587E-03,3.5189E-03,2.7793E-03,
     &2.0376E-03,1.2938E-03,4.8853E-04,6.8745E-05/
      DATA (FRACREFB(IG, 4),IG=1,16) /
     &1.6077E-01,1.5679E-01,1.4648E-01,1.3273E-01,1.1546E-01,9.6779E-02,
     &7.8371E-02,5.8546E-02,3.9611E-02,4.2772E-03,3.5190E-03,2.7793E-03,
     &2.0376E-03,1.2938E-03,4.8853E-04,6.8745E-05/
      DATA (FRACREFB(IG, 5),IG=1,16) /
     &1.6067E-01,1.5608E-01,1.4247E-01,1.2881E-01,1.1449E-01,9.8802E-02,
     &8.0828E-02,6.0977E-02,4.1494E-02,4.5116E-03,3.7290E-03,2.9460E-03,
     &2.1948E-03,1.3778E-03,5.4552E-04,7.9969E-05/

C Minor gas mapping levels:
C     LOWER - N2O, P = 706.272 mbar, T = 278.94 K
C     UPPER - N2O, P = 95.58 mbar, T = 215.7 K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

C     P = 212.725 mb
      REFRAT_PLANCK_A = CHI_MLS(1,9)/CHI_MLS(2,9)

C     P = 95.58 mb
      REFRAT_PLANCK_B = CHI_MLS(1,13)/CHI_MLS(2,13)

C     P = 706.270mb
      REFRAT_M_A = CHI_MLS(1,3)/CHI_MLS(2,3)

C     P = 95.58 mb 
      REFRAT_M_B = CHI_MLS(1,13)/CHI_MLS(2,13)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water vapor 
C     self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)        

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)
         FMN2OMF = MINORFRAC(LAY)*FMN2O
c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/COLDRY(LAY)
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(3) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(3) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(3)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG) - 
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG) - 
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY) *
     &              (N2OM2 - N2OM1)
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000       CONTINUE
         ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(3)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG) - 
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG) - 
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY) *
     &              (N2OM2 - N2OM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O 
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010          CONTINUE
         ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(3)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG) - 
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG) - 
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY) *
     &              (N2OM2 - N2OM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O   
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
        ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 4.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS1) * FAC01(LAY)
         FAC011 = (1. - FS1) * FAC11(LAY)
         FAC101 = FS1 * FAC01(LAY)
         FAC111 = FS1 * FAC11(LAY)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_B*COLCO2(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 4.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)
         FMN2OMF = MINORFRAC(LAY)*FMN2O
c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/COLDRY(LAY)
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_B*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 4.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(3) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(3) + JS1
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG(3)
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            N2OM1 = KB_MN2O(JMN2O,INDM,IG) + FMN2O*
     &           (KB_MN2O(JMN2O+1,INDM,IG)-KB_MN2O(JMN2O,INDM,IG))
            N2OM2 = KB_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &           (KB_MN2O(JMN2O+1,INDM+1,IG)-KB_MN2O(JMN2O,INDM+1,IG))
            ABSN2O = N2OM1 + MINORFRAC(LAY) * (N2OM2 - N2OM1)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &          FAC100 * ABSB(IND0+1,IG) +
     &          FAC010 * ABSB(IND0+5,IG) +
     &          FAC110 * ABSB(IND0+6,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSB(IND1,IG) + 
     &          FAC101 * ABSB(IND1+1,IG) +
     &          FAC011 * ABSB(IND1+5,IG) +
     &          FAC111 * ABSB(IND1+6,IG)) 
     &          + TAUFOR
     &          + ADJCOLN2O*ABSN2O            
            FRACS(LAY,IG) = FRACREFB(IG,JPL) + FPL *
     &          (FRACREFB(IG,JPL+1)-FRACREFB(IG,JPL))
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB4

C     BAND 4:  630-700 cm-1 (low key - H2O,CO2; high key - O3,CO2)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K4/       KA(9,5,13,MG), KB(5,5,13:59,MG) , FORREF(4,MG), 
     &                  SELFREF(10,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(585,MG),ABSB(1175,MG)
      DIMENSION FRACREFA(MG,9), FRACREFB(MG,6)

C Planck fraction mapping level : P = 142.5940 mbar, T = 215.70 K
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.5572E-01,1.4925E-01,1.4107E-01,1.3126E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2393E-02,4.2146E-02,4.5907E-03,3.7965E-03,2.9744E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.5572E-01,1.4925E-01,1.4107E-01,1.3126E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2392E-02,4.2146E-02,4.5906E-03,3.7965E-03,2.9745E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.5572E-01,1.4925E-01,1.4107E-01,1.3126E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2393E-02,4.2146E-02,4.5907E-03,3.7965E-03,2.9745E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.5572E-01,1.4925E-01,1.4107E-01,1.3126E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2393E-02,4.2146E-02,4.5907E-03,3.7964E-03,2.9744E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.5572E-01,1.4925E-01,1.4107E-01,1.3126E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2393E-02,4.2146E-02,4.5907E-03,3.7965E-03,2.9744E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.5572E-01,1.4925E-01,1.4107E-01,1.3126E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2393E-02,4.2146E-02,4.5907E-03,3.7965E-03,2.9744E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.5572E-01,1.4926E-01,1.4107E-01,1.3126E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2393E-02,4.2146E-02,4.5908E-03,3.7964E-03,2.9745E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.5571E-01,1.4926E-01,1.4107E-01,1.3125E-01,1.1791E-01,1.0173E-01,
     &8.2949E-02,6.2393E-02,4.2146E-02,4.5907E-03,3.7964E-03,2.9744E-03,
     &2.2074E-03,1.4063E-03,5.3012E-04,7.4595E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.5952E-01,1.5155E-01,1.4217E-01,1.3077E-01,1.1667E-01,1.0048E-01,
     &8.1511E-02,6.1076E-02,4.1111E-02,4.4432E-03,3.6910E-03,2.9076E-03,
     &2.1329E-03,1.3566E-03,5.2235E-04,7.9935E-05/

C Planck fraction mapping level : P = 95.58350 mb, T = 215.70 K

      DATA (FRACREFB(IG, 1),IG=1,16) /
     &1.5558E-01,1.4931E-01,1.4104E-01,1.3124E-01,1.1793E-01,1.0160E-01,
     &8.3142E-02,6.2403E-02,4.2170E-02,4.5935E-03,3.7976E-03,2.9986E-03,
     &2.1890E-03,1.4061E-03,5.3005E-04,7.4587E-05/
      DATA (FRACREFB(IG, 2),IG=1,16) /
     &1.5558E-01,1.4932E-01,1.4104E-01,1.3124E-01,1.1792E-01,1.0159E-01,
     &8.3142E-02,6.2403E-02,4.2170E-02,4.5935E-03,3.7976E-03,2.9986E-03,
     &2.1890E-03,1.4061E-03,5.3005E-04,7.4587E-05/
      DATA (FRACREFB(IG, 3),IG=1,16) /
     &1.5558E-01,1.4933E-01,1.4103E-01,1.3124E-01,1.1792E-01,1.0159E-01,
     &8.3142E-02,6.2403E-02,4.2170E-02,4.5935E-03,3.7976E-03,2.9986E-03,
     &2.1890E-03,1.4061E-03,5.3005E-04,7.4587E-05/
      DATA (FRACREFB(IG, 4),IG=1,16) /
     &1.5569E-01,1.4926E-01,1.4102E-01,1.3122E-01,1.1791E-01,1.0159E-01,
     &8.3141E-02,6.2403E-02,4.2170E-02,4.5935E-03,3.7976E-03,2.9986E-03,
     &2.1890E-03,1.4061E-03,5.3005E-04,7.4587E-05/
      DATA (FRACREFB(IG, 5),IG=1,16) /
     &1.5947E-01,1.5132E-01,1.4195E-01,1.3061E-01,1.1680E-01,1.0054E-01,
     &8.1785E-02,6.1212E-02,4.1276E-02,4.4424E-03,3.6628E-03,2.8943E-03,
     &2.1134E-03,1.3457E-03,5.1024E-04,7.3998E-05/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     P =   142.5940 mb
      REFRAT_PLANCK_A = CHI_MLS(1,11)/CHI_MLS(2,11)

C     P = 95.58350 mb
      REFRAT_PLANCK_B = CHI_MLS(3,13)/CHI_MLS(2,13)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water 
C     vapor self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(4) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(4) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(4)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
         ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG(4)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
         ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(4)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
        ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         SPECCOMB = COLO3(LAY) + RAT_O3CO2(LAY)*COLCO2(LAY)
         SPECPARM = COLO3(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLO3(LAY) + RAT_O3CO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLO3(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 4.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS1) * FAC01(LAY)
         FAC011 = (1. - FS1) * FAC11(LAY)
         FAC101 = FS1 * FAC01(LAY)
         FAC111 = FS1 * FAC11(LAY)

         SPECCOMB_PLANCK = COLO3(LAY)+REFRAT_PLANCK_B*COLCO2(LAY)
         SPECPARM_PLANCK = COLO3(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 4.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(4) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(4) + JS1

         DO 3000 IG = 1, NG(4)
            TAUG(LAY,IG) =  SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &          FAC100 * ABSB(IND0+1,IG) +
     &          FAC010 * ABSB(IND0+5,IG) +
     &          FAC110 * ABSB(IND0+6,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSB(IND1,IG) + 
     &          FAC101 * ABSB(IND1+1,IG) +
     &          FAC011 * ABSB(IND1+5,IG) +
     &          FAC111 * ABSB(IND1+6,IG)) 
            FRACS(LAY,IG) = FRACREFB(IG,JPL) + FPL *
     &          (FRACREFB(IG,JPL+1)-FRACREFB(IG,JPL))
 3000    CONTINUE

C EMPIRICAL MODIFICATION TO CODE TO IMPROVING STRATOSPHERIC COOLING RATES
C FOR CO2.

         TAUG(LAY,8)=TAUG(LAY,8)*0.92
         TAUG(LAY,9)=TAUG(LAY,9)*0.88
         TAUG(LAY,10)=TAUG(LAY,10)*1.07
         TAUG(LAY,11)=TAUG(LAY,11)*1.1
         TAUG(LAY,12)=TAUG(LAY,12)*0.99
         TAUG(LAY,13)=TAUG(LAY,13)*0.88
         TAUG(LAY,14)=TAUG(LAY,14)*0.83

 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB5

C     BAND 5:  700-820 cm-1 (low key - H2O,CO2; low minor - O3, CCL4)
C                           (high key - O3,CO2)

      PARAMETER (MG=16, MXLAY=203, MAXXSEC=4, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG) 

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY),INDSELF(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K5/       KA(9,5,13,MG), KB(5,5,13:59,MG),
     &                  FORREF(4,MG),SELFREF(10,MG), KA_MO3(9,19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(585,MG),ABSB(1175,MG)
      DIMENSION FRACREFA(MG,9), FRACREFB(MG,5), CCL4(MG)


C Planck fraction mapping level : P = 473.42 mb, T = 259.83
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.4111E-01,1.4222E-01,1.3802E-01,1.3101E-01,1.2244E-01,1.0691E-01,
     &8.8703E-02,6.7130E-02,4.5509E-02,4.9866E-03,4.1214E-03,3.2557E-03,
     &2.3805E-03,1.5450E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.4152E-01,1.4271E-01,1.3784E-01,1.3075E-01,1.2215E-01,1.0674E-01,
     &8.8686E-02,6.7135E-02,4.5508E-02,4.9866E-03,4.1214E-03,3.2558E-03,
     &2.3805E-03,1.5450E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.4159E-01,1.4300E-01,1.3781E-01,1.3094E-01,1.2192E-01,1.0661E-01,
     &8.8529E-02,6.7127E-02,4.5511E-02,4.9877E-03,4.1214E-03,3.2558E-03,
     &2.3805E-03,1.5450E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.4162E-01,1.4337E-01,1.3774E-01,1.3122E-01,1.2172E-01,1.0641E-01,
     &8.8384E-02,6.7056E-02,4.5514E-02,4.9880E-03,4.1214E-03,3.2557E-03,
     &2.3805E-03,1.5450E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.4161E-01,1.4370E-01,1.3770E-01,1.3143E-01,1.2173E-01,1.0613E-01,
     &8.8357E-02,6.6874E-02,4.5509E-02,4.9883E-03,4.1214E-03,3.2558E-03,
     &2.3804E-03,1.5450E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.4154E-01,1.4405E-01,1.3771E-01,1.3169E-01,1.2166E-01,1.0603E-01,
     &8.8193E-02,6.6705E-02,4.5469E-02,4.9902E-03,4.1214E-03,3.2558E-03,
     &2.3804E-03,1.5450E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.4126E-01,1.4440E-01,1.3790E-01,1.3214E-01,1.2153E-01,1.0603E-01,
     &8.7908E-02,6.6612E-02,4.5269E-02,4.9900E-03,4.1256E-03,3.2558E-03,
     &2.3804E-03,1.5451E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.4076E-01,1.4415E-01,1.3885E-01,1.3286E-01,1.2147E-01,1.0612E-01,
     &8.7579E-02,6.6280E-02,4.4977E-02,4.9782E-03,4.1200E-03,3.2620E-03,
     &2.3820E-03,1.5452E-03,5.8423E-04,8.2275E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.4205E-01,1.4496E-01,1.4337E-01,1.3504E-01,1.2260E-01,1.0428E-01,
     &8.4946E-02,6.3625E-02,4.2951E-02,4.7313E-03,3.9157E-03,3.0879E-03,
     &2.2666E-03,1.5193E-03,5.7469E-04,8.1674E-05/

C Planck fraction mapping level : P = 0.2369280 mbar, T = 253.60 K
      DATA (FRACREFB(IG, 1),IG=1,16) /
     &1.4075E-01,1.4196E-01,1.3833E-01,1.3345E-01,1.2234E-01,1.0718E-01,
     &8.8004E-02,6.6308E-02,4.5028E-02,4.9029E-03,4.0377E-03,3.1870E-03,
     &2.3503E-03,1.5146E-03,5.7165E-04,8.2371E-05/
      DATA (FRACREFB(IG, 2),IG=1,16) /
     &1.4081E-01,1.4225E-01,1.3890E-01,1.3410E-01,1.2254E-01,1.0680E-01,
     &8.7391E-02,6.5819E-02,4.4725E-02,4.9121E-03,4.0420E-03,3.1869E-03,
     &2.3504E-03,1.5146E-03,5.7165E-04,8.2371E-05/
      DATA (FRACREFB(IG, 3),IG=1,16) /
     &1.4087E-01,1.4227E-01,1.3920E-01,1.3395E-01,1.2270E-01,1.0694E-01,
     &8.7229E-02,6.5653E-02,4.4554E-02,4.8797E-03,4.0460E-03,3.1939E-03,
     &2.3505E-03,1.5146E-03,5.7165E-04,8.1910E-05/
      DATA (FRACREFB(IG, 4),IG=1,16) /
     &1.4089E-01,1.4238E-01,1.3956E-01,1.3379E-01,1.2284E-01,1.0688E-01,
     &8.7192E-02,6.5490E-02,4.4390E-02,4.8395E-03,4.0173E-03,3.2070E-03,
     &2.3559E-03,1.5146E-03,5.7165E-04,8.2371E-05/
      DATA (FRACREFB(IG, 5),IG=1,16) /
     &1.4091E-01,1.4417E-01,1.4194E-01,1.3457E-01,1.2167E-01,1.0551E-01,
     &8.6450E-02,6.4889E-02,4.3584E-02,4.7551E-03,3.9509E-03,3.1374E-03,
     &2.3226E-03,1.4942E-03,5.7545E-04,8.0887E-05/

C Minor gas mapping level :
C     LOWER - O3, P = 317.34 mbar, T = 240.77 K
C     LOWER - CCL4

      DATA CCL4/
     &     26.1407,  53.9776,  63.8085,  36.1701,
     &     15.4099, 10.23116,  4.82948,  5.03836,
     &     1.75558,0.,0.,0.,
     &     0.,0.,0.,0./

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB
      REAL KA_MO3, MINORFRAC
      REAL O3M1, O3M2


C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 473.420 mb
      REFRAT_PLANCK_A = CHI_MLS(1,5)/CHI_MLS(2,5)

C     P = 0.2369 mb
      REFRAT_PLANCK_B = CHI_MLS(3,43)/CHI_MLS(2,43)

C     P = 317.3480
      REFRAT_M_A = CHI_MLS(1,7)/CHI_MLS(2,7)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the 
C     water vapor self-continuum and foreign continuum is 
C     interpolated (in temperature) separately.

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MO3 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MO3 = COLH2O(LAY)/SPECCOMB_MO3
         IF (SPECPARM_MO3 .GE. ONEMINUS) SPECPARM_MO3 = ONEMINUS
         SPECMULT_MO3 = 8.*SPECPARM_MO3
         JMO3 = 1 + INT(SPECMULT_MO3)
         FMO3 = AMOD(SPECMULT_MO3,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(5) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(5) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(5)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))

               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + WX(1,LAY) * CCL4(IG)
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(5)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF+ TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + WX(1,LAY) * CCL4(IG)
                FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(5)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + WX(1,LAY) * CCL4(IG)
            FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         SPECCOMB = COLO3(LAY) + RAT_O3CO2(LAY)*COLCO2(LAY)
         SPECPARM = COLO3(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLO3(LAY) + RAT_O3CO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLO3(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 4.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS1) * FAC01(LAY)
         FAC011 = (1. - FS1) * FAC11(LAY)
         FAC101 = FS1 * FAC01(LAY)
         FAC111 = FS1 * FAC11(LAY)

         SPECCOMB_PLANCK = COLO3(LAY)+REFRAT_PLANCK_B*COLCO2(LAY)
         SPECPARM_PLANCK = COLO3(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 4.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(5) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(5) + JS1
         
         DO 3000 IG = 1, NG(5)
            TAUG(LAY,IG) =  SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &          FAC100 * ABSB(IND0+1,IG) +
     &          FAC010 * ABSB(IND0+5,IG) +
     &          FAC110 * ABSB(IND0+6,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSB(IND1,IG) + 
     &          FAC101 * ABSB(IND1+1,IG) +
     &          FAC011 * ABSB(IND1+5,IG) +
     &          FAC111 * ABSB(IND1+6,IG)) 
     &          + WX(1,LAY) * CCL4(IG)
            FRACS(LAY,IG) = FRACREFB(IG,JPL) + FPL *
     &          (FRACREFB(IG,JPL+1)-FRACREFB(IG,JPL))
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB6

C     BAND 6:  820-980 cm-1 (low key - H2O; low minor - CO2)
C                           (high key - nothing; high minor - CFC11, CFC12)

      PARAMETER (MG=16, MXLAY=203, MAXXSEC=4, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K6/       KA(5,13,MG), FORREF(4,MG), SELFREF(10,MG), 
     &                  KA_MCO2(19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(65,MG)
      DIMENSION FRACREFA(MG),CFC11ADJ(MG), CFC12(MG)
      REAL KA_MCO2, MINORFRAC
      
C Planck fraction mapping level : P = 473.4280 mb, T = 259.83 K
      DATA FRACREFA /
     &1.4353E-01,1.4774E-01,1.4467E-01,1.3785E-01,1.2376E-01,1.0214E-01,
     &8.1984E-02,6.1152E-02,4.0987E-02,4.5067E-03,4.0020E-03,3.1772E-03,
     &2.3458E-03,1.5025E-03,5.7415E-04,8.2970E-05/

C Minor gas mapping level:
C     LOWER - CO2, P = 706.2720 mb, T = 294.2 K
C     UPPER - CFC11, CFC12

C      DATA CFC11/
C     &     0., 0., 26.5435, 108.850,
C     &     58.7804, 54.0875, 41.1065, 35.6120,
C     &     41.2328, 47.7402, 79.1026, 64.3005,
C     &     108.206, 141.617, 186.565, 58.4782/
C     CFC11 is multiplied by 1.385 to account for the 1060-1107 cm-1 band.

      DATA CFC11ADJ/
     &     0.,  0., 36.7627,    150.757,    
     &     81.4109, 74.9112, 56.9325, 49.3226,  
     &     57.1074, 66.1202, 109.557, 89.0562,  
     &     149.865, 196.140, 258.393, 80.9923/   
      DATA CFC12/
     &     62.8368, 43.2626, 26.7549, 22.2487,
     &     23.5029, 34.8323, 26.2335, 23.2306,
     &     18.4062, 13.9534, 22.6268, 24.2604,
     &     30.0088, 26.3634, 15.8237, 57.5050/


      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and
C     temperature. The water vapor self-continuum and foreign continuum
C     is interpolated (in temperature) separately.  

      HVRTAU = '$Revision: 3.3 $'
 
      DO 2500 LAY = 1, LAYTROP

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.77
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(6) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(6) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         DO 2000 IG = 1, NG(6)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG)))
            ABSCO2 =  (KA_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MCO2(INDM+1,IG) - KA_MCO2(INDM,IG)))
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR
     &           + ADJCOLCO2 * ABSCO2
     &           + WX(2,LAY) * CFC11ADJ(IG)
     &           + WX(3,LAY) * CFC12(IG)
            FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

C     Nothing important goes on above LAYTROP in this band.
      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG(6)
            TAUG(LAY,IG) = 0.0 
     &           + WX(2,LAY) * CFC11ADJ(IG)
     &           + WX(3,LAY) * CFC12(IG)
            FRACS(LAY,IG) = FRACREFA(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB7

C     BAND 7:  980-1080 cm-1 (low key - H2O,O3; low minor - CO2)
C                            (high key - O3; high minor - CO2)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K7/       KA(9,5,13,MG), KB(5,13:59,MG) , FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MCO2(9,19,MG), KB_MCO2(19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(585,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG,9),FRACREFB(MG)
      REAL KA_MCO2, KB_MCO2, MINORFRAC

C Planck fraction mapping level : P = 706.27 mb, T = 278.94 K
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.6312E-01,1.4949E-01,1.4305E-01,1.3161E-01,1.1684E-01,9.9900E-02,
     &8.0912E-02,6.0203E-02,4.0149E-02,4.3365E-03,3.5844E-03,2.8019E-03,
     &2.0756E-03,1.3449E-03,5.0492E-04,7.1194E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.6329E-01,1.4989E-01,1.4328E-01,1.3101E-01,1.1691E-01,9.9754E-02,
     &8.0956E-02,5.9912E-02,4.0271E-02,4.3298E-03,3.5626E-03,2.8421E-03,
     &2.1031E-03,1.3360E-03,4.8965E-04,6.8900E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.6236E-01,1.5081E-01,1.4341E-01,1.3083E-01,1.1684E-01,9.9701E-02,
     &8.0956E-02,5.9884E-02,4.0245E-02,4.3837E-03,3.6683E-03,2.9250E-03,
     &2.0969E-03,1.3320E-03,4.8965E-04,6.8900E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.6096E-01,1.5183E-01,1.4354E-01,1.3081E-01,1.1687E-01,9.9619E-02,
     &8.0947E-02,5.9899E-02,4.0416E-02,4.4389E-03,3.7280E-03,2.9548E-03,
     &2.0977E-03,1.3305E-03,4.8965E-04,6.8900E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.5661E-01,1.5478E-01,1.4414E-01,1.3097E-01,1.1695E-01,9.9823E-02,
     &8.0750E-02,6.0100E-02,4.0741E-02,4.4598E-03,3.7366E-03,2.9521E-03,
     &2.0980E-03,1.3297E-03,4.8965E-04,6.8900E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.4879E-01,1.5853E-01,1.4586E-01,1.3162E-01,1.1729E-01,1.0031E-01,
     &8.0908E-02,6.0460E-02,4.1100E-02,4.4578E-03,3.7388E-03,2.9508E-03,
     &2.0986E-03,1.3288E-03,4.8965E-04,6.8900E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.4117E-01,1.4838E-01,1.4807E-01,1.3759E-01,1.2218E-01,1.0228E-01,
     &8.2130E-02,6.1546E-02,4.1522E-02,4.4577E-03,3.7428E-03,2.9475E-03,
     &2.0997E-03,1.3277E-03,4.8965E-04,6.8900E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.4018E-01,1.4207E-01,1.3919E-01,1.3332E-01,1.2325E-01,1.0915E-01,
     &9.0280E-02,6.5554E-02,4.1852E-02,4.4707E-03,3.7572E-03,2.9364E-03,
     &2.1023E-03,1.3249E-03,4.8965E-04,6.8900E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.4863E-01,1.4926E-01,1.4740E-01,1.3558E-01,1.1999E-01,1.0044E-01,
     &8.1927E-02,6.0989E-02,4.0665E-02,4.4481E-03,3.7369E-03,2.9482E-03,
     &2.0976E-03,1.3281E-03,4.8965E-04,6.8900E-05/

c Planck fraction mapping level : P=95.58 mbar, T= 215.70 K
      DATA FRACREFB /
     &1.5872E-01,1.5443E-01,1.4413E-01,1.3147E-01,1.1634E-01,9.8914E-02,
     &8.0236E-02,6.0197E-02,4.0624E-02,4.4225E-03,3.6688E-03,2.9074E-03,
     &2.0862E-03,1.3039E-03,4.8561E-04,6.8854E-05/

C Minor gas mapping level :
C     LOWER - CO2, P = 706.2620 mbar, T= 278.94 K
C     UPPER - CO2, P = 12.9350 mbar, T = 234.01 K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.

C     P = 706.2620 mb
      REFRAT_PLANCK_A = CHI_MLS(1,3)/CHI_MLS(3,3)

C     P = 706.2720 mb
      REFRAT_M_A = CHI_MLS(1,3)/CHI_MLS(3,3)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately. 

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OO3(LAY)*COLO3(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OO3_1(LAY)*COLO3(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MCO2 = COLH2O(LAY) + REFRAT_M_A*COLO3(LAY)
         SPECPARM_MCO2 = COLH2O(LAY)/SPECCOMB_MCO2
         IF (SPECPARM_MCO2 .GE. ONEMINUS) SPECPARM_MCO2 = ONEMINUS
         SPECMULT_MCO2 = 8.*SPECPARM_MCO2

         JMCO2 = 1 + INT(SPECMULT_MCO2)
         FMCO2 = AMOD(SPECMULT_MCO2,1.0)

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 3.0+(RATCO2-3.0)**0.79
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLO3(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(7) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(7) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(7)

               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY) * (CO2M2 - CO2M1)
               TAUG(LAY,IG) = SPECCOMB *
     &          (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC200 * ABSA(IND0+2,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG) +
     &          FAC210 * ABSA(IND0+11,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC201 * ABSA(IND1+2,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG) +
     &          FAC211 * ABSA(IND1+11,IG)) 
     &          + TAUSELF + TAUFOR
     &          + ADJCOLCO2*ABSCO2
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY) * (CO2M2 - CO2M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(7)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY) * (CO2M2 - CO2M1)
               TAUG(LAY,IG) = SPECCOMB * 
     &         (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG)) 
     &          + TAUSELF + TAUFOR
     &          + ADJCOLCO2*ABSCO2
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020    CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.79
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(7) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(7) + 1
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG(7)
            ABSCO2 = KB_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MCO2(INDM+1,IG) - KB_MCO2(INDM,IG))
            TAUG(LAY,IG) = COLO3(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
     &           + ADJCOLCO2 * ABSCO2
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE

C EMPIRICAL MODIFICATION TO CODE TO IMPROVE STRATOSPHERIC COOLING RATES
C FOR O3 

         TAUG(LAY,8)=TAUG(LAY,8)*0.92
         TAUG(LAY,9)=TAUG(LAY,9)*0.88
         TAUG(LAY,10)=TAUG(LAY,10)*1.07
         TAUG(LAY,11)=TAUG(LAY,11)*1.1
         TAUG(LAY,12)=TAUG(LAY,12)*0.99
         TAUG(LAY,13)=TAUG(LAY,13)*0.88
         TAUG(LAY,14)=TAUG(LAY,14)*0.83

 3500 CONTINUE

      RETURN
      END

*******************************************************************************

      SUBROUTINE TAUGB8

C     BAND 8:  1080-1180 cm-1 (low key - H2O; low minor - CO2,O3,N2O)
C                             (high key - O3; high minor - CO2, N2O)

      PARAMETER (MG=16, MXLAY=203, MAXXSEC=4, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY),SELFFRAC(MXLAY),INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K8/       KA(5,13,MG), KB(5,13:59,MG), FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MCO2(19,MG), KA_MO3(19,MG),
     &                  KA_MN2O(19,MG), KB_MCO2(19,MG), KB_MN2O(19,MG)


      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      REAL KA,KB,KA_MCO2,KA_MO3,KA_MN2O,KB_MCO2,KB_MN2O, MINORFRAC

      DIMENSION ABSA(65,MG),ABSB(235,MG),CFC12(MG),CFC22ADJ(MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

C Planck fraction mapping level : P=473.4280 mb, T = 259.83 K
      DATA FRACREFA /
     &1.6004E-01,1.5437E-01,1.4502E-01,1.3084E-01,1.1523E-01,9.7743E-02,
     &8.0376E-02,6.0261E-02,4.1111E-02,4.4772E-03,3.6511E-03,2.9154E-03,
     &2.1184E-03,1.3048E-03,4.6637E-04,6.5624E-05/

C Planck fraction mapping level : P=95.5835 mb, T= 215.7 K
      DATA FRACREFB /
     &1.4987E-01,1.4665E-01,1.4154E-01,1.3200E-01,1.1902E-01,1.0352E-01,
     &8.4939E-02,6.4105E-02,4.3190E-02,4.5129E-03,3.7656E-03,2.8733E-03,
     &2.0947E-03,1.3201E-03,5.1832E-04,7.7473E-05/

C Minor gas mapping level:
C     LOWER - CO2, P = 1053.63 mb, T = 294.2 K
C     LOWER - O3,  P = 317.348 mb, T = 240.77 K
C     LOWER - N2O, P = 706.2720 mb, T= 278.94 K
C     LOWER - CFC12,CFC11
C     UPPER - CO2, P = 35.1632 mb, T = 223.28 K
C     UPPER - N2O, P = 8.716e-2 mb, T = 226.03 K

      DATA CFC12/
     &     85.4027, 89.4696, 74.0959, 67.7480,
     &     61.2444, 59.9073, 60.8296, 63.0998,
     &     59.6110, 64.0735, 57.2622, 58.9721,
     &     43.5505, 26.1192, 32.7023, 32.8667/
C     Original CFC22 is multiplied by 1.485 to account for the 780-850 cm-1 
C     and 1290-1335 cm-1 bands.
      DATA CFC22ADJ/
     &     135.335, 89.6642, 76.2375, 65.9748,
     &     63.1164, 60.2935, 64.0299, 75.4264,
     &     51.3018, 7.07911, 5.86928, 0.398693,
     &     2.82885, 9.12751, 6.28271, 0./

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water vapor 
C     self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.65
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(8) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(8) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         DO 2000 IG = 1, NG(8)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG)))
            ABSCO2 =  (KA_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MCO2(INDM+1,IG) - KA_MCO2(INDM,IG)))
            ABSO3 =  (KA_MO3(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MO3(INDM+1,IG) - KA_MO3(INDM,IG)))
            ABSN2O =  (KA_MN2O(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MN2O(INDM+1,IG) - KA_MN2O(INDM,IG)))
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR
     &           + ADJCOLCO2*ABSCO2
     &           + COLO3(LAY) * ABSO3
     &           + COLN2O(LAY) * ABSN2O
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)
            FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/COLDRY(LAY)
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.65
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           * COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(8) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(8) + 1
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG(8)
            ABSCO2 =  (KB_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MCO2(INDM+1,IG) - KB_MCO2(INDM,IG)))
            ABSN2O =  (KB_MN2O(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MN2O(INDM+1,IG) - KB_MN2O(INDM,IG)))
            TAUG(LAY,IG) = COLO3(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
     &           + ADJCOLCO2*ABSCO2
     &           + COLN2O(LAY)*ABSN2O 
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB9

C     BAND 9:  1180-1390 cm-1 (low key - H2O,CH4; low minor - N2O)
C                             (high key - CH4; high minor - N2O)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K9/       KA(9,5,13,MG),KB(5,13:59,MG),FORREF(4,MG),
     &                  SELFREF(10,MG),KA_MN2O(9,19,MG),KB_MN2O(19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      REAL KA,KB
      REAL KA_MN2O,KB_MN2O,MINORFRAC,N2OM1,N2OM2
      DIMENSION ABSA(585,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG,9), FRACREFB(MG)

C Planck fractions mapping level : P=212.7250 mb, T = 223.06 K

      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.8129E-01,1.6119E-01,1.3308E-01,1.2342E-01,1.1259E-01,9.7580E-02,
     &7.9176E-02,5.8541E-02,3.9084E-02,4.2419E-03,3.4314E-03,2.6935E-03,
     &1.9404E-03,1.2218E-03,4.5263E-04,6.0909E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.9665E-01,1.5640E-01,1.3101E-01,1.2153E-01,1.1037E-01,9.6043E-02,
     &7.7856E-02,5.7547E-02,3.8670E-02,4.1955E-03,3.4104E-03,2.6781E-03,
     &1.9245E-03,1.2093E-03,4.4113E-04,6.0913E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &2.0273E-01,1.5506E-01,1.3044E-01,1.2043E-01,1.0952E-01,9.5384E-02,
     &7.7157E-02,5.7176E-02,3.8379E-02,4.1584E-03,3.3836E-03,2.6412E-03,
     &1.8865E-03,1.1791E-03,4.2094E-04,4.7410E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &2.0272E-01,1.5963E-01,1.2913E-01,1.2060E-01,1.0820E-01,9.4685E-02,
     &7.6544E-02,5.6851E-02,3.8155E-02,4.0913E-03,3.3442E-03,2.6054E-03,
     &1.8875E-03,1.1263E-03,3.7743E-04,4.7410E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &2.0280E-01,1.6353E-01,1.2910E-01,1.1968E-01,1.0725E-01,9.4112E-02,
     &7.5828E-02,5.6526E-02,3.7972E-02,4.0205E-03,3.3063E-03,2.5681E-03,
     &1.8386E-03,1.0757E-03,3.5301E-04,4.7410E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &2.0294E-01,1.6840E-01,1.2852E-01,1.1813E-01,1.0724E-01,9.2946E-02,
     &7.5029E-02,5.6158E-02,3.7744E-02,3.9632E-03,3.2434E-03,2.5275E-03,
     &1.7558E-03,1.0080E-03,3.5301E-04,4.7410E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &2.0313E-01,1.7390E-01,1.2864E-01,1.1689E-01,1.0601E-01,9.1791E-02,
     &7.4224E-02,5.5500E-02,3.7374E-02,3.9214E-03,3.1984E-03,2.4162E-03,
     &1.6394E-03,9.7275E-04,3.5299E-04,4.7410E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &2.0332E-01,1.7800E-01,1.3286E-01,1.1555E-01,1.0407E-01,9.0475E-02,
     &7.2452E-02,5.4566E-02,3.6677E-02,3.7889E-03,3.0351E-03,2.2587E-03,
     &1.5764E-03,9.7270E-04,3.5300E-04,4.7410E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.9624E-01,1.6519E-01,1.3663E-01,1.1535E-01,1.0719E-01,9.4156E-02,
     &7.6745E-02,5.6987E-02,3.8135E-02,4.1626E-03,3.4243E-03,2.7116E-03,
     &1.7095E-03,9.7271E-04,3.5299E-04,4.7410E-05/

C Planck fraction mapping level : P=3.20e-2 mb, T = 197.92 K

      DATA FRACREFB /
     &2.0914E-01,1.5077E-01,1.2878E-01,1.1856E-01,1.0695E-01,9.3048E-02,
     &7.7645E-02,6.0785E-02,4.0642E-02,4.0499E-03,3.3931E-03,2.6363E-03,
     &1.9151E-03,1.1963E-03,4.3471E-04,5.1421E-05/

C Minor gas mapping level :
C     LOWER - N2O, P = 706.272 mbar, T = 278.94 K
C     UPPER - N2O, P = 95.58 mbar, T = 215.7 K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 212 mb
      REFRAT_PLANCK_A = CHI_MLS(1,9)/CHI_MLS(6,9)

C     P = 706.272 mb 
      REFRAT_M_A = CHI_MLS(1,3)/CHI_MLS(6,3)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCH4(LAY)*COLCH4(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCH4_1(LAY)*COLCH4(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_A*COLCH4(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)

c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/(COLDRY(LAY))
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCH4(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(9) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(9) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(9)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               TAUG(LAY,IG) = SPECCOMB *
     &          (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC200 * ABSA(IND0+2,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG) +
     &          FAC210 * ABSA(IND0+11,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC201 * ABSA(IND1+2,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG) +
     &          FAC211 * ABSA(IND1+11,IG)) 
     &          + TAUSELF + TAUFOR
     &          + ADJCOLN2O*ABSN2O            
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(9)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(9)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/(COLDRY(LAY))
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(9) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(9) + 1
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG(9)
            ABSN2O = KB_MN2O(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MN2O(INDM+1,IG) - KB_MN2O(INDM,IG))
            TAUG(LAY,IG) = COLCH4(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
     &           + ADJCOLN2O*ABSN2O
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE
      RETURN
      END

*******************************************************************************


      SUBROUTINE TAUGB10

C     BAND 10:  1390-1480 cm-1 (low key - H2O; high key - H2O)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K10/      KA(5,13,MG), KB(5,13:59,MG), FORREF(4,MG),
     &                  SELFREF(10,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

C Planck fraction mapping level : P = 212.7250, T = 223.06 K
      DATA FRACREFA/
     & 1.6909E-01, 1.5419E-01, 1.3999E-01, 1.2637E-01,
     & 1.1429E-01, 9.9676E-02, 8.0093E-02, 6.0283E-02,
     & 4.1077E-02, 4.4857E-03, 3.6545E-03, 2.9243E-03,
     & 2.0407E-03, 1.2891E-03, 4.8767E-04, 6.7748E-05/
C Planck fraction mapping level : P = 95.58350 mb, T = 215.70 K
      DATA FRACREFB/
     & 1.7391E-01, 1.5680E-01, 1.4419E-01, 1.2672E-01,
     & 1.0708E-01, 9.7034E-02, 7.8545E-02, 5.9784E-02,
     & 4.0879E-02, 4.4704E-03, 3.7150E-03, 2.9038E-03,
     & 2.1454E-03, 1.2802E-03, 4.8328E-04, 6.7378E-05/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(10) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(10) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         DO 2000 IG = 1, NG(10)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR
            FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(10) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(10) + 1
         INDF = INDFOR(LAY)
         DO 3000 IG = 1, NG(10)
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
     &           + TAUFOR
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

*******************************************************************************

      SUBROUTINE TAUGB11

C     BAND 11:  1480-1800 cm-1 (low - H2O; low minor - O2)
C                              (high key - H2O; high minor - O2)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K11/      KA(5,13,MG), KB(5,13:59,MG) , FORREF(4,MG),
     &                  SELFREF(10,MG), KA_MO2(19,MG), KB_MO2(19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)
      REAL KA_MO2, KB_MO2, MINORFRAC

C Planck fraction mapping level : P=1053.63 mb, T= 294.2 K
      DATA FRACREFA /
     &1.4601E-01,1.3824E-01,1.4240E-01,1.3463E-01,1.1948E-01,1.0440E-01,
     &8.8667E-02,6.5792E-02,4.3893E-02,4.7941E-03,4.0760E-03,3.3207E-03,
     &2.4087E-03,1.3912E-03,4.3482E-04,6.0932E-05/
C Planck fraction mapping level : P=0.353 mb, T = 262.11 K
      DATA FRACREFB /
     &7.2928E-02,1.4900E-01,1.6156E-01,1.5603E-01,1.3934E-01,1.1394E-01,
     &8.8783E-02,6.2411E-02,4.0191E-02,4.4587E-03,3.9533E-03,3.0847E-03,
     &2.2317E-03,1.4410E-03,5.6722E-04,7.7933E-05/

C Minor gas mapping level :
C     LOWER - O2, P = 706.2720 mbar, T = 278.94 K
C     UPPER - O2, P = 4.758820 mbarm T = 250.85 K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(11) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(11) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         SCALEO2 = COLO2(LAY)*SCALEMINOR(LAY)
         DO 2000 IG = 1, NG(11)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG)))
            TAUO2 =  SCALEO2 *
     &           (KA_MO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MO2(INDM+1,IG) - KA_MO2(INDM,IG)))
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG))
     &           + TAUSELF + TAUFOR
     &           + TAUO2
            FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE
      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(11) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(11) + 1
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         SCALEO2 = COLO2(LAY)*SCALEMINOR(LAY)
         DO 3000 IG = 1, NG(11)
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUO2 =  SCALEO2*
     &           (KB_MO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MO2(INDM+1,IG) - KB_MO2(INDM,IG)))
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
     &           + TAUFOR
     &           + TAUO2
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB12

C     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K12/      KA(9,5,13,MG),FORREF(4,MG),SELFREF(10,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(585,MG)
      DIMENSION FRACREFA(MG,9)

C P = 174.1640 mbar, T= 215.78 K
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.3984E-01,1.6809E-01,1.8072E-01,1.5400E-01,1.2613E-01,9.6959E-02,
     &5.9713E-02,3.8631E-02,2.6937E-02,3.1711E-03,2.3458E-03,1.4653E-03,
     &1.0567E-03,6.6504E-04,2.4957E-04,3.5172E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.2745E-01,1.6107E-01,1.6568E-01,1.5436E-01,1.3183E-01,1.0166E-01,
     &6.4506E-02,4.7756E-02,3.4472E-02,3.7189E-03,2.9349E-03,2.1469E-03,
     &1.3746E-03,7.1691E-04,2.8057E-04,5.6242E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.2181E-01,1.5404E-01,1.6540E-01,1.5255E-01,1.3736E-01,9.8856E-02,
     &6.8927E-02,5.1385E-02,3.7046E-02,4.0302E-03,3.0949E-03,2.3772E-03,
     &1.6538E-03,8.9641E-04,4.6991E-04,1.1251E-04/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.1794E-01,1.4864E-01,1.6316E-01,1.5341E-01,1.3986E-01,9.6656E-02,
     &7.2478E-02,5.5061E-02,3.8886E-02,4.3398E-03,3.3576E-03,2.4891E-03,
     &1.7674E-03,1.0764E-03,7.7689E-04,1.1251E-04/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.1635E-01,1.4342E-01,1.5924E-01,1.5670E-01,1.3740E-01,9.7087E-02,
     &7.6250E-02,5.7802E-02,4.0808E-02,4.4113E-03,3.6035E-03,2.6269E-03,
     &1.7586E-03,1.6498E-03,7.7689E-04,1.1251E-04/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.1497E-01,1.3751E-01,1.5587E-01,1.5904E-01,1.3140E-01,1.0159E-01,
     &7.9729E-02,6.1475E-02,4.2382E-02,4.5291E-03,3.8161E-03,2.7683E-03,
     &1.9899E-03,2.0395E-03,7.7720E-04,1.1251E-04/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.1331E-01,1.3015E-01,1.5574E-01,1.5489E-01,1.2697E-01,1.0746E-01,
     &8.4777E-02,6.5145E-02,4.4293E-02,4.7426E-03,3.8383E-03,2.9065E-03,
     &2.8430E-03,2.0401E-03,7.7689E-04,1.1251E-04/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.0993E-01,1.2320E-01,1.4893E-01,1.4573E-01,1.3174E-01,1.1149E-01,
     &9.3326E-02,6.9942E-02,4.6762E-02,4.9309E-03,3.8583E-03,4.1889E-03,
     &3.0415E-03,2.0406E-03,7.7720E-04,1.1251E-04/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.2028E-01,1.2091E-01,1.3098E-01,1.3442E-01,1.3574E-01,1.1739E-01,
     &9.5343E-02,7.0224E-02,5.3456E-02,6.0206E-03,5.0758E-03,4.1906E-03,
     &3.0431E-03,2.0400E-03,7.7689E-04,1.1251E-04/

      EQUIVALENCE (KA,ABSA)
      REAL KA


C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P =   174.164 mb 
      REFRAT_PLANCK_A = CHI_MLS(1,10)/CHI_MLS(2,10)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum adn foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(12) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(12) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(12)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG(12)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(12)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG(12)
            TAUG(LAY,IG) = 0.0
            FRACS(LAY,IG) = 0.0
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB13

C     BAND 13:  2080-2250 cm-1 (low key - H2O,N2O; high minor - O3 minor)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K13/      KA(9,5,13,MG),FORREF(4,MG),SELFREF(10,MG),
     &                  KA_MCO2(9,19,MG), KA_MCO(9,19,MG),KB_MO3(19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(585,MG)
      DIMENSION FRACREFA(MG,9),FRACREFB(MG)
      REAL KA_MCO2,KA_MCO,KB_MO3,MINORFRAC

C Planck fraction mapping level : P=473.4280 mb, T = 259.83 K      
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.7534E-01,1.7394E-01,1.6089E-01,1.3782E-01,1.0696E-01,8.5853E-02,
     &6.6548E-02,4.9053E-02,3.2064E-02,3.4820E-03,2.8763E-03,2.2204E-03,
     &1.5612E-03,9.8572E-04,3.6853E-04,5.1612E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.7489E-01,1.7309E-01,1.5981E-01,1.3782E-01,1.0797E-01,8.6367E-02,
     &6.7042E-02,4.9257E-02,3.2207E-02,3.4820E-03,2.8767E-03,2.2203E-03,
     &1.5613E-03,9.8571E-04,3.6853E-04,5.1612E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.7459E-01,1.7259E-01,1.5948E-01,1.3694E-01,1.0815E-01,8.7376E-02,
     &6.7339E-02,4.9541E-02,3.2333E-02,3.5019E-03,2.8958E-03,2.2527E-03,
     &1.6099E-03,9.8574E-04,3.6853E-04,5.1612E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.7391E-01,1.7244E-01,1.5921E-01,1.3644E-01,1.0787E-01,8.7776E-02,
     &6.8361E-02,4.9628E-02,3.2578E-02,3.5117E-03,2.9064E-03,2.2571E-03,
     &1.6887E-03,1.0045E-03,3.6853E-04,5.1612E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.7338E-01,1.7157E-01,1.5957E-01,1.3571E-01,1.0773E-01,8.7966E-02,
     &6.9000E-02,5.0300E-02,3.2813E-02,3.5470E-03,2.9425E-03,2.2552E-03,
     &1.7038E-03,1.1025E-03,3.6853E-04,5.1612E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.7230E-01,1.7082E-01,1.5917E-01,1.3562E-01,1.0806E-01,8.7635E-02,
     &6.9815E-02,5.1155E-02,3.3139E-02,3.6264E-03,2.9436E-03,2.3417E-03,
     &1.7731E-03,1.1156E-03,4.4533E-04,5.1612E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.7073E-01,1.6961E-01,1.5844E-01,1.3594E-01,1.0821E-01,8.7791E-02,
     &7.0502E-02,5.1904E-02,3.4107E-02,3.5888E-03,2.9574E-03,2.5851E-03,
     &1.9127E-03,1.1537E-03,4.7789E-04,1.0016E-04/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.6700E-01,1.6848E-01,1.5628E-01,1.3448E-01,1.1011E-01,8.9016E-02,
     &7.1973E-02,5.2798E-02,3.5650E-02,3.8534E-03,3.4142E-03,2.7799E-03,
     &2.1288E-03,1.3043E-03,6.2858E-04,1.0016E-04/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.6338E-01,1.5565E-01,1.4470E-01,1.3500E-01,1.1909E-01,9.8312E-02,
     &7.9023E-02,5.5728E-02,3.6831E-02,3.6569E-03,3.0552E-03,2.3431E-03,
     &1.7088E-03,1.1082E-03,3.6829E-04,5.1612E-05/

C Planck fraction mapping level : P=4.758820 mb, T = 250.85 K
      DATA FRACREFB /
     &1.5411E-01,1.3573E-01,1.2527E-01,1.2698E-01,1.2394E-01,1.0876E-01,
     &8.9906E-02,6.9551E-02,4.8240E-02,5.2434E-03,4.3630E-03,3.4262E-03,
     &2.5124E-03,1.5479E-03,3.7294E-04,5.1050E-05/

C Minor gas mapping levels :
C     LOWER - CO2, P = 1053.63 mb, T = 294.2 K
C     LOWER - CO, P = 706 mb, T = 278.94 K
C     UPPER - O3, P = 95.5835 mb, T = 215.7 K

      EQUIVALENCE (KA,ABSA)
      REAL KA


C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 473.420 mb (Level 5)
      REFRAT_PLANCK_A = CHI_MLS(1,5)/CHI_MLS(4,5)

C     P = 1053. (Level 1)
      REFRAT_M_A = CHI_MLS(1,1)/CHI_MLS(4,1)

C     P = 706. (Level 3)
      REFRAT_M_A3 = CHI_MLS(1,3)/CHI_MLS(4,3)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2ON2O(LAY)*COLN2O(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2ON2O_1(LAY)*COLN2O(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MCO2 = COLH2O(LAY) + REFRAT_M_A*COLN2O(LAY)
         SPECPARM_MCO2 = COLH2O(LAY)/SPECCOMB_MCO2
         IF (SPECPARM_MCO2 .GE. ONEMINUS) SPECPARM_MCO2 = ONEMINUS
         SPECMULT_MCO2 = 8.*SPECPARM_MCO2
         JMCO2 = 1 + INT(SPECMULT_MCO2)
         FMCO2 = AMOD(SPECMULT_MCO2,1.0)

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/3.55E-4
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.68
            ADJCOLCO2 = ADJFAC*3.55E-4*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         SPECCOMB_MCO = COLH2O(LAY) + REFRAT_M_A3*COLN2O(LAY)
         SPECPARM_MCO = COLH2O(LAY)/SPECCOMB_MCO
         IF (SPECPARM_MCO .GE. ONEMINUS) SPECPARM_MCO = ONEMINUS
         SPECMULT_MCO = 8.*SPECPARM_MCO
         JMCO = 1 + INT(SPECMULT_MCO)
         FMCO = AMOD(SPECMULT_MCO,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLN2O(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(13) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(13) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(13)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + 
     &              MINORFRAC(LAY) * (CO2M2 - CO2M1)
               COM1 = KA_MCO(JMCO,INDM,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM,IG)-
     &              KA_MCO(JMCO,INDM,IG))
               COM2 = KA_MCO(JMCO,INDM+1,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM+1,IG)-
     &              KA_MCO(JMCO,INDM+1,IG))
               ABSCO = COM1 + 
     &              MINORFRAC(LAY) * (COM2 - COM1)
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))+
     &              SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
     &              + COLCO(LAY)*ABSCO
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(13)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 +
     &              MINORFRAC(LAY) * (CO2M2 - CO2M1)
               COM1 = KA_MCO(JMCO,INDM,IG) + FMCO*
     &           (KA_MCO(JMCO+1,INDM,IG)-
     &              KA_MCO(JMCO,INDM,IG))
               COM2 = KA_MCO(JMCO,INDM+1,IG) + FMCO*
     &           (KA_MCO(JMCO+1,INDM+1,IG)-
     &              KA_MCO(JMCO,INDM+1,IG))
               ABSCO = COM1 +
     &              MINORFRAC(LAY) * (COM2 - COM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
     &              + COLCO(LAY)*ABSCO
                FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(13)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + 
     &              MINORFRAC(LAY) * (CO2M2 - CO2M1)
               COM1 = KA_MCO(JMCO,INDM,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM,IG)-
     &              KA_MCO(JMCO,INDM,IG))
               COM2 = KA_MCO(JMCO,INDM+1,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM+1,IG)-
     &              KA_MCO(JMCO,INDM+1,IG))
               ABSCO = COM1 +
     &              MINORFRAC(LAY) * (COM2 - COM1)
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
     &              + COLCO(LAY)*ABSCO
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG(13)
            ABSO3 = KB_MO3(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MO3(INDM+1,IG) - KB_MO3(INDM,IG))
            TAUG(LAY,IG) = COLO3(LAY)*ABSO3
            FRACS(LAY,IG) =  FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

*******************************************************************************

      SUBROUTINE TAUGB14

C     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K14/      KA(5,13,MG), KB(5,13:59,MG) , FORREF(4,MG),
     &                  SELFREF(10,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

C Planck fraction mapping level : P = 142.5940 mb, T = 215.70 K
      DATA FRACREFA /
     & 1.9360E-01, 1.7276E-01, 1.4811E-01, 1.2238E-01, 
     & 1.0242E-01, 8.6830E-02, 7.1890E-02, 5.4030E-02, 
     & 3.5075E-02, 3.8052E-03, 3.1458E-03, 2.4873E-03,
     & 1.8182E-03, 1.1563E-03, 4.3251E-04, 5.7744E-05/
C Planck fraction mapping level : P = 4.758820mb, T = 250.85 K
      DATA FRACREFB /
     & 1.8599E-01, 1.6646E-01, 1.4264E-01, 1.2231E-01, 
     & 1.0603E-01, 9.2014E-02, 7.5287E-02, 5.6758E-02, 
     & 3.8386E-02, 4.2139E-03, 3.5399E-03, 2.7381E-03,
     & 1.9202E-03, 1.2083E-03, 4.5395E-04, 6.2699E-05/

      Equivalence (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum 
C     and foreign continuum is interpolated (in temperature) separately.  

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(14) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(14) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         DO 2000 IG = 1, NG(14)
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUG(LAY,IG) = COLCO2(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG))
     &           + TAUSELF + TAUFOR
            FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(14) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(14) + 1
         DO 3000 IG = 1, NG(14)
            TAUG(LAY,IG) = COLCO2(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB15

C     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; low minor - N2)
C                              (high - nothing)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K15/      KA(9,5,13,MG), 
     &                  FORREF(4,MG), SELFREF(10,MG), KA_MN2(9,19,MG)

      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(585,MG)
      DIMENSION FRACREFA(MG,9)

C Planck fraction mapping level : P = 1053. mb, T = 294.2 K
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.0689E-01,1.1563E-01,1.2447E-01,1.2921E-01,1.2840E-01,1.2113E-01,
     &1.0643E-01,8.4987E-02,6.0142E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &1.0782E-01,1.1637E-01,1.2290E-01,1.2911E-01,1.2841E-01,1.2113E-01,
     &1.0643E-01,8.4987E-02,6.0142E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &1.0858E-01,1.1860E-01,1.2237E-01,1.2665E-01,1.2841E-01,1.2111E-01,
     &1.0642E-01,8.4987E-02,6.0142E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &1.1022E-01,1.1965E-01,1.2334E-01,1.2383E-01,1.2761E-01,1.2109E-01,
     &1.0642E-01,8.4987E-02,6.0142E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &1.1342E-01,1.2069E-01,1.2360E-01,1.2447E-01,1.2340E-01,1.2020E-01,
     &1.0639E-01,8.4987E-02,6.0142E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &1.1771E-01,1.2280E-01,1.2177E-01,1.2672E-01,1.2398E-01,1.1787E-01,
     &1.0131E-01,8.4987E-02,6.0142E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &1.2320E-01,1.2491E-01,1.2001E-01,1.2936E-01,1.2653E-01,1.1929E-01,
     &9.8955E-02,7.4887E-02,6.0142E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &1.3105E-01,1.2563E-01,1.3055E-01,1.2854E-01,1.3402E-01,1.1571E-01,
     &9.4876E-02,6.0459E-02,5.6457E-02,6.6798E-03,5.5293E-03,4.3700E-03,
     &3.2061E-03,2.0476E-03,7.7366E-04,1.0897E-04/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &1.1375E-01,1.2090E-01,1.2348E-01,1.2458E-01,1.2406E-01,1.1921E-01,
     &1.0802E-01,8.6613E-02,5.8125E-02,6.2984E-03,5.2359E-03,4.0641E-03,
     &2.9379E-03,1.9001E-03,7.2646E-04,1.0553E-04/

      EQUIVALENCE (KA,ABSA)
      REAL KA, KA_MN2, MINORFRAC
      REAL N2M1,N2M2

C Minor gas mapping level : 
C     Lower - Nitrogen Continuum, P = 1053., T = 294.

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.
C     P = 1053. mb (Level 1)
      REFRAT_PLANCK_A = CHI_MLS(4,1)/CHI_MLS(2,1)

C     P = 1053.
      REFRAT_M_A = CHI_MLS(4,1)/CHI_MLS(2,1)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLN2O(LAY) + RAT_N2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLN2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLN2O(LAY) + RAT_N2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLN2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2 = COLN2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2 = COLN2O(LAY)/SPECCOMB_MN2
         IF (SPECPARM_MN2 .GE. ONEMINUS) SPECPARM_MN2 = ONEMINUS
         SPECMULT_MN2 = 8.*SPECPARM_MN2
         JMN2 = 1 + INT(SPECMULT_MN2)
         FMN2 = AMOD(SPECMULT_MN2,1.0)

         SPECCOMB_PLANCK = COLN2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLN2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(15) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(15) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         
         SCALEN2 = COLBRD(LAY)*SCALEMINOR(LAY)
         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(15)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG)-
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG)-
     &              KA_MN2(JMN2,INDM+1,IG))
               TAUN2 = SCALEN2*
     &              (N2M1 + MINORFRAC(LAY) * (N2M2 - N2M1))
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + TAUN2
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG(15)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG)-
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG)-
     &              KA_MN2(JMN2,INDM+1,IG))
               TAUN2 = SCALEN2*
     &              (N2M1 + MINORFRAC(LAY) * (N2M2 - N2M1))
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + TAUN2
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG(15)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG)-
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG)-
     &              KA_MN2(JMN2,INDM+1,IG))
               TAUN2 = SCALEN2 *
     &              (N2M1 + MINORFRAC(LAY) * (N2M2 - N2M1))
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + TAUN2
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG(15)
            TAUG(LAY,IG) = 0.0
            FRACS(LAY,IG) = 0.0
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------


      SUBROUTINE TAUGB16

C     BAND 16:  2600-3250 cm-1 (low key- H2O,CH4; high key - CH4)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)
      COMMON /PLANKG/   FRACS(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFI/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K16/      KA(9,5,13,MG),KB(5,13:59,MG),
     &                  FORREF(4,MG),SELFREF(10,MG)


      COMMON /CVRTAU/    HVRTAU

      CHARACTER*16       HVRTAU

      DIMENSION ABSA(585,MG), ABSB(235,MG)
      DIMENSION FRACREFA(MG,9), FRACREFB(MG)

C Planck fraction mapping level: P = 387.6100 mbar, T = 250.17 K
      DATA (FRACREFA(IG, 1),IG=1,16) /
     &1.1593E-01,2.3390E-01,1.9120E-01,1.3121E-01,1.0590E-01,8.4852E-02,
     &6.4168E-02,4.2537E-02,2.3220E-02,2.1767E-03,1.8203E-03,1.3724E-03,
     &9.5452E-04,5.5015E-04,1.9348E-04,2.7344E-05/
      DATA (FRACREFA(IG, 2),IG=1,16) /
     &2.8101E-01,1.9773E-01,1.4749E-01,1.1399E-01,8.8190E-02,7.0531E-02,
     &4.6356E-02,3.0774E-02,1.7332E-02,2.0054E-03,1.5950E-03,1.2760E-03,
     &9.5034E-04,5.4992E-04,1.9349E-04,2.7309E-05/
      DATA (FRACREFA(IG, 3),IG=1,16) /
     &2.9054E-01,2.1263E-01,1.4133E-01,1.1083E-01,8.5107E-02,6.5247E-02,
     &4.4542E-02,2.7205E-02,1.6495E-02,1.8453E-03,1.5222E-03,1.1884E-03,
     &8.1094E-04,4.9173E-04,1.9344E-04,2.7286E-05/
      DATA (FRACREFA(IG, 4),IG=1,16) /
     &2.9641E-01,2.1738E-01,1.4228E-01,1.0830E-01,8.2837E-02,6.1359E-02,
     &4.4683E-02,2.5027E-02,1.6057E-02,1.7558E-03,1.4193E-03,1.0970E-03,
     &7.8281E-04,4.3260E-04,1.4837E-04,2.2958E-05/
      DATA (FRACREFA(IG, 5),IG=1,16) /
     &2.9553E-01,2.2139E-01,1.4816E-01,1.0601E-01,8.0048E-02,6.0082E-02,
     &4.3952E-02,2.3788E-02,1.5734E-02,1.6586E-03,1.3434E-03,1.0281E-03,
     &7.0256E-04,4.2577E-04,1.2803E-04,1.3315E-05/
      DATA (FRACREFA(IG, 6),IG=1,16) /
     &2.9313E-01,2.2476E-01,1.5470E-01,1.0322E-01,7.8904E-02,5.8175E-02,
     &4.3097E-02,2.3618E-02,1.5385E-02,1.5942E-03,1.2702E-03,9.5566E-04,
     &6.5421E-04,4.0165E-04,1.2805E-04,1.3355E-05/
      DATA (FRACREFA(IG, 7),IG=1,16) /
     &2.9069E-01,2.2823E-01,1.5995E-01,1.0170E-01,7.7287E-02,5.6780E-02,
     &4.1752E-02,2.3899E-02,1.4937E-02,1.4916E-03,1.1909E-03,9.1307E-04,
     &6.3518E-04,3.9866E-04,1.2805E-04,1.3298E-05/
      DATA (FRACREFA(IG, 8),IG=1,16) /
     &2.8446E-01,2.2651E-01,1.7133E-01,1.0299E-01,7.4231E-02,5.6031E-02,
     &4.1368E-02,2.4318E-02,1.4135E-02,1.4216E-03,1.1465E-03,8.9800E-04,
     &6.3553E-04,3.9536E-04,1.2749E-04,1.3298E-05/
      DATA (FRACREFA(IG, 9),IG=1,16) /
     &2.0568E-01,2.5049E-01,2.0568E-01,1.1781E-01,7.5579E-02,5.8136E-02,
     &4.2397E-02,2.6544E-02,1.3067E-02,1.4061E-03,1.1455E-03,8.9408E-04,
     &6.3652E-04,3.9450E-04,1.2841E-04,1.3315E-05/

C Planck fraction mapping level : P=95.58350 mb, T = 215.70 K
      DATA FRACREFB /
     &1.8111E-01,2.2612E-01,1.6226E-01,1.1872E-01,9.9048E-02,8.0390E-02,
     &6.1648E-02,4.1704E-02,2.2976E-02,1.9263E-03,1.4694E-03,1.1498E-03,
     &7.9906E-04,4.8310E-04,1.6188E-04,2.2651E-05/

      EQUIVALENCE (KA,ABSA), (KB,ABSB)
      REAL KA,KB


C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.

C     P = 387. mb (Level 6)
      REFRAT_PLANCK_A = CHI_MLS(1,6)/CHI_MLS(6,6)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature,and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision: 3.3 $'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCH4(LAY)*COLCH4(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCH4_1(LAY)*COLCH4(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCH4(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(16) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(16) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG(16)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG(16)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)
            DO 2020 IG = 1, NG(16)
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF

 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(16) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(16) + 1
         DO 3000 IG = 1, NG(16)
            TAUG(LAY,IG) = COLCH4(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END
