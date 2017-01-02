      SUBROUTINE XS(species,nw,wavl,wav,T,DEN,j,sq,columndepth,zy,IO2)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      character*8 species

* input
!------- EWS - these are the variables that need to be declared here-------!
      dimension columndepth(KJ,NZ)  !testing this as an optional parameter 
      INTEGER nw,j,IO2
      REAL*8 wavl(nw+1), wav(nw)
      REAL*8 T(NZ), DEN(NZ)
      REAL*8 sq(kj,nz,kw)
      REAL*8 zy
      Ja=0
      Jb=1
      Jc=2
      Jd=3
!---------------------------------------------------------------------------!

!---------EWS NOTE: I have removed unnecessary variables to avoid compilation warnings---!
!---------EWS NOTE: Find original call formatting below update---------------------------!
!---------EWS NOTE: nw=nw, wl=wavl, wc=wav, tlev=T, airlev=DEN,jn=j, sq=sq---------------!
!----------------------------Comments Correspond to line below---------------------------!
      !T dependent - also unexplained factor of 1000. (but not using T-dependent data right now)
      IF(species.eq.'HNO3    ') CALL XS_HNO3(nw,wavl,T,j,sq) 
      IF(species.eq.'HO2     ') CALL XS_HO2(nw,wavl,wav,j,sq,Ja)
      !T-dependent (not in use)   
      IF(species.eq.'NO2     ') CALL XS_NO2(nw,wavl,T,j,sq)
      !T-dependent (not in use) 
      IF(species.eq.'H2O2    ') CALL XS_H2O2(nw,wavl,wav,T,j,sq)
      !HAVEN"T INTEGRATED GRACES CODE
      IF(species.eq.'H2O     ') CALL XS_H2O(nw,wavl,wav,j,sq)
      !returns CO+O to j and CO+O1D j+3 !ACK hardcoded needs fixing.. 
      IF(species.eq.'CO2     ') CALL XS_CO2(nw,wavl,j,sq)
      !returns H2, then HCO  
      IF(species.eq.'H2CO    ') CALL XS_H2CO(nw,wavl,j,sq)
      !returns O2 + O(1D), then O2 + O(3P)
      IF(species.eq.'O3      ') CALL XS_O3(nw,wavl,T,j,sq)
      IF(species.eq.'CH4     ') CALL XS_CH4(nw,wavl,j,sq)
      !returns  2 (3)CH2 +H2 then CH4 + (1)CH2  
      IF(species.eq.'C2H6    ') CALL XS_C2H6(nw,wavl,j,sq)
      !2 channels: NO + O2, NO2 + O
      !note NO3 is temperature dependent so would need to change if it does...
      IF(species.eq.'NO3     ') CALL XS_NO3(nw,wavl,T,j,sq)
      !returns N2 + O1D
      IF(species.eq.'N2O     ') CALL XS_N2O(nw,wavl,wav,T,j,sq)
      !returns OH + CL
      IF(species.eq.'HOCL    ') CALL XS_HOCL(nw,wavl,j,sq)
      !returns CL + CL (temp dependent)
      IF(species.eq.'CL2     ') CALL XS_CL2(nw,wavl,T,j,sq)
      !returns CLO + O
      IF(species.eq.'CLOO    ') CALL XS_CLOO(nw,wavl,j,sq)
      !returns CLO + O
      IF(species.eq.'OCLO    ') CALL XS_OCLO(nw,wavl,j,sq)
      !returns CL + NO2
      IF(species.eq.'CLONO   ') CALL XS_CLONO(nw,wavl,j,sq)
      !returns CL + NO
      IF(species.eq.'CLNO    ') CALL XS_CLNO(nw,wavl,j,sq)
      !returns CL + NO2
      IF(species.eq.'CLNO2   ') CALL XS_CLNO2(nw,wavl,j,sq)
      !returns HCO + CL (JPL calls this COHCL)  !products assumed from Yuk's model
      IF(species.eq.'CHCLO   ') CALL XS_CHCLO(nw,wavl,j,sq)
      !returns CL + CH3  (temperature dependent!)
      !there is a far UV branch to H + CH2CL , but no JPL recomended cross section so ignored
      IF(species.eq.'CH3CL   ') CALL XS_CH3CL(nw,wavl,wav,T,j,sq)
      !returns CL + CCL3  (temperature dependent!)
      IF(species.eq.'CCL4    ') CALL XS_CCL4(nw,wavl,wav,T,j,sq)
      !returns CL + CL + CO 
      IF(species.eq.'COCL2   ') CALL XS_COCL2(nw,wavl,j,sq)
      !returns CH3O2 + NO2 
      IF(species.eq.'CH3O2NO2') CALL XS_CH3O2NO2(nw,wavl,j,sq)
      !returns CH3O + CL 
      IF(species.eq.'CH3OCL  ') CALL XS_CH3OCL(nw,wavl,j,sq)
      IF(species.eq.'CH3OOH  ') CALL XS_CH3OOH(nw,wavl,j,sq)
      !returns CL + O1D and CL + O
      IF(species.eq.'CLO     ') CALL XS_CLO(nw,wavl,j,sq)
      !returns CL + NO3, then CLO + NO2 (temperature dependent!)
      IF(species.eq.'CLONO2  ') CALL XS_CLONO2(nw,wavl,wav,T,j,sq)
      !returns, HO2 + NO2, then OH + NO3
      IF(species.eq.'HO2NO2  ') CALL XS_HO2NO2(nw,wavl,j,sq)
      !returns NO3 + NO2, then NO3 + NO + O  (temperature dependent!)
      IF(species.eq.'N2O5    ') CALL XS_N2O5(nw,wavl,wav,T,j,sq)
      !returns CL + CLOO, then CLO + CLO 
      IF(species.eq.'CL2O2   ') CALL XS_CL2O2(nw,wavl,j,sq)
      !returns CL + CLO
      IF(species.eq.'CL2O    ') CALL XS_CL2O(nw,wavl,j,sq)
      !returns CLO + O2  
      IF(species.eq.'CLO3    ') CALL XS_CLO3(nw,wavl,j,sq)
      !returns CLOO + OCLO  
      IF(species.eq.'CL2O4   ') CALL XS_CL2O4(nw,wavl,j,sq)
      !returns H + CL. 
      IF(species.eq.'HCL     ') CALL XS_HCL(nw,wavl,j,sq)
      !returns O + O(1D) then O + O
      IF(species.eq.'O2      ') CALL XS_O2(nw,wavl,T,DEN,j,sq,
     $                                     columndepth,zy,IO2)
      !returns  signo(I,L) for now 
      IF(species.eq.'NO      ') CALL XS_NO(T,DEN,j,columndepth)  
      IF(species.eq.'SO      ') CALL XS_SO(nw,wavl,j,sq)
      IF(species.eq.'OCS     ') CALL XS_OCS(nw,wavl,j,sq)
      !returns SO + O, SO21, SO23
      IF(species.eq.'SO2     ') CALL XS_SO2(nw,wavl,j,sq)
      IF(species.eq.'SO3     ') CALL XS_SO3(nw,wavl,j,sq)
      IF(species.eq.'H2S     ') CALL XS_H2S(nw,wavl,j,sq)
      !returns SS8L then SS8R
      IF(species.eq.'S8      ') CALL XS_S8(nw,wavl,j,sq)
      !note using HO2 cross section for HSO
      IF(species.eq.'HSO     ') CALL XS_HO2(nw,wavl,wav,j,sq,Jb)
      !note this contains an unphysical scaling factor
      IF(species.eq.'S2      ') CALL XS_S2(nw,wavl,j,sq,Ja)
      !note using S2 cross section for S4
      IF(species.eq.'S4      ') CALL XS_S2(nw,wavl,j,sq,Jb)
      !note using S2 cross section for S3
      IF(species.eq.'S3      ') CALL XS_S2(nw,wavl,j,sq,Jc)
!sorg species
      IF(species.eq.'CS2    ') CALL XS_CS2(nw,wavl,j,sq)
      IF(species.eq.'CH3SH  ') CALL XS_CH3SH(nw,wavl,j,sq)

!below here there are calls to XS_simple.
      IF(species.eq.'C2H6S  ')CALL XS_simple(species,nw,wavl,j,sq)   
      IF(species.eq.'C2H6S2 ')CALL XS_simple(species,nw,wavl,j,sq)

!corg species
      IF(species.eq.'C2H2    ') CALL XS_C2H2(nw,wavl,j,sq) 
      IF(species.eq.'C2H4    ') CALL XS_C2H4(nw,wavl,j,sq,Ja)
      IF(species.eq.'CH      ') CALL XS_C2H4(nw,wavl,j,sq,Jd) ! shawn used C2H4 rate here  
      IF(species.eq.'C3H8    ') CALL XS_C3H8(nw,wavl,j,sq)
      IF(species.eq.'CH2CO   ') CALL XS_simple(species,nw,wavl,j,sq)
      IF(species.eq.'CH3CHO  ') CALL XS_CH3CHO(nw,wavl,j,sq)  !shawn used C2H4 rate here with 0.57/0.02/0.05 quantum yields
      IF(species.eq.'C3H6    ') CALL XS_C3H6(nw,wavl,j,sq)
      IF(species.eq.'C2H5CHO ') CALL XS_simple(species,nw,wavl,j,sq)
      IF(species.eq.'C3H3    ') CALL XS_simple(species,nw,wavl,j,sq)    !shawn used C2H4 rate here (check this xsection)
      IF(species.eq.'CH3C2H  ') CALL XS_CH3C2H(nw,wavl,j,sq)
      IF(species.eq.'CH2CCH2 ') CALL XS_CH2CCH2(nw,wavl,j,sq)

c these are just here to remind us that these XS files exist
c       CALL XS_CL2O3(nw,wavl,wav,T,DEN,JCL2O3,sq)  !returns CLO + CLOO  (branch to CLO3??)  

c       CALL XS_CHOCHO(nw,wavl,wav,T,DEN,100,sq,101,102)  !ACK - hardcoded "J-numbers" 
       !returns HCO + HCO, then H2 + CO + CO, then H2CO + CO  (never finished)

      RETURN
      END


       SUBROUTINE XS_HNO3(nw,wl,tlev,jn,sq)

!---------------------------------------------------------------------------------------!
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
*=        HNO3 + hv -> OH + NO2                                              =*
*=  Cross section: Burkholder et al., 1993                                   =*
*=  Kevin's photo.dat data (old JPL recomendation)                           =*
*=                                                                           =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  Jn      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel,plab
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
c      REAL*8 wl(kw), wc(kw)
      REAL*8 wl(nw+1) !EWS wc not needed 
      !Eddie - Currently testing commenting out the below. Compilation warning says it isn't used
      !  but it is. curious. 
      REAL*8 tlev(nz)  ! EWS - used here
      !I am no longer so sure about the above comment - NZ is set by the time we invoke XS's
      !it's not like we are going to change on the fly in the middle of a time-dependent run....
c      REAL*8 airlev(kz) ! currently this isn't used until declared again below

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1, n2
      REAL*8 x1(kdata), x2(kdata)
      REAL*8 y1(kdata), y2(kdata)

* local
      REAL*8 yg1(nw), yg2(nw)
      INTEGER i, iw
      INTEGER ierr
      INTEGER option
      ierr = 0


**************** HNO3 photodissociation
c option 1 -  Burkholder et al. 1993
c option 2 -  Kevin's photo.dat

      option=2

      IF (option .eq. 1) then
* HNO3 cross section parameters from Burkholder et al. 1993

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/HNO3/HNO3_burk.abs',
     &  STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n1 =  83
      n2 = n1
      DO i = 1, n1
         READ(kin,*) y1(i), y2(i)
         x1(i) = 1840. + i*20.         ! wavelengths: from 186 to 350 nm by 2 nm
         x2(i) = x1(i)


      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  !adds 185.9814 to x1, 0 to y1
      CALL addpnt(x1,y1,kdata,n1,               zero,zero) !adds 0,0 to x1,y1
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero) !adds 350.035, 0
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero) !adds 1e36,0

      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_HNO3***'
         STOP
      ENDIF

c note different behavior in addpnt below than above 
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))  !adds 185.914,1.7
      CALL addpnt(x2,y2,kdata,n2,               zero,y2(1))!adds 0, 1.7
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))!adds 3500.35,9.3
      CALL addpnt(x2,y2,kdata,n2,            biggest,y2(n2))!adds 1e36, 9.3
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'something wrong in XS_HNO3'
         STOP
      ENDIF


* quantum yield = 1
* correct for temperature dependence

      DO iw = 1, nw
         DO i = 1, nz
            sq(jn,i,iw) = yg1(iw) * 1.E-20
     $           * exp( yg2(iw)/1.e3*(tlev(i)-298.) ) !why is this divided by 1000 here?
         ENDDO
      ENDDO
      endif  !end option 1


      if (option .eq. 2) then  !use Kevin's data, which is from "some old JPL recommendation"
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/HNO3/HNO3_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF  

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_HNO3***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1

      DO iw = 1, nw
         DO i = 1, nz
            sq(jn,i,iw) = yg1(iw)
         ENDDO
      ENDDO
      endif   !end option 2 

      plab='PHNO3'
      photolabel(jn)=plab
      jn=jn+1
 

      RETURN
      END

       SUBROUTINE XS_HO2(nw,wl,wc,jn,sq,jdum)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HO2 photolysis:    =*
*=          HO2 + hv -> OH + O                                               =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed shape based on work by Lee, 1982; normalized      =*
*=                 to unity at 248 nm                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  Jn      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE
c      INTEGER kin,kj,kw 
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn,jdum
      REAL*8 wl(nw+1), wc(nw) !EWS - wc used here
      !orig was kw, but don't see why this is needed once nw is defined by gridw.f 
      !could go back to kz if I ever want to have a variable altitude grid?
c      REAL*8 tlev(nz)  ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1 
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw), yg(nw)   !orig kw
      REAL*8 qy
      INTEGER i, iw, n, iz
      INTEGER ierr
      INTEGER option
      ierr = 0

**************************************************************
************* HO2 photodissociation
c option 1) cross sections from JPL97 recommendation (identical to 94 recommendation)+ quantum yield
c option 2) from Kevin no qy

      option=2

      if (option.eq.1) then

      OPEN(unit=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/HO2/HO2_jpl94.abs',
     &  STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        x1(i)=x1(i)*10.
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n,          zero,zero)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n,        biggest,zero)
      CALL inter2(nw+1,wl,yg,n,x1,y1,ierr)
       
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr,' Something wrong in XS_HSO2'   !???
        STOP
      ENDIF

**** quantum yield:  absolute quantum yield has not been reported yet, but
****                 Lee measured a quantum yield for O(1D) production at 248
****                 nm that was 15 time larger than at 193 nm
**** here:  a quantum yield of unity is assumed at 248 nm and beyond, for
****        shorter wavelengths a linear decrease with lambda is assumed

      DO iw = 1, nw
         IF (wc(iw) .GE. 248.) THEN
            qy = 1.
         ELSE
            qy = 1./15. + (wc(iw)-193.)*(14./15.)/(248.-193.)
            qy = MAX(qy,0.)
         ENDIF
         DO iz = 1, nz
           sq(jn,iz,iw) = qy * yg(iw)
         ENDDO
      ENDDO

      endif  !end option 1

      if (option .eq. 2) then  !use Kevin's data, which is from "some old JPL recommendation"
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/HO2/HO2_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  35
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_HO2***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1

      DO iw = 1, nw
         DO i = 1, nz
            sq(jn,i,iw) = yg1(iw)
         ENDDO
      ENDDO
      endif   !end option 2 

      if(jdum.eq.0) photolabel(jn)='PHO2'

      if(jdum.eq.1) then
        photolabel(jn)='PHSO'
      endif   

      jn=jn+1

      RETURN
      END
 
       SUBROUTINE XS_NO2(nw,wl,tlev,jn,sq)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for NO2            =*
*=  photolysis:                                                              =*
*=         NO2 + hv -> NO + O(3P)                                            =*
*=  Cross section from JPL94 (can also have Davidson et al.)                 =*  !UPDATE ME
*=  Quantum yield from Gardiner, Sperry, and Calvert                         =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=          see above examples                                               =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel,plab
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (NO2)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata), x2(kdata), x3(kdata)
      REAL*8 y1(kdata), y2(kdata)

* local
      REAL*8 yg1(nw), yg2(nw)  !orig kw
      REAL*8 dum
      INTEGER i, iw, n, idum, ierr,option
      REAL*8 xsno2(nz,nw),wavn
      ierr = 0

* local


**************** NO2 photodissociation

* cross section

*------------NEED TO CHANGE kdata = 1000 FOR DAVIDSON ET AL. DATA---------
*     measurements of Davidson et al. (198x) at 273K
*     from 263.8 to 648.8 nm in approximately 0.5 nm intervals
C     OPEN(UNIT=kin,
C    &  file=PHRT(1:INDEX(PHRT,' ')-1)//'DATAE1/NO2/NO2_ncar_00.abs',
C    &  STATUS='old')
C     n = 750
C     DO i = 1, n
C        READ(kin,*) x1(i), y1(i), dum, dum, idum
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),zero)
C     CALL addpnt(x1,y1,kdata,n,               zero,zero)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),zero)
C     CALL addpnt(x1,y1,kdata,n,           biggest,zero)
C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)
C     IF (ierr .NE. 0) THEN
C        WRITE(*,*) ierr, jlabel(j)
C        STOP
C     ENDIF

* options:
*       1    NO2_jpl94.abs  (same as JPL97)
*       2    HarNO2CS.rxt  from Harder et al.
*       3    Kevin's data from photo.dat

      option = 3

* cross section data from JPL 94 recommendation
* JPL 97 recommendation is identical

      if (option .eq. 1) then

         OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/NO2/NO2_jpl94.abs',
     &  STATUS='old')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO

* read in wavelength bins, cross section at T0 and temperature correction
* coefficient a;  see input file for details.
* data need to be scaled to total area per bin so that they can be used with
* inter3

         DO i = 1, n
            READ(kin,*) x1(i), x3(i), y1(i), dum, y2(i)
            x1(i)=x1(i)*10. !convert to angstroms
            x3(i)=x3(i)*10.
            y1(i) = (x3(i)-x1(i)) * y1(i)*1.E-20
            y2(i) = (x3(i)-x1(i)) * y2(i)*1.E-22
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)

         x1(n+1) = x3(n)
         x2(n+1) = x3(n)
         n = n+1
         n1 = n

         CALL inter3(nw+1,wl,yg1,n,x1,y1,0)
         CALL inter3(nw+1,wl,yg2,n1,x2,y2,0)

* yg1, yg2 are per nm, so rescale by bin widths

         DO iw = 1, nw
            yg1(iw) = yg1(iw)/(wl(iw+1)-wl(iw))
            yg2(iw) = yg2(iw)/(wl(iw+1)-wl(iw))
         ENDDO

         DO iw = 1, nw
            DO i = 1, nz
               xsno2(i,iw) = yg1(iw) + yg2(iw)*(tlev(i)-273.15)
            ENDDO
         ENDDO


      elseif (option .eq. 2) then

         OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/NO2/NO2_Har.abs',
     &  status='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 135
         DO i = 1, n
            READ(kin,*) idum, y1(i)
            x1(i) = FLOAT(idum)*10.  !convert to angstroms
         enddo

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n,               zero,y1(1))
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   zero)
         CALL addpnt(x1,y1,kdata,n,           biggest,   zero)
         CALL inter2(nw+1,wl,yg1,n,x1,y1,ierr)

         DO iw = 1, nw
            DO i = 1, nz
               xsno2(i,iw) = yg1(iw)
            ENDDO
         ENDDO

      elseif (option .eq. 3) then

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/NO2/NO2_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 
     

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_NO2***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

      DO iw = 1, nw
         DO i = 1, nz
            xsno2(i,iw) = yg1(iw)
         ENDDO
      ENDDO

      ENDIF
         

      if (option.eq.1 .or. option.eq.2) then
* quantum yield
* from Gardiner, Sperry, and Calvert

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/NO2/NO2_calvert.yld',
     &  STATUS='old')
      DO i = 1, 8
         READ(kin,*)
      ENDDO
      n = 66
      DO i = 1, n
         READ(kin,*) x1(i),y1(i)
         x1(i)=x1(i)*10.
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               zero,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   zero)
      CALL addpnt(x1,y1,kdata,n,           biggest,   zero)
      CALL inter2(nw+1,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' Something wrong in XS_NO2'
         STOP
      ENDIF

      else if (option.eq.3) then   !Kevin's quantum yield
c - qy needs to be called yg1

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/NO2/NO2qy_zahnle.yld',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

C   CALCULATE NO2 QUANTUM YIELDS
c   from Kevin's original code...

c      DO 47 L=47,60
c      WAVN = WAV(L)/10.     !wav(L) is centered wavelength...
c  47  QYNO2(L) = 1. - 8.E-4*(WAVN - 275.)
C

      do iw=1,n1-1
         if (x1(iw).ge.2941.0 .AND. x1(iw).le.3650.0) then
          wavn =0.5*(x1(iw)+x1(iw+1))/10.
          y1(iw) = 1. - 8.E-4*(wavn - 275.)
         endif
      enddo

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' Something wrong in XS_NO2'
         STOP
      ENDIF
 
      endif    !end quantum yield loop

* combine and return xs*qy

      DO iw = 1, nw
         DO i = 1, nz
            sq(jn,i,iw) = xsno2(i,iw)*yg1(iw)
         ENDDO
      ENDDO

      plab='PNO2'
      photolabel(jn)=plab

c      print *, jn,photolabel

      jn=jn+1

      RETURN
      END

       ! EWS - airlev not used below
c      SUBROUTINE XS_H2O2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_H2O2(nw,wl,wc,tlev,jn,sq)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for H2O2 photolysis =*
*=         H2O2 + hv -> 2 OH                                                 =*
*=  Cross section:  From JPL97, tabulated values @ 298K for <260nm, T-depend.=*  !UPDATE ME
*=                  parameterization for 260-350nm                           =*
*=  Quantum yield:  Assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'

      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel,plab
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc' 
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1), wc(nw) ! EWS - wc needed here (H2O2)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg(nw), yg1(nw)
      REAL*8 qy
      REAL*8 a0, a1, a2, a3, a4, a5, a6, a7
      REAL*8 b0, b1, b2, b3, b4
      REAL*8 xs
      REAL*8 t
      INTEGER i, iw, n, idum
      INTEGER ierr,option
      REAL*8 lambda
      REAL*8 sumA, sumB, chi
      ierr = 0      

**************** H2O2 photodissociation

c options
c     1)temperature-dependent JPL 1994/97
c     2)Kevin's photo.dat data

      option=2

      if (option.eq.1) then
* cross section from JPL94 (identical to JPL97)
* tabulated data up to 260 nm

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/H2O2/H2O2_jpl94.abs',
     &  STATUS='old')
      READ(kin,*) idum,n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)


      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n,               zero,zero)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n,           biggest,zero)
      CALL inter2(nw+1,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' Something wrong in XS_H2O2'
         STOP
      ENDIF

      A0 = 6.4761E+04
      A1 = -9.2170972E+02
      A2 = 4.535649
      A3 = -4.4589016E-03
      A4 = -4.035101E-05
      A5 = 1.6878206E-07
      A6 = -2.652014E-10
      A7 = 1.5534675E-13

      B0 = 6.8123E+03
      B1 = -5.1351E+01
      B2 = 1.1522E-01
      B3 = -3.0493E-05
      B4 = -1.0924E-07

* quantum yield = 1

      qy = 1.

      DO iw = 1, nw

* Parameterization (JPL94)
* Range 260-350 nm; 200-400 K

         IF ((wl(iw) .GE. 2600.) .AND. (wl(iw) .LT. 3500.)) THEN

           lambda = wc(iw)/10.   !convert to nm
           sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda +
     >                  A4)*lambda +A3)*lambda + A2)*lambda +
     >                  A1)*lambda + A0
           sumB = (((B4*lambda + B3)*lambda + B2)*lambda +
     >               B1)*lambda + B0

           DO i = 1, nz
              t = MIN(MAX(tlev(i),200.),400.)
              chi = 1./(1.+EXP(-1265./t))
              xs = (chi * sumA + (1.-chi)*sumB)*1E-21
              sq(j,i,iw) = xs*qy
           ENDDO
         ELSE
           DO i = 1, nz
              sq(j,i,iw) = yg(iw)*qy
           ENDDO
         ENDIF

      ENDDO
      endif  !end option 1
      
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 
C         print*,"HJtest",HJtest

      if (option.eq.2) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/H2O2/H2O2_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)  
      ENDIF

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_H2O2***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)
         ENDDO
      ENDDO
       endif  !end option 2

      plab='PH2O2'
      photolabel(jn)=plab
      jn=jn+1

      RETURN
      END


       SUBROUTINE XS_OCS(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for OCS photolysis  =*
*=         OCS + hv -> CO + S                                                =*
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  Kevin had this as 0.72 at all hieghts/wls                =*
*=                  but what is the other channel?
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (OCS)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** OCS photodissociation

c options
c     1)Kevin's photo.dat data



      option=1

      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/OCS/OCS_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_OCS***'
         STOP
      ENDIF


c Quantum yield (from Kevin's code - no reference)
      qy=0.72

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

         photolabel(jn)='POCS'

      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_SO3(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for SO3 photolysis  =*
*=         SO3 + HV -> SO2 + O                                               =*
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  Assuming 1 with no research                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE
c      INTEGER kin,kj,kw 
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (SO3)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** SO3 photodissociation

c options
c     1)Kevin's photo.dat data

      option=1

      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/SO3/SO3_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_SO3***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

c Quantum yield (no reference nor research)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

         photolabel(jn)='PSO3'

      jn=jn+1


      RETURN
      END
       
       SUBROUTINE XS_S2(nw,wl,jn,sq,jdum)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for S2 photolysis   =*
*=         S2 + HV -> S + S                                                  =*
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=   NOTE - arbitrarily scaled by 6.53 as in Kevin's code --- WHY?
*=  Quantum yield:  Assuming 1 with no research                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)

      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (S2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy, scale
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** S2 photodissociation

c options
c     1)Kevin's photo.dat data

      option=1

      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/S2/S2_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_S2***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

c Quantum yield (no reference nor research)
c scale is Kevin's factor to (and I quote...)
C   SCALE UP S2 CROSS SECTIONS TO OBTAIN CORRECT S2 LIFETIME UP HIGH

      qy=1.0
      scale=6.53

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy*scale
         ENDDO
      ENDDO



      endif  !end option 1
      
c      print *, jn,photolabel

      if(jdum.eq.0) photolabel(jn)='PS2'
      if(jdum.eq.1) photolabel(jn)='PS4'
      if(jdum.eq.2) photolabel(jn)='PS3'

      
      jn=jn+1

      
      RETURN
      END
      
       SUBROUTINE XS_S8(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for S8  photolysis  =*
*=         S8 + HV -> S4 + S4                                                =*
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  Assuming 1 with no research                              =*
*  NOTE - this is pretty kludgy... and should be explained in detail here    =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE 
c      INTEGER kin,kj,kw
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
c orig      REAL*8 wl(nw), wc(nw)   !orig was kw, but don't see why this is needed once nw is defined by gridw.f
      REAL*8 wl(nw+1) ! EWS - wc not needed (S8)


c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=120)
      INTEGER n1
      REAL*8 x1(kdata), x2(kdata)
      REAL*8 y1(kdata), y2(kdata)

* local
      REAL*8 yg1(nw), yg2(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0           

**************** S8 photodissociation
c options
c     1)Kevin's photo.dat data

      option=1

      if (option.eq.1) then  !Kevin's data

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/S8/S8L_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  108
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)


      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_S8***'
         STOP
      ENDIF


      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/S8/S8R_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x2(i), y2(i)
         y2(i)=y2(i)*1.E-17      !data given in units of sigma*1E17
      ENDDO
      CLOSE (kin)
      


      CALL addpnt(x2,y2,kdata,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n1,               zero,zero)
      CALL addpnt(x2,y2,kdata,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg2,n1,x2,y2,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg2,n1,x2,y2,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_S8***'
         STOP
      ENDIF

c      endif

c Quantum yield (no reference nor research)
      qy=1.0

ccc  Kevin's original code
c        do L=1,68    
c          SS8L(L) = MAX(SS8R(L),SS8L(L))  
c        enddo

      do iw=1,68                          !ACK
         yg1(iw)=MAX(yg1(iw),yg2(iw))
      enddo   

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy       !returns SS8L
              sq(jn+1,i,iw) = yg2(iw)*qy     !then SS8R
         ENDDO
      ENDDO
      endif      

      photolabel(jn)='PS8L'
      jn=jn+1

      photolabel(jn)='PS8R'
      jn=jn+1

      photolabel(jn)='PS8'
      jn=jn+1

      !so we are reserving sq(j+3)=0.0 for S8.
      !the actual S8 photorate is calcuated from combination of the S8L and S8R cross sections

      RETURN
      END

       SUBROUTINE XS_H2S(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for H2S photolysis  =*
*=         H2S + HV -> HS + H                                                =*
*=  Cross section:  From Kevin's photo.dat                                   =* 
*=  Cross section:  data from MPI db.  Temperature dependence ignored        =*  
*=  Quantum yield:  Assuming 1 with no research                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE
c      INTEGER kin,kj,kw 
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (H2S)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=21000)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** H2S photodissociation
c options
c     1)Kevin's photo.dat data
c     2)MPI data, high res from 160-260 (not there is some temperature dependence here which is not included), along with far UV data not in photo.dat

      option=2


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/H2S/H2S_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_H2S***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

c Quantum yield (no reference nor research)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      if (option.eq.2) then  !mpi data and extending into the far uv
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/H2S/H2S_mpi.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  20004
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   !inter2 since points

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_H2S***'
         STOP
      ENDIF

c      do i=1,nw
c         print *, wl(i),yg1(i)
c      enddo
c      stop


c Quantum yield (no reference nor research)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 2

      photolabel(jn)='PH2S'

      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_H2O(nw,wl,wc,jn,sq)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield) for H2O photolysis:     !
!           h2o + hv -> h + oh                                                !
!   cross section: taken from three sources                                   !
!     1) jpl94-26 (jpl97-4), 175.5 - 189.3                                                !
!     2) cantrell et al., grl, 24, 17, 2195-2198, 1997,  183.0 - 193.0 nm     !
!     3) yoshino et al.,  chemical physics, 211 (1996) 387-391, 120.38-188.03 !
!                                                                             !
!   quantum yield: is unity between 175.5 and 189.3   (jpl97-4)               !
!-----------------------------------------------------------------------------!

*** NOTE - I haven't integrated the graces code yet ****


c      IMPLICIT NONE
c      INTEGER kin,kj,kw 
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  
       !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc' 
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1), wc(nw) ! EWS - wc needed here (H2O)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* local
      INTEGER ierr,option
      integer kdata
      parameter(kdata=7000)
      real*8 x1(kdata)
      real*8 y1(kdata)
      real*8 x2(kdata)
      real*8 y2(kdata)
      real*8 x3(kdata)
      real*8 y3(kdata)
      real*8 yg(nw), yg1(nw), yg2(nw), yg3(nw)
      real*8 qy
      INTEGER i, n2, n3, idum,iz
      ierr = 0

***************** H2O *************
c options
c     1)from Graces - combination of three (JPL 97 recs)
c     2) Kevin's photo.dat file

      option=2

***NOTE option 1 will not work - need to implement ***


      if (option.eq.1) then

!----------------------------------------------
!     ... jlabel(j) = 'h2o -> prod'
!----------------------------------------------
c      j = j+1
c      jlabel(j) = 'h2o + hv -> h + oh'
c      j = j+1
c      jlabel(j) = 'h2o + hv -> h2 + o(1d)'
c        j = j+1
c      jlabel(j) = 'h2o + hv -> 2h + o(3p)'

!----------------------------------------------------
!     ... cross sections from jpl97 recommendation
!----------------------------------------------------
      open(kin,
     &   file='KDATA/CROSS_SECTIONS/H2O_jpl94.abs',
     & status='old')
      read(kin,*) idum, n1
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n1
        read(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1e-20
      enddo
      close(kin)

      call addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      call addpnt(x1,y1,kdata,n1,          zero,zero)
      call addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      call addpnt(x1,y1,kdata,n1,        biggest,zero)

      call inter2(nw,wl,yg1,n1,x1,y1,ierr)

      if (ierr .ne. 0) then
c        write(*,*) ierr, jlabel(j)
        stop
      endif
!----------------------------------------------------
!     ...cross sections from cantrell et al., 1997
!----------------------------------------------------
      open(kin,
     &  file='KDATA/CROSS_SECTIONS/H2O_cantrell_1996.abs',
     & status='old')
      read(kin,*) idum, n2
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n2
        read(kin,*) x2(i), y2(i)
        y2(i) = y2(i) * 1e-20
      enddo
      close(kin)

      call addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)
      call addpnt(x2,y2,kdata,n2,          zero,zero)
      call addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      call addpnt(x2,y2,kdata,n2,        biggest,zero)

      call inter2(nw,wl,yg2,n2,x2,y2,ierr)

      if (ierr .ne. 0) then
c        write(*,*) ierr, jlabel(j)
        stop
      endif
!----------------------------------------------------
!     ... cross sections from yoshino et al., 1996
!----------------------------------------------------
      open(kin,
     &  file='KDATA/CROSS_SECTIONS/H2O_yoshino_1996.abs',
     & status='old')
      read(kin,*) idum, n3
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n3
        read(kin,*) x3(i), y3(i)
      enddo
      close(kin)

      call addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)
      call addpnt(x3,y3,kdata,n3,          zero,zero)
      call addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      call addpnt(x3,y3,kdata,n3,        biggest,zero)

      call inter2(nw,wl,yg3,n3,x3,y3,ierr)

      if (ierr .ne. 0) then
c        write(*,*) ierr, jlabel(j)
        stop
      endif
!--------------------------------------------------------------------------------
!     ...combine data sets (i.e., yoshino et al., 1996 and cantrell et al., 1997)
!--------------------------------------------------------------------------------
      do i = 1, nw-1
        if (wc(i) .lt. 183.0) then
            yg(i) = yg3(i)
          elseif (wc(i) .le. 194.0) then
            yg(i) = yg2(i)
          else
            yg(i) = 0.
          endif
        enddo
!------------------------------------------------------
!     ... quantum yield assumed to be unity (jpl97-4)
!------------------------------------------------------
!     ... 105 to 145 nm
!         (JPL 1997 which references Stief, L.J., W.A.
!         Payne, and R. B. Klemm, A flash
!         photolysis-resonance fluoresence study of the
!         formation of O(1D) in the photolysis of water
!         and the reaction of O(1D) with H2, Ar, and He,
!         J. Chem. Phys., 62, 4000, 1975.)

      do iw = 1, nw-1

        if (wc(iw) .le. 145.0) then

             do iz = 1, nz
             sq(j-2,iz,iw) = yg(iw) * 0.890
               sq(j-1,iz,iw) = yg(iw) * 0.110
               sq(j,  iz,iw) = yg(iw) * 0.0
           end do

          end if

!     ... > 145nm
!         JPL97
        if (wc(iw) .gt. 145.0) then

             do iz = 1, nz
               sq(j-2,iz,iw) = yg(iw) * 1.0
                 sq(j-1,iz,iw) = yg(iw) * 0.0
                 sq(j,  iz,iw) = yg(iw) * 0.0
           end do

          end if

        end do! end wavelength loop

!     ... Overwrite Lyamn Alpha
!         Slanger, T.G., and G. Black, Photodissociative
!         channels at 1216A for H2O, NH3 and CH4,
!         J. Chem. Phys., 77, 2432, 1982.)
!         **** Caution, Lyman Alpha is always wavelength postion #2 ****
      iw = 2
      do iz = 1, nz
        sq(j-2,iz,iw) = yg(iw) * 0.780
          sq(j-1,iz,iw) = yg(iw) * 0.100
          sq(j,  iz,iw) = yg(iw) * 0.120
      end do

      endif   !end option 1)
      
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 
C         print*,"HJtest",HJtest
      
      if (option.eq.2) then

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/H2O/H2O_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  45
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJTEST.EQ.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)
      ELSE
       	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)     
      ENDIF

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_H2O***'
         STOP
      ENDIF
c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

c Quantum yield (no reference nor research)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif   !end option 2


      photolabel(jn)='PH2O'
      jn=jn+1

      RETURN
      END 
c      end subroutine XS_H2O

       SUBROUTINE XS_SO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for SO photolysis   =*
*=         SO + HV -> S + O                                                  =*
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  Assuming 1 with no research                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE
c      INTEGER kin,kj,kw 
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (SO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** SO photodissociation
c options
c     1)Kevin's photo.dat data

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/SO/SO_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  78
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_SO***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

c Quantum yield (no reference nor research)

      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PSO'


      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_CO2(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CO2 photolysis  =*
*=         CO2 + HV -> CO + O                                                =*
*=         CO2 + HV -> CO + O1D                                              =*
*=    Kevin's code takes everything in the sw loop to O3P                    =*
*=    and everything in Far UV loop to O1D                                   =*

*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  Assuming 1 with no research                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CO2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata), x2(kdata)
      REAL*8 y1(kdata), y2(kdata)

* local
      REAL*8 yg1(nw),yg2(nw)
c      REAL*8 qy !EWS - not used
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** CO2 photodissociation
c options
c     1)Kevin's photo.dat data
c     2)MC temp-dependent merge to zahnle grid (should re-evaluate)
        !this data is for 195K - in the final run should be t-dependent


      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 
C         print*,"HJtest",HJtest

      option=1

      if (option.eq.1) then  !Kevin's data

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CO2/CO2_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  35
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)         
      ENDIF
      
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CO2***'
         STOP
      ENDIF

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CO2/CO2D_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  10
      DO i = 1, n1
         READ(kin,*) x2(i), y2(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n1,               zero,zero)
      CALL addpnt(x2,y2,kdata,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg2,n1,x2,y2,ierr)
      ELSE
      	CALL inter2(nw+1,wl,yg2,n1,x2,y2,ierr)   
      ENDIF

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CO2D***'
         STOP
      ENDIF

c      do i=1,118
c         print *, wl(i),yg1(i),yg2(i)
c      enddo
c      stop

      endif !end option 1

      if (option.eq.2) then  !My gridded data (should perhaps re-evalutate)

!ACK - need to integrate below...

c-mc
c-mc yosh2kz_disw vector from ~/dev/CO2/co2_xsection.pro - 195K cross section
c-mc data from parkinson/yoshino as binned to this wl grid.
c      DATA SCO2ASW/1.4E-19, 7.1E-19, 5.4E-19, 5.2E-19, 5.1E-19, 3.7E-19,
c     2  2.1E-19, 7.4E-20, 3.1E-20, 9.7E-21/ 




      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CO2/CO2_claire1.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  35
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CO2***'
         STOP
      ENDIF



      endif ! end option 2


      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)
              sq(jn+1,i,iw) = yg2(iw)  !CO2_O1D
         ENDDO
      ENDDO

c      print *, jn, photolabel

      photolabel(jn)='PCO2'
      jn=jn+1

      photolabel(jn)='PCO2_O1D'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, wc not used
c      SUBROUTINE XS_H2CO(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_H2CO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for H2CO photolysis =*
*=         H2CO + HV -> H2 + CO          (reac 38)                           =*
*=         H2CO + HV -> HCO + H          (reac 39)                           =*
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  included in Kevin's cross sections                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE 
c      INTEGER kin,kj,kw
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (H2CO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=120)
      INTEGER n1, n2,n3
      REAL*8 x1(kdata), x2(kdata), x3(kdata)
      REAL*8 y1(kdata), y2(kdata), y3(kdata)

* local
      REAL*8 yg1(nw), yg2(nw), yg3(nw)
c      REAL*8 qy !EWS - not used
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** H2CO photodissociation
c options
c     1)Kevin's photo.dat data

      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 
C         print*,"HJtest",HJtest

      option=1

      if (option.eq.1) then  !Kevin's data

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/H2CO/H2CO_zahnle.abs',
     &  STATUS='old')

c pointer

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 = 68
      n2 = 68
      n3 = 68
      DO i = 1, n1
         READ(kin,*) x1(i),y1(i),y2(i),y3(i)
         x2(i)=x1(i)
         x3(i)=x1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)
      ENDIF

!yg1 is H2CO cross section

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_H2CO***'
         STOP
      ENDIF


c Quantum yields:  from Kevin's photo.dat file

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n2,               zero,zero)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,            biggest,zero)
      IF (HJtest.eq.1) THEN
       CALL inter3(nw+1,wl,yg2,n2,x2,y2,ierr)
      ELSE
       CALL inter2(nw+1,wl,yg2,n2,x2,y2,ierr)
      ENDIF

!yg2 is quantum yield for HCO


      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_H2CO***'
         STOP
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)  
      CALL addpnt(x3,y3,kdata,n3,               zero,zero)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,            biggest,zero)
      IF (HJtest.eq.1) THEN
       CALL inter3(nw+1,wl,yg3,n3,x3,y3,ierr)
      ELSE
       CALL inter2(nw+1,wl,yg3,n3,x3,y3,ierr)      
      ENDIF

!yg3 is quantum yield for H2

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_H2CO***'
         STOP
      ENDIF



      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*yg3(iw)   !H2CO -> H2 
              sq(jn+1,i,iw) = yg1(iw)*yg2(iw)     !H2CO -> HCO 

         ENDDO
      ENDDO
      endif   ! end option 1

      photolabel(jn)='PH2CO_H2'
      jn=jn+1

      photolabel(jn)='PH2CO_HCO'
      jn=jn+1

      RETURN
      END


       SUBROUTINE XS_SO2(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for SO2 photolysis  =*
*=         SO2 + HV -> SO + O                                                =*
*=         SO2 + HV -> SO21                                                  =*
*=         SO2 + HV -> SO23                                                  =*
*=  Cross section:  From Kevin's photo.dat or highres if LGRID=1             =*  
*=  Quantum yield:  included in Kevin's cross sections                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE 
c      INTEGER kin,kj,kw
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (SO2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata,kdataHR
      PARAMETER(kdata=120,kdataHR=1200)
      INTEGER n1, n2,n3
      REAL*8 x1(kdata), x2(kdata), x3(kdata)
      REAL*8 y1(kdata), y2(kdata), y3(kdata)

      REAL*8 x32(kdataHR), x33(kdataHR), x34(kdataHR),xdum(kdataHR)
      REAL*8 y32(kdataHR), y33(kdataHR), y34(kdataHR)
      REAL*8 ydum1(kdataHR), ydum2(kdataHR), ydum3(kdataHR)


* local
      REAL*8 yg1(nw), yg2(nw), yg3(nw), ygnew(nw)
c      REAL*8 yg32(nw),yg33(nw), yg34(nw) ! EWS - not used

c      REAL*8 qy !EWS - not used
      real*8 mass  !from PHOTABLOK
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** SO2 photodissociation
C  112)  SO2 + HV -> SO + O
C  144)  SO2 + HV -> SO21
C  145)  SO2 + HV -> SO23


c options
c     1)Kevin's photo.dat data (for all three)
c - note that these cross sections contain quantum yields...
c - JPL -06 recommends looking at Hudson and Kieffer - 1975 natural stratosphere of 1974 CIAP monograph...
c JPL -02 has some better comments - need quantum yields to do this carefully.




      option=1 
      
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 
C         print*,"HJtest",HJtest

      if (option.eq.1) then  !Kevin's data

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/SO2/SO2_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 = 78
      n2 = 78
      n3 = 78
      DO i = 1, n1
         READ(kin,*) x1(i),y1(i),y2(i),y3(i)
         x2(i)=x1(i)
         x3(i)=x1(i)
      ENDDO
      CLOSE (kin)

C  Interpolate the cross sections to the wavelength grid:
c- note these are overwritten by higher resolution data below if LGRID=1

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 
!yg1 is SO + O cross section * quantum yield

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_SO2***'
         STOP
      ENDIF


      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n2,               zero,zero)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,            biggest,zero)
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,ierr)
!yg2 is SO21 cross section * quantum yield

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_SO2***'
         STOP
      ENDIF



      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)  
      CALL addpnt(x3,y3,kdata,n3,               zero,zero)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,            biggest,zero)
      CALL inter2(nw+1,wl,yg3,n3,x3,y3,ierr)
!yg3 is SO23 cross section * quantum yield

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_SO2***'
         STOP
      ENDIF


      if (LGRID.eq.1) then   !on fine wavelength grid, use Danielache data

      newspec=1  !use Sebastian's corrected spectrum

      if (newspec.eq.0) then
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/SO2/2007jd009695-ds01.txt',
     &  STATUS='old') 
      n4 = 3360
!the default used for the sulfur mif calculations in the 2013 paper - from Danielace 2012 JGR paper
      else

      OPEN(UNIT=kin,file=
     & 'PHOTOCHEM/DATA/XSECTIONS
     &/SO2/data-non-systematic-errors-corrections.txt',
     &  STATUS='old')
      n4 = 1184
!Seba said the published versions are incorrect and sent these along in 2013

c longward1.txt is non-systematic + 2007data to 320nm
c      OPEN(UNIT=kin,file=
c     & 'PHOTOCHEM/DATA/XSECTIONS/SO2/longward1.txt',STATUS='old')
c      n4 = 3025
c      endif

c longward2.txt is non-systematic + 2007data to 250nm + 2013data from 250-320
c      OPEN(UNIT=kin,file=
c     & 'PHOTOCHEM/DATA/XSECTIONS/SO2/longward2.txt',STATUS='old')
c      n4 = 4159


      endif




      DO i = 1, 1
         READ(kin,*)  !skip a line...
      ENDDO
      
      
      
      n5 = n4
      n6 = n4
      DO i = 1, n4
       if (newspec.eq.0) then
         READ(kin,*) xdum(i),x32(i),y32(i),ydum1(i),y33(i),ydum2(i),
     $               y34(i),ydum3(i)
       else
         READ(kin,*) x32(i),y32(i),y33(i), y34(i),ydum1(i),
     $               ydum2(i),ydum3(i)
       endif   

         x32(i)=x32(i)*10.   !convert from nm to A
         x33(i)=x32(i)
         x34(i)=x32(i)
      ENDDO
      CLOSE (kin)

      print *, 'using high resolution SO2 cross section'

      CALL addpnt(x32,y32,kdataHR,n4,x32(1)*(1.-deltax),zero)  
      CALL addpnt(x32,y32,kdataHR,n4,               zero,zero)
      CALL addpnt(x32,y32,kdataHR,n4,x32(n1)*(1.+deltax),zero)
      CALL addpnt(x32,y32,kdataHR,n4,            biggest,zero)
      CALL inter2(nw+1,wl,ygnew,n4,x32,y32,0)   


      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_SO2***'
         STOP
      ENDIF

!assuming quantum yield of 1 for <220nm and 0 otherwise

      endif  !end fine wavelength grid

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)    !SO2 + HV -> SO + O 
              sq(jn+1,i,iw) = yg2(iw)     !SO2 + HV -> SO21 
              sq(jn+2,i,iw) = yg3(iw)     !SO2 + HV -> SO23 
         ENDDO
      ENDDO
      endif   ! end option 1


      photolabel(jn)='PSO2'
      jn=jn+1

      photolabel(jn)='PSO2_SO21'
      jn=jn+1

      photolabel(jn)='PSO2_SO23'
      jn=jn+1

      RETURN
      END


       SUBROUTINE XS_O3(nw,wl,tlev,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for O3 photolysis   =*
*=         O3 + HV -> O2 + O(1D)                                             =*
*=         O3 + HV -> O2 + O(3P)                                             =*
*=  Cross section:  From Kevin's photo.dat - plus far UV data from MPI       =*  !UPDATE ME
*=  Quantum yield:  included in Kevin's cross sections                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn!,jdum
      REAL*8 wl(nw+1) ! EWS - wc not needed (O3)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=130)
      INTEGER n1, n2
      REAL*8 x1(kdata), x2(kdata)
      REAL*8 y1(kdata), y2(kdata)

* local
      REAL*8 yg1(nw), yg2(nw)
      REAL*8 o3xs(nz,nw),T273,T203,TI,FR,zero_9, RO3
      REAL*8 TAU,TAU2,TAU3, AT(nz),BT(nz),CT(nz),ALM0(nz),BL
      REAL*8 qy(nz,nw)
      INTEGER i, iw
      INTEGER ierr,option
      INTEGER L
      ierr = 0
**************** O3 photodissociation

C  O3 + HV -> O2 + O(1D)
C  O3 + HV -> O2 + O(3P)

c options
c     1)Kevin's photo.dat data - combined with far uv data from MPI
c - assuming no temperature dependence in the shortwave

      option=1

      if (option.eq.1) then  !Kevin's data

c$$$C below here is temp to bin the far uv data from MPI to our grid...
c$$$      OPEN(UNIT=kin,
c$$$     &  file='PHOTOCHEM/DATA/XSECTIONS/O3/O3_mpi_SW_Mason96.abs',
c$$$     &  STATUS='old')
c$$$c      DO i = 1, 2
c$$$c         READ(kin,*)
c$$$c      ENDDO
c$$$      n1 = 53
c$$$      DO i = 1, n1
c$$$         READ(kin,*) x1(i),y1(i)
c$$$         x1(i)=x1(i)*10.  !convert nm to A
c$$$      ENDDO
c$$$      CLOSE (kin)
c$$$      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
c$$$      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
c$$$      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
c$$$      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
c$$$      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   !inter2 - points to bins
c$$$
c$$$      do i=1,nw
c$$$         print *, wl(i),yg1(i)
c$$$      enddo
c$$$      stop

      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 
C         print*,"HJtest",HJtest

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O3/O3_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 = 118
      n2 = 118
      DO i = 1, n1
         READ(kin,*) x1(i),y1(i),y2(i)
         x2(i)=x1(i)
      ENDDO
      CLOSE (kin)


      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

!yg1 is "ozone 1" from photo.dat, i.e. SO31 in Kevin's code

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_O3***'
         STOP
      ENDIF


      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n2,               zero,zero)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,            biggest,zero)
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,ierr)
!yg2 is "ozone 2" from photo.dat i.e. SO32 in Kevin's original code

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_O3***'
         STOP
      ENDIF


c ok, the first thing kevin does is combine the two cross sections in the intermediate regime

!need to abstract L's, but should be OK for original testing since grids are identical...



      T273 = 273.0  ! needed for double precision
      T203 = 203.0  ! needed for double precision
     

      
      do L=1,nw        


      if (wl(L) .GE. 2632. .AND. wl(L) .LE. 3550.) then     !temperature dependent x-section from Kevin's code 
        do I=1,NZ
         TI = max(tlev(I),T203)
         TI = min(tlev(I),T273)
         FR = (TI - 203.)/70.
         O3xs(I,L) = FR*yg2(L) + (1.-FR)*yg1(L)
        enddo
      else
        do I=1,NZ
         O3xs(I,L) = yg1(L)
        enddo
      endif

      enddo

! so we now have the total temperature-dependent O3 cross section at all heights
! this is what the radiative transfer code wants...
!now, there is also a temperature-dependent quantum yield implemented in the longwave loop


      RO3 = 0.9                   !branching ratio used in shortwave loop
      zero_9 = 0.9

      do L=1,nw
C   THIS SUBROUTINE CALCULATES TEMPERATURE COEFFICIENTS USED TO FIND
C   THE O(1D) QUANTUM YIELD IN O3 PHOTOLYSIS.  (SEE JPL, 1983.)

      IF (wl(L).GE.2532.) THEN    !longwave loop

      do I=1,NZ
       TAU = tlev(I) - 230.
       TAU2 = TAU * TAU
       TAU3 = TAU2 * TAU
       AT(I) =   0.332 + 2.565E-4*TAU + 1.152E-5*TAU2 + 2.313E-8*TAU3
       BT(I) = - 0.575 + 5.590E-3*TAU - 1.439E-5*TAU2 - 3.270E-8*TAU3
       CT(I) =   0.466 + 8.883E-4*TAU - 3.546E-5*TAU2 + 3.519E-7*TAU3
       ALM0(I) = 308.2 + 4.4871E-2*TAU + 6.938E-5*TAU2 - 2.5452E-6*TAU3
       BL = BT(I)*(0.1*WL(L) - ALM0(I))
       RO3 = AT(I)*ATAN(BL) + CT(I)
       RO3 = max(RO3,zero)
       RO3 = min(RO3,zero_9)
       qy(I,L)=RO3                !quantum yield for O(1D)
      enddo

      ELSE     !shortwave loop
       do i=1,nz   
        qy(I,L)=RO3                          !quantum yield for O(1D)
       enddo 
      ENDIF    
      enddo  !loop over wavelength

      

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = o3xs(i,iw)*qy(i,iw)        !return O3 + HV -> O2 + O(1D) cross section * quantum yield
              sq(jn+1,i,iw) = o3xs(i,iw)*(1.0 - qy(i,iw))  !return O3 + HV -> O2 + O(3P) cross section * quantum yield
         ENDDO
      ENDDO
      endif   ! end option 1

c      print *, jn, photolabel

      photolabel(jn)='PO3_O1D'
      jn=jn+1

      photolabel(jn)='PO3'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_HCL(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_HCL(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HCL photolysis  =*
*=               HCL  +  HV  ->  H  +  CL     
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  Assuming 1 with no research                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

c      IMPLICIT NONE
c      INTEGER kin,kj,kw 
c      PARAMETER(kin=33,kj=33,kw=250)    !kin - file unit#  !kj=number of reactions defined - will need to abstract or figure out
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (HCL)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg(nw), yg1(nw)
      REAL*8 qy
      INTEGER i
      INTEGER ierr,option
      ierr = 0      

**************** HCL photodissociation
c options
c     1)Kevin's photo.dat data
c     2)JPL-06 recomendation (more or less the same as above, but over wider grid.
        !HCL absorbs strongly shortward of 135nm but not included in rec because not likely to be important in any high CO2 atmosphere...
      option=2


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/HCL/HCL_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 
C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  34
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_HCL***'
         STOP
      ENDIF

c      print *, x1
c      print *, y1
c      print *, ''
c      print *, yg1
c      stop

c Quantum yield (no reference nor research)

      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(j,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1



      if (option.eq.2) then  !JPL-06
      OPEN(UNIT=kin, 
     & file='PHOTOCHEM/DATA/XSECTIONS/HCL/HCL_JPL06.abs',STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  31
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.   !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg,n1,x1,y1,0)   !inter2 is grid-bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_HCL***'
         STOP
      ENDIF

c      print *, wl
c      print *, y1
c      print *, ''
c      print *, yg     
c       stop

c Quantum yield ( - there are some channel data for excited states of CL but whatever...)

      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 2

      photolabel(jn)='PHCL'
      jn=jn+1


      RETURN
      END



       !EWS - wc not used
c      SUBROUTINE XS_O2(nw,wl,wc,tlev,airlev,jn,sq,columndepth,zy,IO2)
       SUBROUTINE XS_O2(nw,wl,tlev,airlev,jn,sq,columndepth,zy,IO2)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for O2 photolysis   =*
*=   returns O2 + hv -> O + O(1D) then O2 + hv -> O + O
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  from kevin's code                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*=               this sub also contains the O2 column depth and zenith angle =*
*=               needed for calculating the Shumman-Runge correction         =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'      
      SAVE/PBLOK/
* input
      INTEGER jn,nw
      INTEGER IO2
      REAL*8 wl(nw+1) ! EWS - wc not needed here (O2)
      REAL*8 tlev(nz) ! EWS - used here
      REAL*8 airlev(nz) ! - EWS - used here

* weighting functions
      REAL*8 sq(kj,nz,kw)

      
* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata), x2(kdata)
      REAL*8 y1(kdata), y2(kdata)

      INTEGER kdata2
      PARAMETER(kdata2=100000)
      real*8 x5(kdata2),y5(kdata2)

* local
      REAL*8 yg1(nw),yg2(nw),yg5(nw)
c      REAL*8 qy, scale !EWS - not used
      INTEGER i
      INTEGER ierr,option
      real*8 PLOG(NZ),SIGL(NZ),BIGX(NZ),KA(17),SIG0(NZ,17),CL(NZ),KB(17)
      real*8 A(17,9), B(17,5),SRO2(NZ,17),TO2L(NZ)
      real*8 columndepth(KJ,NZ) 
      real*8 BK,C,SD,AM,e_19,zy,PI,ZYR,U0
      integer KMAX,L,K
      BK = 1.38E-16 !Boltzmann constant - in erg/K
      PI = 3.14159
      ZYR = ZY*PI/180. ! note ZY is passed in subroutine call ! - solar angle in radians
      U0 = COS(ZYR)
      AM = 1./U0
      !1 don't know why kevin had this in here. something about double precision     
      e_19=2e-19
      ierr =0 

**************** O2 photodissociation
c options
c     1)Kevin's photo.dat data
c     2) High resolution O2 cross section

      if(IO2.le.1) option=1
      if(IO2.eq.2) option=2 

      if (option.eq.1) then  !Kevin's data

c this option calculates the O2 cross section for the modern atmosphere using
c the Schumman-Runge corrections from Allen and Frederick 1982.  This cross section
c is completely overwritten in SUBROUTINES/Photo.f when IO2=1 (i.e most of the time).
c  This cross section changes based on pressure,temperature
c and zenith angle, and so needs to be recomputed each time if these change
c BUT - if I'm never going to use it, I shouldn't bother recomputing each time
c when I go to a time-dependent code - rather I should focus on what needs to change with
c the exponential sums...         

c for starters, I'm not going to make the S-R corrections work for a variable grid, although SO2HZ will...

C - actually, option 1 should only be called for the LR grid. it's easier that way.
c - eventually, i should enforce this somewhere in the main code.

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O2/O2D_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  10
      DO i = 1, n1
         READ(kin,*) x2(i), y2(i)
      ENDDO
      CLOSE (kin)
      
C      DO i = 1, n1
C         print*,'x2(i),wl(i),y2(i),yg2(i)',x2(i),wl(i),y2(i),yg2(i)
C      ENDDO
      
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      CALL addpnt(x2,y2,kdata,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n1,               zero,zero)
      CALL addpnt(x2,y2,kdata,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg2,n1,x2,y2,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg2,n1,x2,y2,ierr)   
      ENDIF
C      DO i = 1, n1
C         print*,'x2(i),wl(i),y2(i),yg2(i)',x2(i),wl(i),y2(i),yg2(i)
C      ENDDO

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_O2***'
         STOP
      endif

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O2/O2_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  34
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

C      DO i = 1, n1
C         print*,'x1(i),wl(i),y1(i),yg1(i)',x1(i),wl(i),y1(i),yg1(i)
C      ENDDO  
      
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF
            
C      DO i = 1, n1
C         print*,'x1(i),wl(i),y1(i),yg1(i)',x1(i),wl(i),y1(i),yg1(i)
C      ENDDO    

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_O2***'
         STOP
      ENDIF

c   so now, yg1 is SO2HZ
c      print *, yg1
c      stop

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O2/AllenFrederick1982.coeff',
     &  STATUS='old')

      DO i = 1, 4
         READ(kin,*)
      ENDDO
      DO L = 1,17
         READ(kin,*)  KA(L), (A(L,K),K=1,9)
      ENDDO
      CLOSE (kin)

C   REPEAT THIS SECTION ONLY IF PRESSURE AND TEMPERATURE VARY WITH TIME 
c (i.e. call this subroutine when it needs to be called...)

      do I=1,NZ
       PLOG(I) = LOG10(airlev(I)*BK*TLEV(I)/1.E3)   
      enddo

      do L=1,17  !ack hardcoded wavelength grid

        do I=1,NZ
         SIGL(I) = 0.
        enddo

        do I=1,NZ
         IF (L .LT. 15) BIGX(I) = PLOG(I) 
         IF (L .GE. 15) BIGX(I) = TLEV(I)
        enddo 

        KMAX = INT(KA(L)) !number of coefficients

        do K=1,KMAX
         do I=1,NZ
          SIGL(I) = SIGL(I) + A(L,K)*BIGX(I)**(K-1)
         enddo
        enddo 

        do I=1,NZ
          SIG0(I,L) = 10.**SIGL(I)
        enddo
     
      enddo  !end loop over L for O2


C            REPEAT THIS SECTION ONLY IF SOLAR ZENITH ANGLE OR O2 VARIES
C            WITH TIME
c (i.e. call this subroutine when it needs to be called...)

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O2/AllenFrederick1982.coeff',
     &  STATUS='old')

      DO i = 1, 26
         READ(kin,*)
      ENDDO
      DO L = 1,17  !ack
         READ(kin,*)  KB(L), (B(L,K),K=1,5)
      ENDDO
      CLOSE (kin)


 
        do I=1,NZ
         TO2L(I) = LOG10(columndepth(j+1,I))
        enddo

        do L=1,17     !note this will only work for Kevin's original grid given 17 elements of B and KB
                         !unless we abstract B and KB using inter3.  don't know how this would work...
                      !it's OK to leave as is as long as the corrections are applied to the proper wavelength
           do I=1,NZ     
            CL(I) = 0.
           enddo

          KMAX = INT(KB(L))   !number of coefficients

           do K=1,KMAX
            do I=1,NZ
             CL(I) = CL(I) + B(L,K)*TO2L(I)**(K-1)
            enddo
           enddo 

c           if (L.eq.1) print *, (CL(I),I=1,NZ)
c           print *, CL(80),AM**(-1.0*(10.**CL(80)))


! there is some weirdness here at the upper boundary that causes a 0 where (I presume)
! a 1 should be for L=14.  should look at allen and frederick for confirmation...
!jim's original code had identical behavoir          

C     compute Shumman-Runge correction for O2 cross section at each height over each wavelength
          do I=1,NZ
            C = 10.**CL(I)
            SD = SIG0(I,L) * AM**(-C)
            SRO2(I,L) = min(SD,e_19)  ! e_19 = 2e-19
c            if (L.eq.1) print *, SD,min(SD,e_19)            
c            if (L.eq.1) print *, C,AM**(-C),AM            
         enddo

       enddo  !end loop over wavelength
c       print *, ''
c            print *, (SRO2(I,1),I=1,NZ)

c           stop    !uncomment this and pritn CL() statement above to investigate L=14 wierdness..


!ACK - below will only work with Kevin's original grid...
!unless I calculate based on kevin's grid and then interpolate
!but would need to know how it works first - perhaps wait until
! I use correlated k-s
!note especially the L-10 in the SRO2 vector which works because kevin's
!original shortwave grid has 10 elements

C      DO i = 1, 10
C       print*,'(After Stuff) x2,wl,y2,yg2',x2(i),wl(i),y2(i),yg2(i)
C      ENDDO
C      DO i = 1, 34
C       print*,'(After Stuff) x1,wl,y1,yg1',x1(i),wl(i),y1(i),yg1(i)
C      ENDDO

      do L=1,nw
       do i=1,nz
        sq(jn,i,L) = yg2(L)      !for O1D ...
C           if(i.eq.1)print*,"L,jn,sq(jn,1,L)",L,jn,sq(jn,i,L)
        
        if (HJtest.EQ.1.) then
         sq(jn+1,I,L)=0.0 
        else
         sq(jn+1,I,L)=yg1(L)        !for O2

         if (wl(L) .GE. 1754. .AND. wl(L) .LE. 2041.) then   
            sq(jn+1,I,L)=sq(jn+1,I,L) + SRO2(I,L-10) !add in Schuman-Runge correction  
c            sq(j+1,I,L)=sq(j+1,I,L) + SRO2(I,L) !add in Schuman-Runge correction  
         endif
        endif
C           if(i.eq.1)print*,"L,jn+1,sq(jn+1,1,L)",L,jn+1,sq(jn+1,i,L)
       enddo
      enddo 

      endif  !end option 1



      if (option.eq.2) then  !high resolution O2 data, same data for O1D+O
         print *,'using high resolution O2 cross section'

!still use the old shortwave cross-section for O2 + hv -> O1D+O
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O2/O2D_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  10
      DO i = 1, n1
         READ(kin,*) x2(i), y2(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n1,               zero,zero)
      CALL addpnt(x2,y2,kdata,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg2,n1,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_O2***'
         STOP
      endif



!but use high resolution for the S-R bands...
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O2/Yoshino92.abs', 
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  94380
      DO i = 1, n1
         READ(kin,*) x5(i), y5(i)
         x5(i)=x5(i)*10.  !convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x5,y5,kdata2,n1,x5(1)*(1.-deltax),zero)  
      CALL addpnt(x5,y5,kdata2,n1,               zero,zero)
      CALL addpnt(x5,y5,kdata2,n1,x5(n1)*(1.+deltax),zero)
      CALL addpnt(x5,y5,kdata2,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg5,n1,x5,y5,0)   



***Yoshino goes from 1792-2025*** so,i need to put the rest of the grid back on in order to better compare these photo rates
claire grid is same as old to 1770-1786 and after 2248-2273
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/O2/O2_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  34
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)  

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_O2***'
         STOP
      ENDIF

      do i=1,nw
      if(wl(i).le.1770. .OR. wl(i).ge.2026.) yg5(i)=yg1(i)   
c      print *, wl(i),yg5(i),yg1(i)
      enddo
c      stop

 

      DO L = 1, nw
         DO i = 1, nz
        sq(jn,i,L) = yg2(L)      !for O1D ...
        sq(jn+1,I,L)=yg5(L)        !for O2
         ENDDO
      ENDDO
      endif  !end option 2

c      print *, jn, photolabel

      photolabel(jn)='PO2_O1D'
      jn=jn+1

      photolabel(jn)='PO2_O3P'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CH4(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CH4(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH4 photolysis  =*
*=                    CH4 + HV  ->  ^1CH2 + H2                               =*
*=                    CH4 + HV  ->  CH3 +  H                                 =*
*=                    CH4 + HV  ->  ^3CH2 + H + H                            =*
*=  Cross section:  from photo.dat                                           =* 
*=  Quantum yield:  1 for 1st reac, except at Ly a                           =* 
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'      
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CH4)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used


* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CH4 photodissociation
c options
c     1)Kevin's photo.dat data

      option=1
      
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest
      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH4/CH4_zahnle.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  10
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
      	CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF 


      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH4***'
         STOP
      ENDIF

c Quantum yields and reactions depend on wavelength and presence of hydrocarbons

      DO iw = 1, nw
         if (wl(iw).eq.1216) then  !ack hardcoded wavelength
c-mab: November 2016: Switched qy2 and qy3 values from original due to comparison with a previous person.
c-mab: A discussion within the GSFC group (following discrepancy with Kopparapu et al 2012 code) suggested they must've gotten switched somewhere accidentally? 
            qy=0.24
            qy2=0.51
            qy3=0.25
         else  !first reaction dominates everywhere but Ly a
          qy=1.0
          qy2=0.0
          qy3=0.0
         endif   

         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
              sq(jn+2,i,iw) = yg1(iw)*qy3
         ENDDO
      ENDDO
      endif  !end option 1

c      print *, jn, photolabel

      photolabel(jn)='PCH4_1CH2'
      jn=jn+1

      photolabel(jn)='PCH4_CH3'
      jn=jn+1

      photolabel(jn)='PCH4_3CH2'
      jn=jn+1

      RETURN
      END


       ! EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_C2H6(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_C2H6(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for C2H6 photolysis =*
*=           C2H6 + HV  ->  2 3CH2 + H2   (IN PLACE OF C2H2 AND C2H4)        =*
*=           C2H6 + HV  ->  CH4 + 1CH2                                       =*
*=                                                                           =*
*=  next are the reactions in Jim and Shawn's organic haze code              =*
*=  the first reac above is zeroed while the second is modified              =*
*=           C2H6 + HV  ->  C2H2 + H2 + H2                                   =*
*=           C2H6 + HV  ->  C2H4 + H + H                                     =*
*=           C2H6 + HV  ->  C2H4 + H2                                        =*
*=           C2H6 + HV  ->  CH3 + CH3                                        =*
*=                                                                           =*
*=  Cross section:  From Kevin's photo.dat, verified - featureless sw        =* 
*=  Quantum yield:  Assuming Kevin's vals(0.8/0.2) with no research          =*
*=                  using's shawn's values for the other four                =*
*=                  -shawn's values vary at Ly a                             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      CHARACTER*8 ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (C2H6)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** C2H6 photodissociation
c options
c     1)Kevin's photo.dat data

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/C2H6/C2H6_zahnle.abs',
     &  STATUS='old')
     
      HJtest=0 !HOT JUPITER TEST
      do i=1,nsp
C         print*,i,ISPEC(i)
         if (ISPEC(i).eq.'HE') HJtest=1  !temp solution to get 3 reactions for hot jupiters at Ly alpha
      enddo 

C         print*,"HJtest",HJtest

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  10
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      IF (HJtest.eq.1) THEN
        CALL inter3(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ELSE
      	CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
      ENDIF  

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_C2H6***'
         STOP
      ENDIF

c Quantum yield's depend on if there are higher order hydrocarbon and wavelength
      HCtest=0
      do i=1,nsp
         if (ISPEC(i).eq.'C2H4') HCtest=1
      enddo   

      DO iw = 1, nw
            if (wl(iw).eq.1216) then !wl dependent quantum yield just at ly a
              qy  = 0.0
              qy2 = 0.25
              qy3  = 0.25
              qy4 = 0.30
              qy5  = 0.12
              qy6  = 0.08
            else    !different quantum yields everywhere but ly a
              qy  = 0.0
              qy2 = 0.02
              qy3  = 0.27
              qy4 = 0.14
              qy5  = 0.56
              qy6  = 0.01
            endif   

!override the above if we aren't including hydrocarbon aerosols

      if (HCtest.eq.0) then
       qy  = 0.8
       qy2 = 0.2
c       print *, 'note to future self - XS_C2H6 not formally tested in'
c       print *, 'non hc case. should be fine, but check to be sure'
c       print *, 'upon seeing this note for the first time...'
      endif 


         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
              if (HCtest.eq.1) then
              sq(jn+2,i,iw) = yg1(iw)*qy3
              sq(jn+3,i,iw) = yg1(iw)*qy4
              sq(jn+4,i,iw) = yg1(iw)*qy5
              sq(jn+5,i,iw) = yg1(iw)*qy6
              endif   

         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      photolabel(jn)='PC2H6_1'
      jn=jn+1
      photolabel(jn)='PC2H6_2'
      jn=jn+1

      if (HCtest.eq.1.0) then
      photolabel(jn)='PC2H6_3'
      jn=jn+1
      photolabel(jn)='PC2H6_4'
      jn=jn+1
      photolabel(jn)='PC2H6_5'
      jn=jn+1
      photolabel(jn)='PC2H6_6'
      jn=jn+1
      endif

      RETURN
      END

       ! EWS - wl, wc, sq, and nw not used
c      SUBROUTINE XS_NO(nw,wl,wc,tlev,airlev,jn,sq,columndepth)
       SUBROUTINE XS_NO(tlev,airlev,jn,columndepth)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for NO photolysis   =*
*=   returns NO + hv -> N + O
*=  Cross section:  From Kevin's photo.dat                                   =*  !UPDATE ME
*=  Quantum yield:  from kevin's code                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/QBLOK.inc'
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER jn
c      REAL*8 wl(nw+1), wc(nw) ! EWS - not used
      REAL*8 tlev(nz) ! EWS - used here
      REAL*8 airlev(nz) ! - EWS - used here

* weighting functions
c      REAL*8 sq(kj,nz,kw) !EWS - not used

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
c      INTEGER n1, n2 ! EWS - not used
c      REAL*8 x1(kdata), x2(kdata), x3(kdata) ! EWS - not used
c      REAL*8 y1(kdata), y2(kdata) ! EWS - not used

* local
c      REAL*8 yg1(nw),yg2(nw) ! EWS - not used
c      REAL*8 qy ! EWS - not used
      INTEGER i
      INTEGER ierr,option
      real*8 columndepth(KJ,NZ),PLOG(NZ),signol(NZ) 
      real*8 BK,SD,AM,zy,PI,ZYR,U0
      integer L,K
      real*8 ANO(9,2),BNO(5,2),LLNO(35),CNO(NZ)
      ierr = 0
      BK = 1.38E-16               !Boltzmann constant - in erg/K
      PI = 3.14159
      ZYR = ZY*PI/180.            ! note ZY is passed in subroutine call - solar angle in radians
      U0 = COS(ZYR)
      AM = 1./U0
      e_15=1.e-15                 !1 don't know why kevin had this in here. something about double precision      


C   NO PREDISSOCIATION COEFFICIENTS (ALLEN AND FREDERICK, 1982)
C
      DATA ANO/-1.790868E+1, -1.924701E-1, -7.217717E-2, 5.648282E-2,
     2  4.569175E-2, 8.353572E-3, 3*0.,
     3  -1.654245E+1, 5.836899E-1, 3.449436E-1, 1.700653E-1,
     4  -3.324717E-2, -4.952424E-2, 1.579306E-2, 1.835462E-2,
     5  3.368125E-3/
C
      DATA BNO/7.836832E+3, -1.549880E+3, 1.148342E+2, -3.777754E+0,
     2  4.655696E-2, 1.297581E+4, -2.582981E+3, 1.927709E+2,
     3  -6.393008E+0, 7.949835E-2/
C
      DATA LLNO/3*0, 2*2, 3*0, 2*1, 25*0/





**************** NO photodissociation
c options
c     1)Kevin's photo.dat data
        
      option=1


      if (option.eq.1) then  !Kevin's data

C   REPEAT THIS SECTION ONLY IF PRESSURE AND TEMPERATURE VARY WITH TIME 
c (i.e. call this subroutine when it needs to be called...)

      do I=1,NZ
       PLOG(I) = LOG10(airlev(I)*BK*TLEV(I)/1.E3)   
      enddo

C
C          COEFFICIENTS FOR NITROUS OXIDE (NO)
        do L=1,2

          do I=1,NZ
           SIGNOL(I) = 0.
          enddo
 
          do K=1,9
           do I=1,NZ
            SIGNOL(I) = SIGNOL(I) + ANO(K,L)*PLOG(I)**(K-1)
           enddo
          enddo
 
          do I=1,NZ
           SIGNO0(I,L) = 10.**SIGNOL(I)
          enddo
        enddo  !end loop over L for NO


C            REPEAT THIS SECTION ONLY IF SOLAR ZENITH ANGLE OR O2 VARIES
C            WITH TIME
c (i.e. call this subroutine when it needs to be called...)


        do I=1,NZ
         TO2L(I) = LOG10(columndepth(2,I))  !ACK - hardcoded to O2 which is position 2 in photoreac
        enddo

C            COEFFICIENTS FOR NO
        do L=1,2   !again - this is OK provided the corrections are applied at the proper wavelengths

          do I=1,NZ
           CNO(I) = 0.
          enddo
 
          do K=1,5
           do I=1,NZ
            CNO(I) = CNO(I) + BNO(K,L)*TO2L(I)**(K-1)
           enddo
          enddo

          do I=1,NZ
           SD = SIGNO0(I,L) * AM**CNO(I)
           SIGNO(I,L) = min(SD,e_15) ! e_15 = 1e-15
          enddo

        enddo  !end loop over L for NO coefficients

!NOTE - we are not returning a cross section for NO - the photolysis rate
! is computed in the big wavelength loop in Photo.f
!this also means that NO is not counted in the absorbtion scheme

      endif  !end option 1

      photolabel(jn)='PNO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev and wc not used
c      SUBROUTINE XS_NO3(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_NO3(nw,wl,tlev,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for NO3 photolysis  =*
*=         NO3 + hv -> NO + O2                                               =*
*=                  -> NO2 + O
*=  Cross section:  From JPL-06 (via max plank                               =*
*=  Quantum yield:  T-dependent from Johnston et al. 1996 (JPL-06 rec)       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (NO3)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata,kdata2
      PARAMETER(kdata=300,kdata2=100)
      INTEGER n1, n1o
      REAL*8 x1(kdata), x2(kdata2),x2o(kdata2)
      REAL*8 y1(kdata), y2(kdata2),y3(kdata2),y4(kdata2),y5(kdata2)
      REAL*8 y6(kdata2), y7(kdata2)

* local
      REAL*8 yg1(nw),yg2(nw),yg3(nw),yg4(nw),yg5(nw)
      REAL*8 yg6(nw),yg7(nw)
      REAL*8 qy1(nw,nz),qy2(nw,nz)
      REAL*8 T1,T2,T3
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** NO3 photodissociation

c options
c     1)JPL-06 data

      option=1

      if (option.eq.1) then  !JPL data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/NO3/NO3.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  289
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   ! inter2 is grid - bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_NO3***'
         STOP
      ENDIF



c Quantum yield (is temperature dependent)

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/NO3/NO3QY.dat',
     &  STATUS='old')

      DO i = 1, 7
         READ(kin,*)    !skip header
      ENDDO
      n1o =  55
      DO i = 1, n1o
         READ(kin,*) x2o(i),y2(i),y3(i),y4(i),y5(i),y6(i),y7(i)
         x2o(i)=x2o(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

!interpolate the quantum yield data to working wavelength grid
!test to see if I interpolate each group if the numbers sum to 1000.
      n1=n1o
      x2=x2o
      CALL addpnt(x2,y2,kdata2,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata2,n1,               zero,zero)
      CALL addpnt(x2,y2,kdata2,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata2,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg2,n1,x2,y2,ierr)   ! inter2 is grid - bins
      x2=x2o
      n1=n1o
      CALL addpnt(x2,y3,kdata2,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y3,kdata2,n1,               zero,zero)
      CALL addpnt(x2,y3,kdata2,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y3,kdata2,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg3,n1,x2,y3,ierr)   ! inter2 is grid - bins
      x2=x2o
      n1=n1o
      CALL addpnt(x2,y4,kdata2,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y4,kdata2,n1,               zero,zero)
      CALL addpnt(x2,y4,kdata2,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y4,kdata2,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg4,n1,x2,y4,ierr)   ! inter2 is grid - bins
      x2=x2o
      n1=n1o
      CALL addpnt(x2,y5,kdata2,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y5,kdata2,n1,               zero,zero)
      CALL addpnt(x2,y5,kdata2,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y5,kdata2,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg5,n1,x2,y5,ierr)   ! inter2 is grid - bins
      x2=x2o
      n1=n1o
      CALL addpnt(x2,y6,kdata2,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y6,kdata2,n1,               zero,zero)
      CALL addpnt(x2,y6,kdata2,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y6,kdata2,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg6,n1,x2,y6,ierr)   ! inter2 is grid - bins
      x2=x2o
      n1=n1o
      CALL addpnt(x2,y7,kdata2,n1,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y7,kdata2,n1,               zero,zero)
      CALL addpnt(x2,y7,kdata2,n1,x2(n1)*(1.+deltax),zero)
      CALL addpnt(x2,y7,kdata2,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg7,n1,x2,y7,ierr)   ! inter2 is grid - bins
 
      do l=1,nw
c         print *, wl(l),yg2(l),yg3(l),yg4(l)
c         print *, wl(l),yg5(l),yg6(l),yg7(l)
      enddo   

      T1=298.
      T2=230.
      T3=190.

      do l=1,nw
       do i=1,nz
          
          if (tlev(i) .ge. T1) then
             qy1(l,i)=yg2(l)
             qy2(l,i)=yg5(l)
         else if (tlev(i) .lt. T1 .and. tlev(i) .ge. T2) then
             w1=(tlev(i)-T2)/(T1-T2)
             w2=1.-w1
             qy1(l,i)=w1*yg2(l)+w2*yg3(l)
             qy2(l,i)=w1*yg5(l)+w2*yg6(l)
          else if (tlev(i) .lt. T2 .and. tlev(i) .ge. T3) then
             w1=(tlev(i)-T3)/(T2-T3)
             w2=1.-w1
             qy1(l,i)=w1*yg3(l)+w2*yg4(l)
             qy2(l,i)=w1*yg6(l)+w2*yg7(l)
          else if (tlev(i) .lt. T3) then
             qy1(l,i)=yg4(l)
             qy2(l,i)=yg7(l)
          endif   
       enddo   
      enddo

      
      DO iw = 1, nw
c         print *, iw, wl(iw)
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy1(iw,i)/1000.             !for NO + O2
              sq(jn+1,i,iw) = yg1(iw)*qy2(iw,i)/1000.           !for NO2 + O
              
c              if (iw .eq. 94) print *, qy1(iw,i),qy2(iw,i)
         ENDDO
      ENDDO
      endif  !end option 1

c      print *, jn, photolabel

c      write(14,*) jn, photolabel      

      photolabel(jn)='PNO3_NO'
      jn=jn+1

      photolabel(jn)='PNO3_NO2'
      jn=jn+1

      RETURN
      END



       !EWS - airlev not used
c      SUBROUTINE XS_N2O(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_N2O(nw,wl,wc,tlev,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for N2O photolysis  =*
*=         N2O + hv -> N2 + O1D                                              =*
*=  Cross section:  From JPL-06 (via max plank)                              =*
*=  Quantum yield:  T-dependent formula from JPL-06 rec                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1), wc(nw) ! EWS - wc needed here (N2O)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      REAL*8 lambda
c      REAL*8 T1,T2,T3 !EWS - not used
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** N2O photodissociation

c options
c     1)JPL-06 data

      option=1

      if (option.eq.1) then  !JPL data

!just use formula below..

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/N2O/N2O_JPL06.abs',
     &  STATUS='old')

      DO i = 1, 2
         READ(kin,*)
      ENDDO
      n1 =  69
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   ! inter2 is grid - bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_N2O***'
         STOP
      ENDIF


c cross section is temperature dependent - see JPL-06

      A0 = 68.21023
      A1 = -4.071805
      A2 = 4.301146E-02
      A3 = -1.777846E-04
      A4 = 2.520672E-07

      B0 = 123.4014
      B1 = -2.116255
      B2 = 1.111572E-02
      B3 = -1.881058E-05

**** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
**** Ravishankara), so quantum yield of O(1D) is assumed to be unity
      qy = 1.


      DO iw = 1, nw
         lambda = wc(iw)/10. !convet to nm

         IF (lambda .GE. 173. .AND. lambda .LE. 240.) THEN
           DO iz = 1, nz
             t = MAX(194.,MIN(tlev(iz),320.))    !320 or 302K ??? - 320 from graces...
             A = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
             B = (((B3*lambda+B2)*lambda+B1)*lambda+B0)
             B = (t-300.)*EXP(B)
             sq(jn,iz,iw) = qy * EXP(A+B)
c             print *, A
c             print *, B
           ENDDO
         ELSE
           DO iz = 1, nz
             sq(jn,iz,iw) = 0.
           ENDDO
         ENDIF
      ENDDO

      endif  !end option 1

      photolabel(jn)='PN2O'
      jn=jn+1

      RETURN
      END


       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CLO(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CLO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CLO photolysis  =*
*=                CLO + HV  ->  CL + O1D                                     =*
*=                CLO + HV  ->  CL + O                                       =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1 for O1D < 265 nm  1 for O above                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CLO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CLO photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CLO/CLO_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  72
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLO***'
         STOP
      ENDIF

c Quantum yield (1 for O1D < 265 nm then 1 for O longward

      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
            if (wl(iw) .le. 2650.) then
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*0.0
            else
               sq(jn,i,iw) = yg1(iw)*0.0
               sq(jn+1,i,iw) = yg1(iw)*qy
            endif
  
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCLO_O1D'
      jn=jn+1

      photolabel(jn)='PCLO_O3P'
      jn=jn+1

      RETURN
      END


       !EWS airlev, tlev, and wc not needed
c      SUBROUTINE XS_HOCL(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_HOCL(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HOCL photolysis =*
*=                HOCL + HV  ->  OH + CL                                     =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  using 1 given that second branch (HCL + 0) is <2%        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (HOCL)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=150)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** HOCL photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/HOCL/HOCL_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  111
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_HOCL***'
         STOP
      ENDIF

c Quantum yield (JPL-06 (i.e <2% O and >95% OH) are cited so just using 1 branch 
c      (i.e. ignoring HOCL + HV -> HCL + O)

      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PHOCL'
      jn=jn+1

      RETURN
      END

       !EWS - airlev and wc not used
c      SUBROUTINE XS_CL2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CL2(nw,wl,tlev,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HOCL photolysis =*
*=                CL2 + HV  ->  CL + CL                                      =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1 (JPL-06)                                               =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CL2)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=350)
c      INTEGER n1, n2 ! EWS - not used
c      REAL*8 x1(kdata), x2(kdata), x3(kdata) !EWS - not used
c      REAL*8 y1(kdata), y2(kdata) !EWS - not used

* local
c      REAL*8 yg1(nw) !EWS - not used
      REAL*8 qy, a
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** HOCL photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06

c using temperature dependent formula from JPL-06

c Quantum yield is 1 (jpl-06)


      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
          if ((wl(iw) .ge. 2500.) .AND. (wl(iw) .le. 5500.)) then
            a=TANH(402.7/tlev(I))
               sq(jn,i,iw) = qy*1.0E-20*(a**0.5)*
     $           (27.3*exp(-99.0*a*(log(329.5/wl(iw)*10.))**2 ) + 
     $           0.932*exp(-91.5*a*(log(406.5/wl(iw)*10.))**2))
          else 
             sq(jn,i,iw)=0.0
          endif   
         ENDDO 
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCL2'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CLOO(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CLOO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CLOO photolysis =*
*=                CLOO + HV  ->  O + CLO                                     =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CLOO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=50)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CLOO photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CLOO/CLOO_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  31
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLOO***'
         STOP
      ENDIF

c Quantum yield is 1  (JPL-06)


      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCLOO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_OCLO(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_OCLO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for OCLO photolysis =*
*=                OCLO + HV  ->  O + CLO                                     =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (OCLO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=250)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** OCLO photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/OCLO/OCLO_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  225
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
         y1(i)=y1(i)*1E-20  !convert
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_OCLO***'
         STOP
      ENDIF

c Quantum yield is 1  (JPL-06)
c some uncertainty about temperature dependence of bands AND about the possibily
c of extra branches (CL + O2 and CLOO), but following rec and using 1 and 1 branch

      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='POCLO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CLONO(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CLONO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CLONO photolysis=*
*=                CLONO + HV  ->  CL + NO2                                   =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CLONO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=50)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CLONO photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CLONO/CLONO_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  34
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLONO***'
         STOP
      ENDIF

c Quantum yield is 1  (JPL-06)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCLONO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev not used
c      SUBROUTINE XS_CLONO2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CLONO2(nw,wl,wc,tlev,jn,sq)
*-------------------------------------------------------------------------------*
*=  PURPOSE:                                                                   =*
*=  Provide product of (cross section) x (quantum yield) for CLONO2 photolysis =*
*=                CLONO2 + HV  ->  CL + NO3                                    =*
*=                CLONO2 + HV  ->  CLO + NO2                                   =*
*=  Cross section:  From JPL-06                                                =* 
*=  Quantum yield:  temperature dependent from jpl-06                          =*
*-------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                         =*
*-------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1), wc(nw) ! EWS - wc needed here (CLONO2)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=150)
      INTEGER n1, n2,n3
      REAL*8 x1(kdata), x2(kdata), x3(kdata)
      REAL*8 y1(kdata), y2(kdata), y3(kdata)

* local
      REAL*8 yg1(nw),yg2(nw),yg3(nw)
      REAL*8 qy1,qy2
      INTEGER i, iw,iz
      INTEGER ierr,option
      ierr = 0      

**************** CLONO2 photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CLONO2/CLONO2_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n =  119
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
         x1(i)=x1(i)*10.       !convert nm to A
        y1(i) = y1(i) * 1E-20  !convert
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE (kin)

      n1=n
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLONO2***'
         STOP
      ENDIF

      n2 = n
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,          zero,zero)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,        biggest,zero)
      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLONO2***'
         STOP
      ENDIF

      n3 = n
      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,          zero,zero)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,        biggest,zero)
      CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLONO2***'
         STOP
      ENDIF



c Quantum yields  (JPL-06)


      DO iw = 1, nw-1

*** quantum yields (from jpl97)

         IF( wc(iw) .LT. 3080.) THEN
            qy1 = 0.6
         ELSEIF( (wc(iw) .GE. 3080.) .AND. (wc(iw) .LE. 3640.) ) THEN
            qy1 = 7.143e-3 * wc(iw)/10. - 1.6
         ELSEIF( wc(iw) .GT. 3640. ) THEN
            qy1 = 1.0
         ENDIF
         qy2 = 1. - qy1

* compute T-dependent cross section

         DO iz = 1, nz
            xs = yg1(iw)*( 1. +
     $           yg2(iw)*(tlev(iz)-296) +
     $           yg3(iw)*(tlev(iz)-296)*(tlev(iz)-296))
            sq(jn,iz,iw) = qy1 * xs
            sq(jn+1,iz,iw) = qy2 * xs

c            if(iz .eq. 1) write(*,333) wc(iw), xs
c 333        format(0pf8.3,1pe11.4)

         ENDDO
      ENDDO

      endif  !end option 1

      photolabel(jn)='PCLONO2_CL'
      jn=jn+1

      photolabel(jn)='PCLONO2_CLO'
      jn=jn+1


      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CLNO(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CLNO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CLNO photolysis =*
*=                CLNO + HV  ->  CL + NO                                     =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CLNO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=150)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CLNO photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CLNO/CLNO_JPL06.abs',
     &  STATUS='old')

      DO i = 1, 4
         READ(kin,*)
      ENDDO
      n1 =  140
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLNO***'
         STOP
      ENDIF

c Quantum yield is 1  (JPL-06)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCLNO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CLNO2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CLNO2(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CLNO2 photolysis =*
*=                CLNO2 + HV  ->  CL + NO2                                   =*
*=     (second branch to CLONO + O is <2% so ignored)                        =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CLNO2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=50)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CLNO2 photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CLNO2/CLNO2_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  26
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CLNO2***'
         STOP
      ENDIF

c Quantum yield is 1  (JPL-06)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCLNO2'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CHCLO(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CHCLO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CHCLO photolysis=*
*=                CHCLO + HV  ->  HCO + CL                                   =*
*=     (JPL06 calls this species COHCL - formyl chloride)                    =*
*=      (products assumed from Yuk's earth model                             =*
*=  Cross section:  From JPL-06                                              =* 
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CHCLO)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CHCLO photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CHCLO/CHCLO_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  68
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CHCLO***'
         STOP
      ENDIF

c Quantum yield is 1  (JPL-06)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCHCLO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev not used
c      SUBROUTINE XS_CH3CL(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CH3CL(nw,wl,wc,tlev,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH3CL photolysis=*
*=         CH3CL + hv -> CL + CH3  
*=           another branch to H + CH2CL exists in far UV, but no recomended =*
*=           cross section exists - ignoring for now...                      +*
*=  Cross section:  298K data from JPL-06 (via max plank)                    =*
*=                  T-dependent formula from JPL-06 rec                      =*
*=  Quantum yield:  needs some work ACK
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1),wc(nw) ! EWS - wc needed here
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=50)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)


* local
      REAL*8 yg1(nw)
      REAL*8 qy
      REAL*8 lambda
c      REAL*8 T1,T2,T3 !EWS - not used
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CH3CL photodissociation

c options
c     1)JPL-06 data

      option=1

      if (option.eq.1) then  !JPL data

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH3CL/CH3CL_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  32
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   ! inter2 is grid - bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH3CL***'
         STOP
      ENDIF


c cross section is temperature dependent - see JPL-06

      A0 = -299.80
      A1 = 5.1047
      A2 = -3.3630E-02
      A3 = 9.5805E-05
      A4 = -1.0135E-07

      B0 = -7.1727
      B1 = 1.4837E-01
      B2 = -1.1463E-03
      B3 = 3.9188E-06
      B4 = -4.9994E-09

**** there is a quantum yield for H formation in the far UV, but no JPL recomended
**** cross section.

      qy = 1.


      DO iw = 1, nw
         lambda = wc(iw)/10. !convet to nm

         IF (lambda .GE. 174. .AND. lambda .LE. 216.) THEN
           DO iz = 1, nz
             t = MAX(210.,MIN(tlev(iz),300.))    ! 210<T<300
             A = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
             B = (((B4*lambda+B3)*lambda+B2)*lambda+B1)*lambda+B0
             B = (t-273.)*B
             sq(jn,iz,iw) = qy * 10**(A+B)
           ENDDO
         ELSE
           DO iz = 1, nz
             sq(jn,iz,iw) = yg1(iw)*qy
           ENDDO
         ENDIF
      ENDDO


c      print *, (wc(iw),sq(j,1,iw),iw=1,nw)
c      stop

      endif  !end option 1

      photolabel(jn)='PCH3CL'
      jn=jn+1


      RETURN
      END

       !EWS - airlev not used
c      SUBROUTINE XS_CCL4(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CCL4(nw,wl,wc,tlev,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CCL4 photolysis =*
*=         CCL4 + hv -> CCL3 + CL                                            =*
*=  Cross section:  T-dependent formula from JPL-06 rec                      =*
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1), wc(nw) !EWS - wc needed here
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=50)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      REAL*8 lambda
c      REAL*8 T1,T2,T3 !EWS - not used
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CCL4 photodissociation

c options
c     1)JPL-06 data

      option=1

      if (option.eq.1) then  !JPL data

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CCL4/CCL4_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  44
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   ! inter2 is grid - bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CCL4***'
         STOP
      ENDIF


c cross section is temperature dependent - see JPL-06

      A0 = -37.104
      A1 = -5.8218E-01
      A2 = 9.9974E-03
      A3 = -4.6765E-05
      A4 = 6.8501E-08

      B0 = 1.0739
      B1 = -1.6275E-02
      B2 = 8.8141E-05
      B3 = -1.9811E-07
      B4 = 1.5022E-10

**** there is a quantum yield for CCL2 formation in the far UV, but no JPL recomended
**** cross section.

      qy = 1.


      DO iw = 1, nw
         lambda = wc(iw)/10.   !convet to nm

         IF (lambda .GE. 194. .AND. lambda .LE. 250.) THEN
           DO iz = 1, nz
             t = MAX(210.,MIN(tlev(iz),300.))    ! 210<T<300
             A = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
             B = (((B4*lambda+B3)*lambda+B2)*lambda+B1)*lambda+B0
             B = (t-273.)*B
             sq(jn,iz,iw) = qy * 10**(A+B)
           ENDDO
         ELSE
           DO iz = 1, nz
             sq(jn,iz,iw) = yg1(iw)*qy
           ENDDO
         ENDIF
      ENDDO

c      print *, (wc(iw),sq(j,1,iw),iw=1,nw)
c      stop

      endif  !end option 1

      photolabel(jn)='PCCL4'
      jn=jn+1
   
      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_COCL2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_COCL2(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for COCL2 photolysis=*
*=                COCL2 + HV  ->  CL + CL + CO                               =*
*=   (products assumed from Yuk's earth model/JPL rec yeild between 200-280  =*
*=  Cross section:  From JPL-06 (temp-dependent at wings but not included    =* 
*=  Quantum yield:  1                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (COCOL2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=60)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** COCL2 photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/COCL2/COCL2_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  53
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_COCL2***'
         STOP
      ENDIF

c Quantum yield is 1  (JPL-06)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCOCL2'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CL2O2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CL2O2(nw,wl,jn,sq)
*-------------------------------------------------------------------------------*
*=  PURPOSE:                                                                   =*
*=  Provide product of (cross section) x (quantum yield) for CL2O2 photolysis  =*
*=                CL2O2 + HV  ->  CL + CLOO                                    =*
*=                CL2O2 + HV  ->  CLO + CLO                                    =*
*=  Cross section:  From JPL-06       (CLOOCL                                  =* 
*=  Quantum yield:  1 for CL < 300 nm  0.9 for CL above                        =*
*-------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                         =*
*-------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CL2O2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=150)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CL2O2 photodissociation
c options
c     1)JPL-06 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CL2O2/CL2O2_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 2
c         READ(kin,*)
c      ENDDO
      n1 =  131
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CL2O2***'
         STOP
      ENDIF

c Quantum yield (1 for CL < 390 nm then 0.9 longward
cNOTE - the second branch is uncertain - I should try to vary this when tuning 


      qy=1.0
      qy2=0.9

      DO iw = 1, nw
         DO i = 1, nz
            if (wl(iw) .le. 3000.) then
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*0.0
            else
               sq(jn,i,iw) = yg1(iw)*qy2
               sq(jn+1,i,iw) = yg1(iw)*(1-qy2)
            endif
  
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCL2O2_CL'
      jn=jn+1

      photolabel(jn)='PCL2O2_CLO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CH3O2NO2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CH3O2NO2(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CH3O2NO2 photolysis =*
*=                CH3O2NO2 + HV -> CH3O2 + NO2                                   =*
*=  Cross section:  From IUPAC 2006 recomendation                                =* 
*=  Quantum yield:  assuming 1 without research                                  =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CH3O2NO2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=35)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CH3O2NO2 photodissociation
c options
c     1)IUPAC 2006 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH3O2NO2/CH3O2NO2_IUPAC2006.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  26
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH3O2NO2***'
         STOP
      ENDIF

c assuming Quantum yield is 1 without research
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCH3O2NO2'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CH3OCL(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CH3OCL(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CH3OCL photolysis   =*
*=                CH3OCL + HV -> CH3O + CL                                       =*
*=  Cross section:  From IUPAC 2008 recomendation                                =* 
*=  Quantum yield:  alt. channel to CH2O + HCL < 1% so using 1 for above (JPL-06)=*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CH3OCL)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=35)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CH3OCL photodissociation
c options
c     1)IUPAC 2008 data

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH3OCL/CH3OCL_IUPAC2008.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  27
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH3OCL***'
         STOP
      ENDIF

c assuming Quantum yield is 1 (other channel is negligble (plus we don't do CH2O anyway...)
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCH3OCL'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CH3OOH(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CH3OOH(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CH3OCL photolysis   =*
*=                CH3OOH + HV -> CH3O + OH                                       =*
*=  Cross section:  JPL-06                                                       =* 
*=  Quantum yield:  1 (JPL-06)                                                   =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CH3OOH)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=40)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CH3OOH photodissociation
c options
c     1)JPL-06 rec 

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH3OOH/CH3OOH_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  32
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH3OOH***'
         STOP
      ENDIF

c Quantum yield is 1 
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCH3OOH'
      jn=jn+1

      RETURN
      END
       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_HO2NO2(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_HO2NO2(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for HO2NO2 photolysis   =*
*=                HO2NO2 + HV -> HO2 + NO2                                       =*
*=                HO2NO2 + HV -> OH  + NO3                                       =*
*=  Cross section:  JPL-06                                                       =* 
*=  Quantum yield:  wl dependent (JPL-06)                                        =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (HO2NO2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=65)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** HO2NO2 photodissociation
c options
c     1)JPL-06 rec 

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/HO2NO2/HO2NO2_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  54
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_HO2NO2***'
         STOP
      ENDIF

c Quantum yield depends on wavelength

      DO iw = 1, nw
            if (wl(iw) .le. 2000.) qy=0.7  !QY for HO2+NO2 when < 200nm
            if (wl(iw) .gt. 2000.) qy=0.8  !QY for HO2+NO2 when > 200nm
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*(1-qy)
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PHO2NO2_NO2'
      jn=jn+1

      photolabel(jn)='PHO2NO2_NO3'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CL2O(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CL2O(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CL2O photolysis     =*
*=                CL2O + HV -> CL + CLO                                          =*
*=  Cross section:  JPL-06                                                       =* 
*=  Quantum yield:  1 (JPL-06) There may be more branches at shorter wavelengths =*
*=    but no rec. yuk uses just this branch, which is good enough for me.        =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CL2O)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=55)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CL2O photodissociation
c options
c     1)JPL-06 rec 

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CL2O/CL2O_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  45
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_Cl20***'
         STOP
      ENDIF

c Quantum yield is 1, ignoring potential other branches at short wavelengths... 
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCL2O'
      jn=jn+1

      RETURN
      END

       !EWS - airlev not used
c      SUBROUTINE XS_N2O5(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_N2O5(nw,wl,wc,tlev,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for N2O5 photolysis     =*
*=                N2O5 + HV -> NO3 + NO2                                         =*
*=                N2O5 + HV -> NO3 + NO + O   (O3P)                              =*
*=  Cross section:  combination of data and t-dependent parameterization JPL-06  =* 
*=  Quantum yield:  1 (JPL-06) There may be more branches at shorter wavelengths =*
*=    but no rec. previous recs had lots of O yeild at short wl, but 06 rec      =*
*=    calls this into question.  keeping here but zeroing out, as yuk did.       =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1), wc(nw) ! EWS - wc needed here (N2O5)
      REAL*8 tlev(nz) ! EWS - used here
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=150)
      INTEGER n1, n2,n3
      REAL*8 x1(kdata), x2(kdata), x3(kdata)
      REAL*8 y1(kdata), y2(kdata), y3(kdata),A(kdata),B(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      REAL*8 lambda
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** N2O5 photodissociation
c options
c     1)JPL-06 rec 

      option=1


      if (option.eq.1) then  !JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/N2O5/N2O5_JPL06_298.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  34    !after 258 nm, the xs becomes temperature dependent, so am not reading the whole file...
      DO i = 1, n1 
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

c open/read in temperature-dependent parameters
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/N2O5/N2O5td.abs',
     &  STATUS='old')

      n2=20
      n3=n2

      DO i = 1, 1
         READ(kin,*)
      ENDDO

      do i = 1,n2
         READ(kin,*) x2(i), y2(i),y3(i)
         x2(i)=x2(i)*10.  !convert nm to A
         x3(i)=x2(i)
      enddo   


      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_N2O5***'
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n2,               zero,zero)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,            biggest,zero)
      CALL inter2(nw+1,wl,A,n2,x2,y2,1)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_N2O5***'
         STOP
      ENDIF



      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)  
      CALL addpnt(x3,y3,kdata,n3,               zero,zero)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,            biggest,zero)
      CALL inter2(nw+1,wl,B,n3,x3,y3,1)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_N2O5***'
         STOP
      ENDIF

c ok now have data read into yg1, and A and B factors interpolated on to the grid



c      do i=1,nw
c       print *, wl(i), yg1(i),A(i),B(i)
c      enddo

c      stop
c see N2O5thoughts.txt for details
c also test the sq's below.
c - make sure there are 0's in the proper place.
c - also check 420 nm.

      qy = 1.  !for now.  it is likely more complicated (see N2O5thoughts.txt)


      DO iw = 1, nw
         lambda = wc(iw)/10. !convet to nm

         IF (lambda .GE. 260. .AND. lambda .LE. 420.) THEN
           DO iz = 1, nz
             t = MAX(233.,MIN(tlev(iz),295.))    ! param valid from 233<T<295
             sq(jn,iz,iw) = qy *10**(A(iw)+1000*B(iw)/t)
             sq(jn+1,iz,iw) = (1-qy)*10**(A(iw)+1000*B(iw)/t)
           ENDDO
         ELSE
           DO iz = 1, nz
             sq(jn,iz,iw) = yg1(iw)*qy
c             print *, jn
c             print *, iz
c             print *, iw
c             print *, qy
             sq(jn+1,iz,iw) = yg1(iw)*(1-qy)
           ENDDO
         ENDIF
      ENDDO
      endif  !end option 1

      photolabel(jn)='PN2O5_NO2'
      jn=jn+1

      photolabel(jn)='PN2O5_NO'
      jn=jn+1

      RETURN
      END

       !EWS - airlev, tlev, and wc not used
c      SUBROUTINE XS_CLO3(nw,wl,wc,tlev,airlev,jn,sq)
       SUBROUTINE XS_CLO3(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CLO3 photolysis     =*
*=                CLO3 + HV -> CLO + O2                                          =*
*=  Cross section:  taken by eye from Kopitzky et al 2002                        =* 
*=  Quantum yield:  assuming 1                                                   =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CLO3)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=55)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CLO3 photodissociation
c options
c     1)extraction of Kopitzky et al 2002 by Mark Claire 

      option=1


      if (option.eq.1) then  !
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CLO3/CLO3_Kopitzky_MC.abs',
     &  STATUS='old')

      DO i = 1, 1
         READ(kin,*)
      ENDDO
      n1 =  39
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
         y1(i)=y1(i)*1e-18  !convert to proper units in the cross section
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_Cl03***'
         STOP
      ENDIF

c assuming Quantum yield is 1, 
      qy=1.0
c      print *, yg1
c      stop
      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCLO3'
      jn=jn+1

      RETURN
      END


       SUBROUTINE XS_CL2O3(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CL2O3 photolysis    =*
*=                CL2O3 + HV -> CLO + CLOO                                       =*
*=  Cross section:  JPL-06                                                       =* 
*=  Quantum yield:  assuming 1                                                   =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CL2O3)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=30)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CL2O3 photodissociation
c options
c     1)JPL-06

      option=1


      if (option.eq.1) then  !
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CL2O3/CL2O3_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 1
c         READ(kin,*)
c      ENDDO
      n1 =  21
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_Cl203***'
         STOP
      ENDIF

c assuming Quantum yield is 1, 
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCL2O3'
      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_CL2O4(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CL2O4 photolysis    =*
*=                CL2O4 + HV -> CLOO + OCLO                                      =*
*=  Cross section:  JPL-06                                                       =* 
*=  Quantum yield:  assuming 1                                                   =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CL2O4)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=85)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CL2O4 photodissociation
c options
c     1)JPL-06

      option=1


      if (option.eq.1) then  !
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CL2O4/CL2O4_JPL06.abs',
     &  STATUS='old')

c      DO i = 1, 1
c         READ(kin,*)
c      ENDDO
      n1 =  76
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_Cl204***'
         STOP
      ENDIF

c assuming Quantum yield is 1, 
      qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
               sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCL2O4'
      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_CS2(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CS2 photolysis  =*
*=         CS2 + HV -> CS + S                                                =*
*=         CS2 + HV -> CS2X                                                  =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  included in cross section file. 0/1                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed (CS2)
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=1200)
      INTEGER n1, n2
      REAL*8 x1(kdata), x2(kdata)
      REAL*8 y1(kdata), y2(kdata)

* local
      REAL*8 yg1(nw), yg2(nw)
c      REAL*8 qy  !EWS - not used
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CS2 photodissociation
c options
c     1)from MPI database

      option=1

      if (option.eq.1) then  !MPI

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CS2/CS2_mpi.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*) !ignore header lines
      ENDDO
      n1 = 158
      DO i = 1, n1
         READ(kin,*) x1(i),y1(i)
         x1(i)=x1(i)*10.  !convert to Angstroms
      ENDDO
      CLOSE (kin)

      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CS2/CS2X_mpi.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n2 = 908
      DO i = 1, n2
         READ(kin,*) x2(i),y2(i)
         x2(i)=x2(i)*10.  !convert to Angstroms
      ENDDO
      CLOSE (kin)



      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   !inter2 is used for discrete points -> bins
!yg1 is CS2 cross section* quantum yield for CS + S



      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CS2***'
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata,n2,               zero,zero)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,            biggest,zero)
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,ierr)   !inter2 is used for discrete points -> bins
!yg2 is CS2 cross section* quantum yield for CS2X

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CS2***'
         STOP
      ENDIF

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)   !CS2 -> CS + S
              sq(jn+1,i,iw) = yg2(iw)     !CS2 -> CS2X 
         ENDDO
      ENDDO

      endif   ! end option 1

       photolabel(jn)='PCS2_CS'
       jn=jn+1
       photolabel(jn)='PCS2_CS2X'
       jn=jn+1


      RETURN
      END

       SUBROUTINE XS_CH3SH(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH3SH photo.    =*
*=           CH3SH + HV  ->  H + CH3S                                        =*
*=           CH3SH + HV  ->  HS + CH3                                        =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  Assuming Shawn's vals (0.93/0.07)                        =* <- Check me!
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=1000)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** CH3SH photodissociation
c options
c     1)MPI data

      option=1


      if (option.eq.1) then  !MPI
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH3SH/CH3SH_mpi.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  235
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i) * 10. ! <- To convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)  !inter2 is used for discrete points -> bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH3SH***'
         STOP
      ENDIF

c Quantum yield (fix this) - sddg 

      qy  = 0.93
      qy2 = 0.07

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
         ENDDO
      ENDDO
      endif  !end option 1

      photolabel(jn)='PCH3SH_H'
      jn=jn+1

      photolabel(jn)='PCH3SH_HS'
      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_simple(species,nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for generic photo.  =*
*=  Cross section:  MPI database format, read in units of nm and cm^-2       =*
*=  Quantum yield:  1 for now, but we could add in a read for this later.    =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  Jn      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel,plab
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      INTEGER n1
      REAL(kind=8),dimension(:),allocatable :: x1,y1

* local
      REAL*8 yg1(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      CHARACTER*90 XsecFile
      CHARACTER*8 species
      ierr = 0

**************** photodissociation

c options
c     1)MPI data file

      option=1
      write(XsecFile,'(a,a,a,a,a)') 
     &  'PHOTOCHEM/DATA/XSECTIONS/',
     &  trim(species),'/',trim(species),'_mpi.abs'

      if (option.eq.1) then  !MPI
      OPEN(UNIT=kin,file=XsecFile,Status='old')


c Count lines in input file, allocate size of x1 and y1 accordingly
      kdata=0
      DO i = 1, 100000
         READ(kin,*,end=10)
         kdata=kdata+1
      ENDDO
 10   n1=kdata-3 ! for three header lines
      kdata=kdata+4 ! for interpolates
      ALLOCATE (x1(kdata)) ! set length of x1
      ALLOCATE (y1(kdata)) ! set length of y1
      CLOSE(kin)
   

      OPEN(UNIT=kin,file=XsecFile,Status='old')
      DO i = 1,3
         READ(kin,*)
      ENDDO

      DO i = 1, kdata
         READ(kin,*,end=11) x1(i), y1(i)
         x1(i)=x1(i)*10.   !default MPI unit is nm - convert to Angstroms
      ENDDO
 11   CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)   
      CALL addpnt(x1,y1,kdata,n1,               zero,zero) 
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)

      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   !inter2 is used for discrete points -> bins


      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_simple*** '
         STOP
      ENDIF
      
c      if (species.eq.'CH2CO') then
c         do i=1,nw
c         print *, wl(i),yg1(i)
c         enddo
c         stop
c      endif   

c Quantum yield (this could be read in from the reactions.rx file, eventually)
c but for now, xs_simple is restricted to molecules with one channel.
      qy=1.00

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
         ENDDO
      ENDDO

      endif  !end option 1

      write(plab,'(a,a)') 'P',trim(species)
      photolabel(jn)=plab
      jn=jn+1

      DEALLOCATE(x1)
      DEALLOCATE(y1)

      RETURN
      END

       SUBROUTINE XS_CHOCHO(nw,wl,jn,sq)
*---------------------------------------------------------------------------------*
*=  PURPOSE:                                                                     =*
*=  Provide product of (cross section) x (quantum yield) for CHOCHO photolysis   =*
*=                CHOCHO + HV -> HCO + HCO                                       =*
*=                CHOCHO + HV -> H2  + CO  + CO                                  =*
*=                CHOCHO + HV -> H2CO + CO                                       =*
*=  Cross section:  JPL-06+ IUPAC-05 hybrid                                      =* 
*=    NOT GUARANTEED TO WORK. NOT THE BEST TEMPLATE FOR A NEW XS                 =*
*=  Quantum yield:  wl dependent table (JPL-06)                                  =*
*---------------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                           =*
*---------------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=305,kdata2=50)
      INTEGER n1,n2,n3,n4,n5
      REAL*8 x1(kdata), x2(kdata2),x3(kdata2),x4(kdata2),x5(kdata2)
      REAL*8 y1(kdata), y2(kdata2),y3(kdata2),y4(kdata2),y5(kdata2)

* local
      REAL*8 yg1(nw),yg2(nw),yg3(nw),yg4(nw),yg5(nw)
      REAL*8 qy
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0      

**************** CHOCHO photodissociation
c options
c     1)JPL-06 rec 

      option=1

      if (option.eq.1) then  !hybrid
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CHOCHO/CHOCHO_hybrid.abs',
     &  STATUS='old')

c      DO i = 1, 4
c         READ(kin,*)
c      ENDDO
      n1 =  297
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i)*10.  !convert nm to A
      ENDDO
      CLOSE (kin)
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CHOCHO***'
         STOP
      ENDIF


       
c Quantum yield depends on wavelength and is read in from a table in JPL-06
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CHOCHO/CHOCHO.qy',
     &  STATUS='old')

      DO i = 1, 8
         READ(kin,*)
      ENDDO
      n2 =  46
      n3=n2
      n4=n2
      n5=n2
      DO i = 1, n2
         READ(kin,*) x2(i),y2(i),y3(i),y4(i),y5(i)
         x2(i)=x2(i)*10.  !convert nm to A
         x3(i)=x2(i)
         x4(i)=x3(i)
         x5(i)=x4(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata2,n2,x2(1)*(1.-deltax),zero)  
      CALL addpnt(x2,y2,kdata2,n2,               zero,zero)
      CALL addpnt(x2,y2,kdata2,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata2,n2,            biggest,zero)
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CHOCHO***'
         STOP
      ENDIF

      CALL addpnt(x3,y3,kdata2,n3,x3(1)*(1.-deltax),zero)  
      CALL addpnt(x3,y3,kdata2,n3,               zero,zero)
      CALL addpnt(x3,y3,kdata2,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata2,n3,            biggest,zero)
      CALL inter2(nw+1,wl,yg3,n3,x3,y3,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CHOCHO***'
         STOP
      ENDIF


      CALL addpnt(x4,y4,kdata2,n4,x4(1)*(1.-deltax),zero)  
      CALL addpnt(x4,y4,kdata2,n4,               zero,zero)
      CALL addpnt(x4,y4,kdata2,n4,x4(n4)*(1.+deltax),zero)
      CALL addpnt(x4,y4,kdata2,n4,            biggest,zero)
      CALL inter2(nw+1,wl,yg4,n4,x4,y4,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CHOCHO***'
         STOP
      ENDIF

      CALL addpnt(x5,y5,kdata2,n5,x5(1)*(1.-deltax),zero)  
      CALL addpnt(x5,y5,kdata2,n5,               zero,zero)
      CALL addpnt(x5,y5,kdata2,n5,x5(n5)*(1.+deltax),zero)
      CALL addpnt(x5,y5,kdata2,n5,            biggest,zero)
      CALL inter2(nw+1,wl,yg5,n5,x5,y5,0)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CHOCHO***'
         STOP
      ENDIF



      do L=1,38
       print 100,wl(L),yg1(L),yg2(L),yg3(L),yg4(L),yg5(L)
      enddo
      stop

 100  format(3X, F8.1, 5(2X,1PE10.3))


c this is not finished, only used for interpolative purposes...



      DO iw = 1, nw
            if (wl(iw) .le. 2000.) qy=0.7  !QY for HO2+NO2 when < 200nm
            if (wl(iw) .gt. 2000.) qy=0.8  !QY for HO2+NO2 when > 200nm
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*(1-qy)  !ACK - will need to update if this is ever used
         ENDDO
      ENDDO
      endif  !end option 1
      
      photolabel(jn)='PCHOCHO_1'
      jn=jn+1

      photolabel(jn)='PCHOCHO_2'
      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_C2H2(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for C2H2 photolysis =*
*=           C2H2 + HV  ->  C2H + H                                          =*
*=           C2H2 + HV  ->  C2 + H2                                          =*
*=                                                                           =*
*=           also used for C2H4 (Jdum=1), C3H8 (Jdum=2), CH (Jdum=3)         =*
*=                                                                           =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  Assuming Jim's vals(0.06/0.16) with no research          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/

* input
      INTEGER nw,jn,jdum
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100000)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** C2H2 photodissociation
c options
c     1)MPI cross sections, Jim's branching ratios

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/C2H2/C2H2_mpi.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  61743
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i) * 10. ! <- To convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   !discrete points into bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_C2H2***'
         STOP
      ENDIF

c      do i=1,nw
c       print *, wl(i),yg1(i)
c      enddo   
c      stop

c Quantum yields
      DO iw = 1, nw
       if (wl(iw).ge.1216 .and. wl(iw).le.1754) then  !for C2H2 in Shawns "shortwave" loop,which has some absorbtion in lya where ours doesn't
          qy=0.3
          qy2=0.1
       else  !for 1754-2532 
        qy=0.06 
        qy2=0.10
       endif 
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      if (Jdum.eq.3) photolabel(jn)='PCH'
      if (Jdum.eq.0) photolabel(jn)='PC2H2_H'

      jn=jn+1

      if (Jdum.eq.0) then
         photolabel(jn)='PC2H2_H2'
         jn=jn+1
      endif   


      RETURN
      END
     

       SUBROUTINE XS_CH3CHO(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH3CHO photo.   =*
*=           CH3CHO + HV  ->  CH3 + HCO                                      =*
*=           CH3CHO + HV  ->  CH4 + CO                                       =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  Assuming Jim's vals (0.50/0.50)                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/

* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=150)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** CH3CHO photodissociation
c options
c     1)MPI data

      option=1


      if (option.eq.1) then  !mpi data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH3CHO/CH3CHO_mpi.abs',
     &  STATUS='old')

c      DO i = 1, 3
c         READ(kin,*)
c      ENDDO
      n1 =  105
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i) * 10. ! <- To convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH3CHO***'
         STOP
      ENDIF

c Quantum yield (fix this)

      qy  = 0.50
      qy2 = 0.50

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      photolabel(jn)='PCH3CHO_1'
      jn=jn+1

      photolabel(jn)='PCH3CHO_2'
      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_C3H6(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for C3H6 photolysis =*
*=           C3H6 + HV  ->  C2H2 + CH3 + H                                   =*
*=           C3H6 + HV  ->  CH2CCH2 + H2                                     =*
*=           C3H6 + HV  ->  C2H4 + CH23                                      =*
*=           C3H6 + HV  ->  C2H + CH4 + H                                    =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  Assuming Jim's vals (0.34/0.57/0.02/0.05)                =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/

* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** C3H6 photodissociation
c options
c     1)MPI data

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/C3H6/C3H6_mpi.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  18
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i) * 10. ! <- To convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_C3H6***'
         STOP
      ENDIF

c Quantum yield (fix this)

      qy  = 0.34
      qy2 = 0.57
      qy3 = 0.02
      qy4 = 0.05

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
              sq(jn+2,i,iw) = yg1(iw)*qy3
              sq(jn+3,i,iw) = yg1(iw)*qy4
         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      photolabel(jn)='PC3H6_1'
      jn=jn+1

      photolabel(jn)='PC3H6_2'
      jn=jn+1

      photolabel(jn)='PC3H6_3'
      jn=jn+1

      photolabel(jn)='PC3H6_4'
      jn=jn+1

      RETURN
      END


       SUBROUTINE XS_CH3C2H(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH3C2H photo.   =*
*=           CH3C2H + HV  ->  C3H3 + H                                       =*
*=           CH3C2H + HV  ->  C3H2 + H2                                      =*
*=           CH3C2H + HV  ->  CH3 + C2H                                      =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  Assuming Jim's vals (0.40/0.15/0.02)                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=1000)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** CH3C2H photodissociation
c options
c     1)MPI data

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH3C2H/CH3C2H_mpi.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  222
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i) * 10. ! <- To convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH3C2H***'
         STOP
      ENDIF
      

c Quantum yield (fix this)

      qy  = 0.40
      qy2 = 0.15
      qy3 = 0.02

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
              sq(jn+2,i,iw) = yg1(iw)*qy3
         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      photolabel(jn)='PCH3C2H_1'
      jn=jn+1

      photolabel(jn)='PCH3C2H_2'
      jn=jn+1

      photolabel(jn)='PCH3C2H_3'
      jn=jn+1

      RETURN
      END

       SUBROUTINE XS_CH2CCH2(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH2CCH2 photo.  =*
*=           CH2CCH2 + HV  ->  C3H3 + H                                      =*
*=           CH2CCH2 + HV  ->  C2H2 + CH23                                   =*
*=           CH2CCH2 + HV  ->  CH3 + C2H                                     =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  Assuming Jim's vals (0.40/0.15/0.06)                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/

* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** CH2CCH2 photodissociation
c options
c     1)MPI data

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/CH2CCH2/CH2CCH2_kasting.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  10
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i) * 10. ! <- To convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_CH2CCH2***'
         STOP
      ENDIF

c Quantum yield (fix this)

      qy  = 0.40
      qy2 = 0.15
      qy3 = 0.06

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
              sq(jn+2,i,iw) = yg1(iw)*qy3
         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      photolabel(jn)='PCH2CCH2_1'
      jn=jn+1

      photolabel(jn)='PCH2CCH2_2'
      jn=jn+1

      photolabel(jn)='PCH2CCH2_3'
      jn=jn+1

      RETURN
      END
 

       SUBROUTINE XS_C2H4(nw,wl,jn,sq,Jdum)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for C2H4 photolysis =*
*=           C2H4 + HV  ->  C2H2 + H2                                        =*
*=           C2H4 + HV  ->  C2H2 + H                                         =*
*=  Cross section:  From Jim's code                                          =*
*=  Quantum yield:  Assuming Jim's vals(0.49/0.51) with no research          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=100)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)
c
* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw,Jdum
      INTEGER ierr,option
      ierr = 0

**************** C2H4 photodissociation
c options
c     1)Kasting data

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/C2H4/C2H4_kasting.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  45
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   


      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_C2H4***'
         STOP
      ENDIF

c Quantum yield (fix this)

      qy  = 0.49
      qy2 = 0.51
      
      if (Jdum.eq.3) qy=1.0

      DO iw = 1, nw
         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
          if (Jdum.eq.0) sq(jn+1,i,iw) = yg1(iw)*qy2
         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      photolabel(jn)='PC2H4_H2'
      if (Jdum.eq.3) photolabel(jn)='PCH'
      jn=jn+1

      if (Jdum.eq.0) then
      photolabel(jn)='PC2H4_H+H'
      jn=jn+1
      endif

      RETURN
      END


       SUBROUTINE XS_C3H8(nw,wl,jn,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for C3H8 photolysis =*
*=           C3H8 + HV  ->  C3H6 + H2                                        =*
*=           C3H8 + HV  ->  C2H6 + CH21                                      =*
*=           C3H8 + HV  ->  C2H4 + CH4                                       =*
*=           C3H8 + HV  ->  C2H5 + CH3                                       =*
*=  Cross section:  From MPI database                                        =*
*=  Quantum yield:  Assuming Jim's vals(0.49/0.51) with no research          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:  see above subroutines                                       =*
*-----------------------------------------------------------------------------*
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      CHARACTER*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      SAVE/PBLOK/
* input
      INTEGER nw,jn
      REAL*8 wl(nw+1) ! EWS - wc not needed
c      REAL*8 tlev(nz) ! EWS - not used
c      REAL*8 airlev(nz) ! - EWS - not used

* weighting functions
      REAL*8 sq(kj,nz,kw)

* data arrays
      INTEGER kdata
      PARAMETER(kdata=200)
      INTEGER n1
      REAL*8 x1(kdata)
      REAL*8 y1(kdata)

* local
      REAL*8 yg1(nw)
      REAL*8 qy,qy2
      INTEGER i, iw
      INTEGER ierr,option
      ierr = 0

**************** C3H8 photodissociation
c options
c     1)MPI data

      option=1


      if (option.eq.1) then  !Kevin's data
      OPEN(UNIT=kin,
     &  file='PHOTOCHEM/DATA/XSECTIONS/C3H8/C3H8_mpi.abs',
     &  STATUS='old')

      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n1 =  175
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         x1(i)=x1(i) * 10. ! <- To convert from nm to Angstroms
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n1,x1,y1,ierr)   
!ok since bin to bin
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in XS_C3H8***'
         STOP
      ENDIF

c Quantum yields done in wl loop. diff at ly a then everywhere else

      DO iw = 1, nw

       if (wl(iw).eq.1216) then  !for C3H8 at Ly a
        qy1 = 0.33 
        qy2 = 0.08
        qy3 = 0.39
        qy4 = 0.20
       else  !c3h8 quantum yields everywhere but Ly a
        qy1 = 0.94 
        qy2 = 0.0
        qy3 = 0.0
        qy4 = 0.06
       endif   

         DO i = 1, nz
              sq(jn,i,iw) = yg1(iw)*qy
              sq(jn+1,i,iw) = yg1(iw)*qy2
              sq(jn+2,i,iw) = yg1(iw)*qy3
              sq(jn+3,i,iw) = yg1(iw)*qy4
         ENDDO
      ENDDO
      endif  !end option 1

!      print *, jn, photolabel

      photolabel(jn)='PC3H8_C3H6'
      jn=jn+1

      photolabel(jn)='PC3H8_C2H6'
      jn=jn+1

      photolabel(jn)='PC3H8_C2H4'
      jn=jn+1

      photolabel(jn)='PC3H8_C2H5'
      jn=jn+1

      RETURN
      END
