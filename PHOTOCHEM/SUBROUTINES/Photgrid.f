      SUBROUTINE PHOTGRID
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      real *8 skip
C ***** SET UP THE VERTICAL GRID ZS *****
c-mc for now, since we are invoking constant dz, we can choose dzgrid based on NZ
      !!! Modified below so that all possible NZs have a corresponding dzgrid - Eddie !!!
      if (NZ .le. 100) dzgrid = 1.0e5
      if (NZ .le. 100.AND.PLANET.EQ.'WASP12B') dzgrid = 1.2949E+07
      if (NZ .gt. 100 .and. NZ .lt. 300) dzgrid = 0.5e5
      if (NZ .ge. 300 .and. NZ .le. 640) dzgrid = 0.25e5
      if (NZ .gt. 640 .and. NZ .lt. 800) dzgrid = 0.125e5
      if (NZ .ge. 800) dzgrid = 0.0625e5
      !dzgrid is constant stepsize for the troposphere and stratosphere


      do I=1,NZ
       Z(I) = (I - 0.5)*dzgrid
      enddo

      DZ(1)=Z(1)+0.5*dzgrid    !for now, hacking in the DZ(1) to be the same as the others. should it be 1/2 this?
      do I=2,NZ             !just to see if its working, though I will leave it like this.
        DZ(I)=Z(I)-Z(I-1) 
      enddo  
      RETURN
      END


      SUBROUTINE gridw(nw,wl,wc,wu,LGRID)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the wavelength grid for all interpolations and radiative transfer =*
*=  calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.   =*
*=  No gaps are allowed within the wavelength grid.                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW  - INTEGER, number of wavelength grid _points_                     (O)=*
*=  WL  - REAL, vector carrying the lower limit of each wavel. interval   (O)=*
*=  WC  - REAL, vector carrying the center wavel of each wavel. interval  (O)=*
*=              (wc(i) = 0.5*(wl(i)+wu(i), i = 1..NW-1)                      =*
*=  WU  - REAL, vector carrying the upper limit of each wavel. interval   (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
C R. F. Esswein 020214 Change all REAL declarations to REAL*8
C R. Esswein 020221 Change file references to full path names.
C R. Esswein 030407 Choose extension to lower wavelengths.
C M. Claire  060802  Integrating into Kevin's code
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      character*60 string
 
* input:
      INTEGER LGRID
* output:

      REAL*8 wl(kw), wc(kw), wu(kw)
      INTEGER nw

* local:

      REAL*8 wincr
      INTEGER iw,I
      LOGICAL ok
      INTEGER idum
      REAL*8 dum
      INTEGER mopt


*_______________________________________________________________________
      DO iw = 1, kw-1
         wl(iw) = 0.
         wu(iw) = 0.
         wc(iw) = 0.
      ENDDO


**** chose wavelengths

* some pre-set options
*     mopt = 1    equal spacing
*     mopt = 2    Isaksen's grid
*     mopt = 3    combined Kockarts/Isaksen grid + Lyman-Alpha
*     mopt = 4    user-defined
*     mopt = 5    (Notes below)
! grid from Kevin's code,used in Zahnle.flx/.grid
! (This is also Allen Grid for S-R + old JPL grid)
C-mab The present stellar flux files for Hot Jupiters also use this grid.
*     mopt = 6    grid from Jim's climate code  (entirly a hack right now for interpolative purposes)
*     mopt = 7    Jim's old grid, but high resolution from 175-220

!note - before using, make sure that nw is the number of wavelengths to loop over and
!that wl(nw+1)=wu(nw) 
! this is needed for the interpolations to work correctly
      if (LGRID.eq.0) mopt = 5
      if (LGRID.eq.1) mopt = 7


      IF (mopt .EQ. 1) GO TO 1
      IF (mopt .EQ. 2) GO TO 2
      IF (mopt .EQ. 3) GO TO 3
      IF (mopt .EQ. 4) GO TO 4
      IF (mopt .EQ. 5) GO TO 5
      IF (mopt .EQ. 6) GO TO 6
      IF (mopt .EQ. 7) GO TO 7

 1    CONTINUE
c      nw = 140 + 1
      nw = 108
      wincr = 1.0
      DO 10, iw = 1, nw-1
         wl(iw) = 280. + wincr*FLOAT(iw-1)
        wu(iw) = wl(iw) + wincr
         wc(iw) = ( wl(iw) + wu(iw) )/2.
   10 CONTINUE

c         print *, wl, wl(nw),wu(nw-1)

      wl(nw) = wu(nw-1) !why is this here? some sort of check digit, i guess
      !081006 - the final value is needed for the interpolater
      


      GO TO 9

 2    CONTINUE
      nw = 0
c orig      OPEN(unit=kin,
c orig    &  file=PHRT(1:INDEX(PHRT,' ')-1)//'DATAE1/GRIDS/isaksen.grid',
c orig    &  status='old')

      OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/isaksen.grid',status='old')
      DO iw = 1, 2 
         READ(kin,*)
      ENDDO
      DO iw = 1, 130
c      DO iw = 1, 108 !Isaken grid will only work with this line until we redo the definition of WAVL(108) to WAVL(NW) or something similar
         nw = nw + 1
         READ(kin,*) idum, dum, wc(nw), wl(nw), wu(nw)
      ENDDO
      CLOSE(kin)
      nw = nw + 1
      wl(nw) = wu(nw-1)
      GO TO 9

*** grid for strat photolysis calculations, extended at short wavelengths

 3    CONTINUE
      nw = 1
* include Lyman-Alpha wavelengths ([120.,121.4],[121.4,121.9],[123.,-])

      wl(nw) = 120.0
      wu(nw) = 121.4
      wc(nw) = (wl(nw)+wu(nw))*0.5

      nw = nw+1
      wl(nw) = wu(nw-1)
      wu(nw) = 121.9
      wc(nw) = (wl(nw)+wu(nw))*0.5

      nw = nw+1
      wl(nw) = wu(nw-1)
      wu(nw) = 122.3
      wc(nw) = (wl(nw)+wu(nw))*0.5

      nw = nw+1
      wl(nw) = wu(nw-1)
      wu(nw) = 123.1
      wc(nw) = (wl(nw)+wu(nw))*0.5

      nw = nw+1
      wl(nw) = wu(nw-1)
      wu(nw) = 123.8
      wc(nw) = (wl(nw)+wu(nw))*0.5

      nw = nw+1
      wl(nw) = wu(nw-1)
      wu(nw) = 124.6
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 125.4
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 126.2
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 127.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 128.6
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 129.4
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 130.3
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 132.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 135.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 137.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 145.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 155.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 165.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)
      wu(nw) = 170.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1

      wl(nw) = wu(nw-1)

      OPEN(kin,file='PHOTOCHEM/DATA/GRIDS/kockarts.grid',status='old')
      DO iw = 1, 16
         nw = nw + 1
         READ(kin,*) wl(nw), wu(nw)
         wc(nw) = ( wl(nw) + wu(nw) ) / 2.
         IF (iw .eq. 1) THEN
             wu(nw-1) = wl(nw)
             wc(nw-1) = (wl(nw-1) + wu(nw-1))*0.5
         ENDIF
      ENDDO
      CLOSE(kin)

      OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/isaksen.grid',status='old')

      DO iw = 1, 2 + 10
         READ(kin,*)
      ENDDO
      DO iw = 11, 130   !note that this is too many elements for now...X
c      DO iw = 11, 70   !note that this is too many elements for now...X
         nw = nw + 1
         READ(kin,*) idum, dum, wc(nw), wl(nw), wu(nw)
      ENDDO
      CLOSE(kin)
      nw = nw + 1
      wl(nw) = wu(nw-1)
      GO TO 9

 4    CONTINUE

C define wavelength intervals of width 1 nm from 150 - 420 nm:
      nw = 1
      wl(1) = 150.
      DO iw = 151, 420
        wu(nw) = Float(iw)
        wc(nw) = (wl(nw) + wu(nw))/2.
        nw = nw+1
        wl(nw) = Float(iw)
      ENDDO
C define wavelength intervals of width 10 nm from 420 - 700 nm:
      DO iw = 430, 700, 10
        wu(nw) = Float(iw)
        wc(nw) = (wl(nw) + wu(nw))/2.
        nw = nw+1
        wl(nw) = Float(iw)
      ENDDO


c define wavelength intervals from 178 - 730 nm (17 bins):

c     nw=1
c     wl(nw) = 178.6 
c     wu(nw) = 180.0
c     wc(nw) = (wl(nw) + wu(nw))/2      

c     nw=nw+1
c     wl(nw) = 180.0
c     wu(nw) = 183.0
c     wc(nw) = (wl(nw) + wu(nw))/2     

c     nw=nw+1
c     wl(nw) = 183.0
c     wu(nw) = 187.0
c     wc(nw) = (wl(nw) + wu(nw))/2

c     nw=nw+1
c     wl(nw) = 187.0
c     wu(nw) = 192.0
c     wc(nw) = (wl(nw) + wu(nw))/2     

c     nw=nw+1
c     wl(nw) = 192.0
c     wu(nw) = 198.0
c     wc(nw) = (wl(nw) + wu(nw))/2

c     nw=nw+1
c     wl(nw) = 198.0
c     wu(nw) = 205.0
c     wc(nw) = (wl(nw) + wu(nw))/2

c     nw=nw+1
c     wl(nw) = 205.0
c     wu(nw) = 213.0
c     wc(nw) = (wl(nw) + wu(nw))/2     

c     nw=nw+1
c     wl(nw) = 213.0
c     wu(nw) = 222.0
c     wc(nw) = (wl(nw) + wu(nw))/2

c     nw=nw+1
c     wl(nw) = 222.0
c     wu(nw) = 231.0
c     wc(nw) = (wl(nw) + wu(nw))/2         

c     nw=nw+1
c     wl(nw) = 241.0
c     wu(nw) = 289.9
c     wc(nw) = (wl(nw) + wu(nw))/2     

c     nw=nw+1
c     wl(nw) = 289.9
c     wu(nw) = 305.5
c     wc(nw) = (wl(nw) + wu(nw))/2     

c     nw=nw+1
c     wl(nw) = 305.5
c     wu(nw) = 313.5
c     wc(nw) = (wl(nw) + wu(nw))/2     

c     nw=nw+1
c     wl(nw) = 313.5
c     wu(nw) = 337.5
c     wc(nw) = (wl(nw) + wu(nw))/2     

c     nw=nw+1
c     wl(nw) = 337.5
c     wu(nw) = 420.0
c     wc(nw) = (wl(nw) + wu(nw))/2    

c     nw=nw+1 
c     wl(nw) = 420.0
c     wu(nw) = 475.0
c     wc(nw) = (wl(nw) + wu(nw))/2    

c     nw=nw+1
c     wl(nw) = 475.
c     wu(nw) = 729.0
c     wc(nw) = (wl(nw) + wu(nw))/2    

c     nw=nw+1
c     wl(nw) = 729.0
c     wu(nw) = 743.6
c     wc(nw) = (wl(nw) + wu(nw))/2    
c     nw=nw+1

c     wl(nw) = 743.6
      GO TO 9

 5    CONTINUE 
c-mc read in Kevin's grid here
c-mab Presently, only grid option for hot Jup stellar fluxes.
      nw = 118
       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/zahnle.grid',status='old')

c      nw = 108
c      OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/zahnle.grid.orig',status='old')


      DO i = 1, 2
         READ(kin,*)  !skip header
      ENDDO

      do L=1,nw
       READ(kin,*) WL(L),WU(L)
       wc(L) = ( wl(L) + wu(L) )/2.
      enddo

      wl(nw+1) = wu(nw)  !final point for interpolative array

      CLOSE (kin)
      GO TO 9


 6    CONTINUE 
c-mc read in Jim's grid here ( this won't work to drive the code)
c-mc this is just used as a convienient interpolator

!solar grid
c$$$      nw = 38
c$$$       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/kastingGridSolar.dat',status='old')
c$$$
c$$$c      DO i = 1, 2
c$$$c         READ(kin,*)  !skip header
c$$$c      ENDDO
c$$$      do L=1,nw
c$$$       READ(kin,*) WL(L)
c$$$       wl(L)=wl(L)*1e4  !convert from microns to angstroms
c$$$      enddo
c$$$
c$$$      nw=37

!IR grid
      nw = 55
       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/kastingIRgrid.dat',
     &         status='old')

c      DO i = 1, 2
c         READ(kin,*)  !skip header
c      ENDDO
      do L=1,nw
       READ(kin,*) WU(L),WC(L)    !note - wu here is wavenumbers, wc is wavelength (in reverse order)
       wc(L)=wc(L)*1e4  !convert from microns to angstroms
      enddo

      do L=1,nw
         I = nw - L+1     
         wl(L)=wc(I)
      enddo   

      nw=54

      CLOSE (kin)
      GO TO 9

 7    CONTINUE 



c      nw = 556  !Kevin's grid here,but with 1Angstrom resolution between 175-220
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/claire.grid',status='old')
c      print *, 'using 1A grid between 175 and 220nm'


c - BELOW is the grid I used for most of the original highres SO2 calculations
c      nw = 1017  !Kevin's grid here,but with 0.5 Angstrom resolution between 175-220
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/0.5A.grid',status='old')
c      print *, 'using 0.5A grid between 175 and 220nm'


c      nw = 5228  !Kevin's grid here,but with 0.5 Angstrom resolution between 175-220 and 0.05A in S-R
c      print *,'using 0.5A grid between 175-220nm with 0.05A in SR bands'
c      OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/test2.dat',status='old')

c      nw = 12245  !Kevin's grid here,but with 0.5 Angstrom resolution between 175-220 and 0.02A in S-R
c      print *,'using 0.5A grid between 175-220nm with 0.02A in SR bands'
c      OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/test3.dat',status='old')

c THIS IS THE MAIN GRID FOR THE SULFUR MIF PAPER
      nw = 2889  !Kevin's grid here,but with 0.5 Angstrom resolution between 175-220 and 0.1A in S-R
      string='using 0.5A grid between 175-220nm with 0.1A in SR bands'

      OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/test.dat',status='old')


C This grid is for testing the photoexcitation bands...
c      nw = 4719  !Kevin's grid here,but with 0.5 Angstrom resolution between 175-320 and 0.1A in S-R
c      string='using 0.5A grid between 175-320nm with 0.1A in SR bands'
c      OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/testPE.dat',status='old')

c      nw = 534  !Kevin's grid here,but with 0.5 Angstrom resolution between 175-200
c       !this is a test to see if Halevy et al did something wrong. they sure did..
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/0.5A.gridTEST',status='old')
c      print *, 'using 0.5A grid between 175 and 200nm - HALEVY TEST'


c      nw = 1580  !Kevin's grid here,but with 1 Angstrom resolution between 175-330
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/long1A.grid',status='old')
c      print *, 'using 1A grid between 175 and 330nm'


c the low res grid - this is used only for testing the interpolation of the new Yosihino O2 xsection vs. JFK's exponential sums
c      nw = 118
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/zahnle.grid',status='old')
c      DO i = 1, 2
c         READ(kin,*)  !skip header
c      ENDDO

c the low res grid - with 2 extra points around ly a.  testing thuillier ly flux  (orig 10A to get total flux)
c      nw = 120
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/zahnle2.grid',status='old')
c      DO i = 1, 2
c         READ(kin,*)  !skip header
c      ENDDO



c 50 A grid - this is used only for testing the interpolation of the Thuilier et al fluxes versus Photo.dat
c      nw = 120
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/test4.dat',status='old')
c      print *, 'using a 50 A grid everywhere'


c using test grid between 380 and 382 nm just for testing...
c      nw = 4
c       OPEN(kin, file='PHOTOCHEM/DATA/GRIDS/gridtest',status='old')
c      print *, 'using test grid 380-382nm'

      print *, string
      write (14,*) ''
      write (14,*) string



      do L=1,nw
       READ(kin,*) WL(L),WU(L)
       wc(L) = ( wl(L) + wu(L) )/2.
      enddo

      wl(nw+1) = wu(nw)  !final point for interpolative array

      CLOSE (kin)
      GO TO 9


 9    CONTINUE


c-mc should probably print these out to output file rather than screen
      print *, 'NW = ',nw,'   WAVELENGTH GRID:',wl(1),' -',wu(nw),
     $  ' Angstroms'
      write(14,*)'NW = ',nw,'   WAVELENGTH GRID:',wl(1),' -',wu(nw),
     $  ' Angstroms'
* check grid for assorted improprieties:

      CALL gridck(kw,nw,wl,ok)

      IF (.NOT. ok) THEN
         print *,'STOP in GRIDW:  The w-grid does not make sense'
         STOP
      ENDIF

*_______________________________________________________________________

      RETURN
      END


      SUBROUTINE readflux(nw,wl,f,timega,msun)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read and re-grid extra-terrestrial flux data.                            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*=  timega - solar system age in Ga                                          =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax))                                               =*
*=  05/96  Put in different preset options                                   =*
C R. F. Esswein 020214 Change all REAL declarations to REAL*8
C    Use parameters "zero" and "largest" as
C      parameters to "addpnt" calls.
C R. Esswein 020221 Change file references to full path names.
c       M. Claire       091306 integrating into Kevin's code
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/QBLOK.inc'
      integer kdata
      parameter(kdata=26500) 

* input: (wavelength grid)
      INTEGER nw
      REAL*8 wl(kw)
      INTEGER iw

* output: (extra terrestrial solar flux)
      REAL*8 f(kw)

* INTERNAL:

* work arrays for input data files:
      REAL*8 x1(kdata), x2(kdata), x3(kdata)
      REAL*8 y1(kdata), y2(kdata), y3(kdata)
      INTEGER nhead, n, i, ierr

* data gridded onto wl(kw) grid:

      REAL*8 yg1(kw)
      REAL*8 yg2(kw)
      REAL*8 yg3(kw)
      REAL*8 hc
      PARAMETER(hc = 6.62E-34 * 2.998E8)

!      INTEGER msun
      REAL*8 refrac
      EXTERNAL refrac


*_______________________________________________________________________
* select desired extra-terrestrial solar irradiance, using msun:

c      msun = 12    !Gj 581 from Lucianne 
c      msun = 13    !high resolution solar data from ATLAS1/3 (Thullier et al 2004)
c      msun = 14    !kevin's data from photo.dat
c      msun = 15    !AD Leo from VPL climate DB
c      msun = 16    !AD LEO from VPL website
c      msun = 17    !Gj 581 from Lucianne 

!ACK - implementing YOUNGSUN.f into msun=13 for now, but it could apply to 14 if I ever wanted to use it again...
!but not for Mdwarfs...



      IF (msun .EQ. 13) THEN
         print *, "MSUN IS 13 (our sun)"
         call sleep(1)
         nhead = 0
         ierr = 0


!ATLAS 1 was closer to the peak of solar activity  (early 1992ish)
         n = 26150
         OPEN(UNIT=kin,file='PHOTOCHEM/DATA/FLUX/composite.atl1_1',
     &                 STATUS='old')
         print *, 'using Atlas 1 spectrum'

!ATLAS 3 was closer to the mean of solar activity  (late 1994ish)
c         n = 26182
c         OPEN(UNIT=kin,file='PHOTOCHEM/DATA/FLUX/composite.atl3',STATUS='old')

         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)  !this flux in mw/m2/nm, but is sampled at subangstom resolution
            x1(i)=x1(i)*10e0   !convert wavelength from nm to Angstoms
            x2(i)=x1(i)      ! x2 also angstroms
            x3(i)=x1(i)      ! x3 also angstroms
            y3(i)=y1(i)/10e0   ! for y3, convert thuillier flux to mw/m2/A
         ENDDO
         CLOSE (kin)

!iImplementing youngsun on Thuillier grid and apply to flux in energy space 
!BEFORE converting to photochemical units...

! I'm keeping x1,y1 as high res orig just to test integrated energy against the interpolations done with x3/y3

!so send youngsun x1 and n along with timega, relative flux returned in y2
         CALL youngsun(n,timega,x1,y2)

!OK have y2 which are flux multipliers on thuillier grid. apply to flux, then convert and done!

         do L=1,n
c            print *, x1(L),y3(L),y2(L),y3(L)*y2(L)
            y3(L)=y3(L)*y2(L)  !add in youngsun correction
         enddo
c         stop

         n3=n
         n2=n
         ierr=0

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,          zero,zero)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,        biggest,zero)
      CALL inter2(nw+1,wl,yg3,n3,x3,y3,ierr)  !inter2 is points to bins
!so yg3 is flux corrected for solar age on the model wavelength grid


      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,          zero,zero)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kdata,n2,        biggest,zero)
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,ierr)  !inter2 is points to bins
!yg 2 is youngsun correction only on the model wavelength grid



         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,'  Something wrong in Grid.f/readflux'
            STOP
         ENDIF         

c so yg3 is the average values for flux in mw/m2/A at any point on the working wavelength grid
c to convert to actual flux in the bin (mw/m2), need to multiply by bin width 

c-before we convert, some tests of energy conservation across the interpolation
c-note - this are conserved across various wavelength grids.
c         call cspint(n,x1,y1/10d0,2000d0,3550d0,y,e,work,poop)  
c         call cspint(n,x1,y1/10d0,3500d0,5000d0,y,e,work,poop)  
c         print *, poop  !ok this is the right answer for 200-350 (55.8 W/m2)  this is the thuillier flux

c         call cspint(nw,wl,yg3,2000d0,3550d0,y,e,work,poop)  
c         call cspint(nw,wl,yg3,3500d0,5000d0,y,e,work,poop)  
c         print *, poop !flux interpolated to new grid 

c         stop

c-convert to total energy in bin, then convert to photons/cm2/A (factor of 1/10 since formula derived for per nm)


c - NOTE:  One remaining problem here.  The ly A bin is now integrated over large interval and is not just ly a.  need to figure this out
c-modern day Ly a is about 5e11 in quantized units according to photo.dat

c- integrating Thuillier from 1210-1220A yields 5.6e11 photons/cm2/s (the "average" sun is 5e11).  We should just use this value.

c - the youngsun correction will be overwritten and should be done in energy space (rather than quantized state), but if done just at the wave center it's works fine

         indexLa=minloc(x1,1,x1.ge.1216)

         DO iw = 1, nw
            if (wl(iw).eq.1216) then
               f(iw)=5.6e11*y2(indexLa)
               relflux(iw)=y2(indexLa)
            else 
               f(iw) = yg3(iw)*(wl(iw+1)-wl(iw))*5.039e8*wl(iw)/10. !convert to photons/cm2/s
               relflux(iw)=yg2(iw)
            endif
c            print *, wl(iw), yg3(iw)
         ENDDO


        ! do i=1,nw
        ! print *, wl(i),yg3(i)
        ! enddo

c         do i=1,nw
c         print *, wl(i),f(i),relflux(i)
c         enddo
c         stop

      ENDIF  !end msun=13

      IF (msun .EQ. 14) THEN

         nhead = 2
         ierr = 0
c         n = 108
c         OPEN(UNIT=kin,file='PHOTOCHEM/DATA/FLUX/zahnle.flx.orig',STATUS='old')

         n = 118
         OPEN(UNIT=kin,file='PHOTOCHEM/DATA/FLUX/zahnle.flx',
     &                 STATUS='old')
         print *, 'using photo.dat solar flux'

         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, 10
            READ(kin,*) x1(i), y1(i)
            y1(i)=y1(i)*50.   !kevin's origin
         ENDDO
         DO i = 11, n
            READ(kin,*) x1(i), y1(i)
         ENDDO




         CLOSE (kin)
         
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n,               zero,zero)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n,            biggest,zero)
      CALL inter3(nw+1,wl,yg1,n,x1,y1,0)   !inter3 doesn't have any error checking at the moment
      !inter3 is bins to bins. 

         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,'  Something wrong in readflux.f'
            STOP
         ENDIF         

         DO iw = 1, nw
               f(iw) = yg1(iw)
         ENDDO
c         do i=1,nw
c         print *, wl(i),f(i)
c         enddo
c         stop

      ENDIF

      IF (msun .EQ. 15) THEN

         nhead = 1
         ierr = 0

         n = 4821
         OPEN(UNIT=kin,file='PHOTOCHEM/DATA/FLUX/dMV.flx',
     &                 STATUS='old')
         print *, 'using ADLeo spectrum - NEEDS CHECKING'

         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)   !wl in nm, flux in W/m2/nm
            x1(i)=x1(i)*10  !convert to angstroms 
            y1(i)=y1(i)*1e-4/1.98468e-16*x1(i)/10.     !convert to photons/cm2/s 
         ENDDO




         CLOSE (kin)
         
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n,               zero,zero)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n,            biggest,zero)
      CALL inter3(nw+1,wl,yg1,n,x1,y1,0)   !inter3 doesn't have any error checking at the moment
!ACK- I think this should be inter2

         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,'  Something wrong in readflux.f'
            STOP
         ENDIF         

         DO iw = 1, nw
               f(iw) = yg1(iw)
         ENDDO
c         do i=1,nw
c         print *, wl(i),f(i)
c         enddo
c         stop


      ENDIF



      IF (msun .EQ. 12) THEN

         nhead = 0
         ierr = 0

         n = 7649
         n = 7621  !removed three negative points near 3340 A that were a wavelength overlap
! also removed  74551.5, 75156.4, 75761.2, that were a wl overlap
!also duplicate points at:       192346.  9.42869e-05,       194040.  0.000102508       195734.  9.71178e-05,       197427.  9.13277e-05,       199121.  8.77170e-05,       199121.  9.64108e-05,       200814.  9.24294e-05,       200814.  8.79063e-05,       202508.  8.70496e-05
!      202508.  8.55883e-05,       204201.  8.55937e-05,       204201.  8.06153e-05,       205895.  7.24943e-05
!      205895.  8.30841e-05,       207588.  8.47565e-05,       207588.  7.70679e-05,       209282.  6.40397e-05
!      209282.  7.66739e-05,       210975.  6.94243e-05
!      210975.  7.56043e-05,       212669.  4.08846e-05
!      212669.  7.75429e-05,       214362.  6.88654e-05,       216056.  3.20911e-05


         OPEN(UNIT=kin,
     &      file='PHOTOCHEM/DATA/FLUX/gl581_scaled_revised.txt',
     &         STATUS='old')
         print *, 'using Gl581 spectrum'
         !NOTE: this is converting to the flux distance of GJ 581g -- use for THAT planet, not earthlike ones
         !unless you edit the stuff below...

         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)   !wl in Angstroms, flux in ergs/cm^2/s/A, at 1 AU
c ergs/cm2/s is equivalent to mw/m^2           
            if (y1(i).lt.0.0) y1(i)=0.0  !ignore negative fluxes...
            y1(i)=y1(i)*(1./0.146)**2 !convert to distance of Gj581g which is 0.146 AU
         ENDDO



         CLOSE (kin)
         
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n,               zero,zero)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n,            biggest,zero)
      CALL inter2(nw+1,wl,yg1,n,x1,y1,0)   !inter2 is points to bins

         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,'  Something wrong in readflux.f'
            STOP
         ENDIF         

c so yg1 is the average values for flux in mW/m2/A at any point on the working wavelength grid
c to convert to actual flux in the bin (mW/m2), need to multiply by bin width 

c - NOTE:  One remaining problem here.  The ly A bin is now integrated over large interval and is not just ly a.  need to figure this out
c- integrating Lucianne's spectrum (at 0.146 AU) from 1205-1225A yields 108.881 mW/m2/A or 2177.62 mW/m2 which converts to 
c 1.33432e+14 in photons/cm2/s 


         DO iw = 1, nw
            if (wl(iw).eq.1216) then
               f(iw)=1.33432e14    
            else 
               f(iw) = yg1(iw)*(wl(iw+1)-wl(iw))*5.039e8*wl(iw)/10.
            endif
         ENDDO

         close (kin)


      ENDIF


!~~~~~~~~~~~ GJ 876 added by Eddie as Used in Domagal-Goldman et al. 2014 ~~~~~~~~~~~~!
!~~~~~~~~~~~ Other stars incorporated by Giada into this same piece of code ~~~~~~~~~~!
!-------------the stars in this section all have the same input units of:-------------!
!-------wavelength = nm; flux = mW/m2/nm; they have the equivalent flux at 1 AU-------!
!-------such that you can keep the planets at 1 AU and have the correct flux----------!

      IF (msun .GT. 15.AND.msun.NE.22) THEN
c-mab: Recall msun = 22 is using Kevin's grid. We don't need the conversions below.. 
                       
      IF (msun .EQ. 76) THEN
         n = 26035
         print *, "MSUN IS 76!"
         call sleep(1)
         nhead = 0
         ierr = 0
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/GJ876_atlas_units_grid.dat',
     &                STATUS='old')
          print *, 'using GJ 876 Spectrum - Giddyup and Ride on Cowboy!'
      ENDIF

      IF (msun .EQ. 77) THEN
         n = 26035
         print *, "MSUN IS 76 - scaled to F star!"
         call sleep(1)
         nhead = 0
         ierr = 0
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/gj876_77_units.txt',
     &                STATUS='old')
          print *, 'using GJ 876 Spectrum - Giddyup and Ride on Cowboy!'
      ENDIF

      IF (msun .EQ. 78) THEN
         n = 26035
         print *, "MSUN IS 76 - scaled to AD Leo!"
         call sleep(1)
         nhead = 0
         ierr = 0
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/gj876_78_units.txt',
     &                STATUS='old')
          print *, 'using GJ 876 Spectrum - Giddyup and Ride on Cowboy!'
      ENDIF


      IF (msun .EQ. 16) THEN
         n = 26035
         print *, "MSUN IS 16 (AD Leo)!"
         call sleep(1)
         nhead = 0
         ierr = 0
         ! EWS Note: Fixed to spectrum added in Jan 2016 to include correct Lyman-Alpha Flux
         ! See Table 2 of Domagal-Goldman et al. 2014, ApJ, 792:90, and references therein 
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/adleo_dat_units_fixed.txt',
     &                STATUS='old') !
          print *, 'Suzanne Hawley would be proud of you!'
      ENDIF

      IF (msun .EQ. 20) THEN
         n = 26035
         print *, "MSUN IS 20 (AD Leo w/ UV flux scaled to GJ 876)!"
         call sleep(1)
         nhead = 0
         ierr = 0
         ! EWS Note: Fixed to spectrum added in Jan 2016 to include correct Lyman-Alpha Flux
         ! See Table 2 of Domagal-Goldman et al. 2014, ApJ, 792:90, and references therein 
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/adleo_20_units_fixed.txt',
     &                STATUS='old')
          print *, 'Suzanne Hawley would be proud of you!'
      ENDIF

      IF (msun .EQ. 17) THEN
         n = 26141
         print *, "MSUN IS 17 (T3200)!"
         call sleep(1)
         nhead = 0
         ierr = 0
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/T3200_units.txt',
     &                STATUS='old')
          print *, 'T3200 is a modeled inactive M dwarf.'
      ENDIF


      IF (msun .EQ. 18) THEN
         n = 26034
         print *, "MSUN IS 18 (K2V)!"
         call sleep(1)
         nhead = 0
         ierr = 0
         ! EWS Note: Fixed to spectrum added in Jan 2016 to include correct Lyman-Alpha Flux
         ! See Table 2 of Domagal-Goldman et al. 2014, ApJ, 792:90, and references therein 
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/K2V_units_fixed.txt',
     &                STATUS='old')
          print *, 'K, you have a K star.'
      ENDIF

      IF (msun .EQ. 19) THEN
         n = 26035
         print *, "MSUN IS 19 (F2V)!"
         call sleep(1)
         nhead = 0
         ierr = 0
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/F2V_units.txt',
     &                STATUS='old')
          print *, 'An F star a day keeps the doctor away'
      ENDIF

      IF (msun .EQ. 20) THEN
         n = 26150
         print *, "MSUN IS 20 (M8V)!"
         call sleep(1)
         nhead = 0
         ierr = 0
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/M8_active_photogrid.txt',
     &                STATUS='old')
          print *, 'Here is an M8 dwarf!'
      ENDIF


      IF (msun .EQ. 21) THEN
         n = 26024
         print *, "MSUN IS 21 (Proxima Centauri)!"
         call sleep(1)
         nhead = 0
         ierr = 0
         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/Proxima_units_rev1.txt',
     &                STATUS='old')
          print *, 'Proxima Centauri! Yay for nearby red stars!'
      ENDIF



        ! Do same unit conversions as found in msun=13 option - Eddie
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)  !this flux in mw/m2/nm, but is sampled at subangstom resolution
            x1(i)=x1(i)*10e0   !convert wavelength from nm to Angstoms
            x2(i)=x1(i)      ! x2 also angstroms
            x3(i)=x1(i)      ! x3 also angstroms
            y3(i)=y1(i)/10e0   ! for y3, convert thuillier flux to mw/m2/A
         ENDDO
         CLOSE (kin)

! We are not bothering with Youngsun stuff because it is a different star
         n3=n
         ierr=0

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,          zero,zero)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kdata,n3,        biggest,zero)
      CALL inter2(nw+1,wl,yg3,n3,x3,y3,ierr)  !inter2 is points to bins
!so yg3 is flux on the model wavelength grid

         !error check for call to inter2
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,'  Something wrong in Grid.f/readflux'
            STOP
         ENDIF     

! NOTE: This explicitly assumes the correct Lyman-Alpha flux *at the planet* has been included 
!      in the input spectrum. It would benefit you to ensure that this is really the case! 
        DO iw = 1, nw
               f(iw) = yg3(iw)*(wl(iw+1)-wl(iw))*5.039e8*wl(iw)/10. !convert to photons/cm2/s
               print *, wl(iw), yg3(iw)
        ENDDO


      ENDIF  !msun != 13
!c-mab msun = 22 added below for the GOV spectra used in wasp12b Hot Jupiter run...
      IF (msun .EQ. 22) THEN

         nhead = 2
         ierr = 0

         n = 118

         OPEN(UNIT=kin,
     &    file='PHOTOCHEM/DATA/FLUX/totalwasp12_G0V_forAtmos.dat',
     &                STATUS='old')
       
         print *, 'Using Ravis GOV spectrum for WASP12b'
c-mab: Ravis wasp12b star fluxes are already converted to photons/cm^2/s, distance unsure
c-mab: to convert from ergs to photons to create above file wc was used not wl
c-mab: is over grid identical to zanhle.grid, so gridding done with LGRID = 0 earlier
                  
         DO i = 1, nhead
            READ(kin,*)
         ENDDO

         DO i = 1, nw
            READ(kin,*) x1(i), skip, y1(i)  
c-mab: wavl,wavu,flux given  
c-mab: flux is already earth-equivalent in photons/s/cm^2/bin       
            if (y1(i).lt.0.0) y1(i)=0.0  !ignore negative fluxes...
         ENDDO
         
             CLOSE (kin)
                     
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kdata,n,               zero,zero)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n,            biggest,zero)
      CALL inter3(nw+1,wl,yg1,n,x1,y1,0)   
!inter3 doesn't have any error checking at the moment
      !inter3 is bins to bins. 

         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr,'  Something wrong in readflux.f'
            STOP
         ENDIF         

         DO iw = 1, nw
               f(iw) = yg1(iw)
         ENDDO
c         do i=1,nw
c         print *, wl(i),x1(i),y1(i),f(i)
c         enddo
c         stop

      ENDIF
      
*_______________________________________________________________________

      RETURN
      END


      SUBROUTINE gridck(k,n,x,ok)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Check a grid X for various improperties.  The values in X have to comply =*
*=  with the following rules:                                                =*
*=  1) Number of actual points cannot exceed declared length of X            =*
*=  2) Number of actual points has to be greater than or equal to 2          =*
*=  3) X-values must be non-negative                                         =*
*=  4) X-values must be unique                                               =*
*=  5) X-values must be in ascending order                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  K  - INTEGER, length of X as declared in the calling program          (I)=*
*=  N  - INTEGER, number of actual points in X                            (I)=*
*=  X  - REAL, vector (grid) to be checked                                (I)=*
*=  OK - LOGICAL, .TRUE. -> X agrees with rules 1)-5)                     (O)=*
*=                .FALSE.-> X violates at least one of 1)-5)                 =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
C R. F. Esswein 020214 Change all REAL declarations to REAL*8
C       M.C.            060802  Integrated into Kevin's code
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER k, n
      REAL*8 x(k)

* output:
      LOGICAL ok

* local:
      INTEGER i
*_______________________________________________________________________

      ok = .TRUE.

* check if dimension meaningful and within bounds

      IF (n .GT. k) THEN
         ok = .false.
         print *,'Number of data exceeds dimension'
         print *, k,n
         RETURN
       ENDIF         

      IF (n .LT. 2) THEN
         ok = .FALSE.
         print *, 'Too few data, number of data points must be >= 2'
         RETURN
      ENDIF

* disallow negative grid values

      IF(x(1) .LT. 0.) THEN
         ok = .FALSE.
         print *,'Grid cannot start below zero'
         RETURN
      ENDIF

* check sorting

      DO 10, i = 2, n
         IF( x(i) .LE. x(i-1)) THEN
            ok = .FALSE.
            print *,'Grid is not sorted or contains multiple values'
            print *, i, x(i),x(i-1)
            RETURN
         ENDIF
   10 CONTINUE
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew )

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
*=  ascending order                                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
*=  Y    - REAL vector of length LD, y-values                            (IO)=*
*=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
*=         program                                                           =*
*=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
*=         N < LD.  On exit, N is incremented by 1.                          =*
*=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
*=  YNEW - REAL, y-value of point to be added                             (I)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

C calling parameters

      INTEGER ld, n
      REAL*8 x(ld), y(ld)
      REAL*8 xnew, ynew
      INTEGER ierr

C local variables

      INTEGER insert
      INTEGER i

C-----------------------------------------------------------------------

* initialize error flag

      ierr = 0

* check n<ld to make sure x will hold another point

      IF (n .GE. ld) THEN
         WRITE(0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
         WRITE(0,*) '                        All elements used.'
         STOP
      ENDIF

      insert = 1
      i = 2

* check, whether x is already sorted.
* also, use this loop to find the point at which xnew needs to be inserted
* into vector x, if x is sorted.

 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
           WRITE(0,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//
     >                'in ascending order!', i,x(i),x(i-1)
           STOP
        ELSE
           IF (xnew .GT. x(i)) insert = i + 1
        ENDIF
        i = i+1
        GOTO 10
      ENDIF

* if <xnew,ynew> needs to be appended at the end, just do so,
* otherwise, insert <xnew,ynew> at position INSERT

      IF ( xnew .GT. x(n) ) THEN
 
         x(n+1) = xnew
         y(n+1) = ynew
  
      ELSE

* shift all existing points one index up

         DO i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         ENDDO

* insert new point

         x(insert) = xnew
         y(insert) = ynew
  
      ENDIF

* increase total number of elements in x, y

      n = n+1

      END

      SUBROUTINE inter1(ng,xg,yg, n,x,y)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on single, discrete points, onto a discrete target  =*
*=  grid.                                                                    =*
*=  The original input data are given on single, discrete points of an       =*
*=  arbitrary grid and are being linearly interpolated onto a specified      =*
*=  discrete target grid.  A typical example would be the re-gridding of a   =*
*=  given data set for the vertical temperature profile to match the speci-  =*
*=  fied altitude grid.                                                      =*
*=  Some caution should be used near the end points of the grids.  If the    =*
*=  input data set does not span the range of the target grid, the remaining =*
*=  points will be set to zero, as extrapolation is not permitted.           =*
*=  If the input data does not encompass the target grid, use ADDPNT to      =*
*=  expand the input array.                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG  - INTEGER, number of points in the target grid                    (I)=*
*=  XG  - REAL, target grid (e.g. altitude grid)                          (I)=*
*=  YG  - REAL, y-data re-gridded onto XG                                 (O)=*
*=  N   - INTEGER, number of points in the input data set                 (I)=*
*=  X   - REAL, grid on which input data are defined                      (I)=*
*=  Y   - REAL, input y-data                                              (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/95  Loop 10 restructured                                              =*
C R. F. Esswein 020214 Change all REAL declarations to REAL*8
C $Id$
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER n, ng
      REAL*8 xg(ng)
      REAL*8 x(n), y(n)

* output:
      REAL*8 yg(ng)

* local:
      REAL*8 slope
      INTEGER jsave, i, j
*_______________________________________________________________________

      jsave = 1

      DO 20, i = 1, ng
         yg(i) = 0.
         j = jsave
   10    CONTINUE
            IF ((x(j) .GT. xg(i)) .OR. (xg(i) .GE. x(j+1))) THEN
               j = j+1
            IF (j .LE. n-1) GOTO 10
*        ---- end of loop 10 ----
            ELSE
               slope = (y(j+1)-y(j)) / (x(j+1)-x(j))
               yg(i) = y(j) + slope * (xg(i) - x(j))
               jsave = j
             ENDIF
   20 CONTINUE
*_______________________________________________________________________


      RETURN
      END



      SUBROUTINE inter2(ng,xg,yg,n,x,y,ierr)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on single, discrete points onto a set of target     =*
*=  bins.                                                                    =*
*=  The original input data are given on single, discrete points of an       =*
*=  arbitrary grid and are being linearly interpolated onto a specified set  =*
*=  of target bins.  In general, this is the case for most of the weighting  =*
*=  functions (action spectra, molecular cross section, and quantum yield    =*
*=  data), which have to be matched onto the specified wavelength intervals. =*
*=  The average value in each target bin is found by averaging the trapezoi- =*
*=  dal area underneath the input data curve (constructed by linearly connec-=*
*=  ting the discrete input values).                                         =*
*=  Some caution should be used near the endpoints of the grids.  If the     =*
*=  input data set does not span the range of the target grid, an error      =*
*=  message is printed and the execution is stopped, as extrapolation of the =*
*=  data is not permitted.                                                   =*
*=  If the input data does not encompass the target grid, use ADDPNT to      =*
*=  expand the input array.                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
*=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
*=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
*=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
*=        bin i (i = 1..NG-1)                                                =*
*=  N   - INTEGER, number of points in input grid                         (I)=*
*=  X   - REAL, grid on which input data are defined                      (I)=*
*=  Y   - REAL, input y-data                                              (I)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER ng, n
      REAL*8 x(n), y(n), xg(ng)
      INTEGER ierr
* output:
      REAL*8 yg(ng)

* local:
      REAL*8 area, xgl, xgu
      REAL*8 darea, slope
      REAL*8 a1, a2, b1, b2
      INTEGER ngintv
      INTEGER i, k, jstart


*_______________________________________________________________________

*  test for correct ordering of data, by increasing value of x

      DO 10, i = 2, n
         IF (x(i) .LE. x(i-1)) THEN
            ierr = 1
            WRITE(*,*)'data not sorted', i,x(i),x(i-1)
            RETURN
         ENDIF
   10 CONTINUE     


      DO i = 2, ng
        IF (xg(i) .LE. xg(i-1)) THEN
           ierr = 2
          WRITE(0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
          RETURN
        ENDIF
      ENDDO




* check for xg-values outside the x-range

      IF ( (x(1) .GT. xg(1)) .OR. (x(n) .LT. xg(ng)) ) THEN
          WRITE(0,*) '>>> ERROR (inter2) <<<  Data do not span '//
     >               'grid.  '
          WRITE(0,*) '                        Use ADDPNT to '//
     >               'expand data and re-run.'
          STOP
      ENDIF

*  find the integral of each grid interval and use this to 
*  calculate the average y value for the interval      
*  xgl and xgu are the lower and upper limits of the grid interval

      jstart = 1
      ngintv = ng - 1
      DO 50, i = 1,ngintv

* initalize:

            area = 0.0
            xgl = xg(i)
            xgu = xg(i+1)

*  discard data before the first grid interval and after the 
*  last grid interval
*  for internal grid intervals, start calculating area by interpolating
*  between the last point which lies in the previous interval and the
*  first point inside the current interval

            k = jstart

            IF (k .LE. n-1) THEN

*  if both points are before the first grid, go to the next point
   30         CONTINUE
                IF (x(k+1) .LE. xgl) THEN
                   jstart = k - 1
                   k = k+1
                   IF (k .LE. n-1) GO TO 30
                ENDIF

*  if the last point is beyond the end of the grid, complete and go to the next
*  grid
   40         CONTINUE
                 IF ((k .LE. n-1) .AND. (x(k) .LT. xgu)) THEN          

                    jstart = k-1

* compute x-coordinates of increment

                    a1 = MAX(x(k),xgl)
                    a2 = MIN(x(k+1),xgu)

*  if points coincide, contribution is zero

                    IF (x(k+1).EQ.x(k)) THEN
                       darea = 0.e0
                    ELSE
                       slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                       b1 = y(k) + slope*(a1 - x(k))
                       b2 = y(k) + slope*(a2 - x(k))
                       darea = (a2 - a1)*(b2 + b1)/2.
c                       print *,a2,a1,k,y(k),slope,b2,b1,darea
                    ENDIF


*  find the area under the trapezoid from a1 to a2

                    area = area + darea

* go to next point
              
                    k = k+1
                    GO TO 40

                ENDIF

            ENDIF

*  calculate the average y after summing the areas in the interval
            yg(i) = area/(xgu - xgl)
            

   50 CONTINUE
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE inter3(ng,xg,yg, n,x,y, FoldIn)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on a set of bins onto a different set of target     =*
*=  bins.                                                                    =*
*=  The input data are given on a set of bins (representing the integral     =*
*=  of the input quantity over the range of each bin) and are being matched  =*
*=  onto another set of bins (target grid).  A typical example would be an   =*
*=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
*=  vals, that has to be matched onto the working wavelength grid.           =*
*=  The resulting area in a given bin of the target grid is calculated by    =*
*=  simply adding all fractional areas of the input data that cover that     =*
*=  particular target bin.                                                   =*
*=  Some caution should be used near the endpoints of the grids.  If the     =*
*=  input data do not span the full range of the target grid, the area in    =*
*=  the "missing" bins will be assumed to be zero.  If the input data extend =*
*=  beyond the upper limit of the target grid, the user has the option to    =*
*=  integrate the "overhang" data and fold the remaining area back into the  =*
*=  last target bin.  Using this option is recommended when re-gridding      =*
*=  vertical profiles that directly affect the total optical depth of the    =*
*=  model atmosphere.                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
*=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
*=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
*=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
*=           y-value for bin i (i = 1..NG-1)                                 =*
*=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
*=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
*=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
*=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
*=           y-value for bin i (i = 1..N-1)                                  =*
*=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
*=           FoldIn = 0 -> No folding of "overhang" data                     =*
*=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
*=                         last target bin                                   =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  06/96  Added FoldIn switch                                               =*
C R. F. Esswein 020214 Change all REAL declarations to REAL*8
C                               Use generic names of intrinsic functions.
C       M. Claire       060806  integrated into Kevin's code
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      
* input:
      INTEGER n, ng
      REAL*8 xg(ng)
      REAL*8 x(n), y(n)

      INTEGER FoldIn

* output:
      REAL*8 yg(ng)

* local:
      REAL*8 a1, a2, sum
      REAL*8 tail
      INTEGER jstart, i, jl, k
*_______________________________________________________________________

* check whether flag given is legal
      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
         WRITE(0,*) '                        Must be 0 or 1'
         STOP
      ENDIF

* do interpolation

      jstart = 1

      DO 30, i = 1, ng - 1

         yg(i) = 0.
         sum = 0.
         jl = jstart

         IF (jl .LE. n-1) THEN

   20      CONTINUE

             IF (x(jl+1) .LT. xg(i)) THEN
                jstart = jl
                jl = jl+1
                IF (jl .LE. n-1) GO TO 20
             ENDIF               

   25      CONTINUE

             IF ((x(jl) .LE. xg(i+1)) .AND. (jl .LE. n-1)) THEN

                a1 = MAX(x(jl),xg(i))
                a2 = MIN(x(jl+1),xg(i+1))

                sum = sum + y(jl) * (a2-a1)/(x(jl+1)-x(jl))
                jl = jl+1
                GO TO 25

             ENDIF

           yg(i) = sum 

         ENDIF

   30 CONTINUE


* if wanted, integrate data "overhang" and fold back into last bin

      IF (FoldIn .EQ. 1) THEN

         jl = jl-1
         a1 = xg(ng)     ! upper limit of last interpolated bin
         a2 = x(jl+1)     ! upper limit of last input bin considered

*        do folding only if grids don't match up and there is more input 
         IF ((a2 .GT. a1) .OR. (jl+1 .LT. n)) THEN
           tail = y(jl) * (a2-a1)/(x(jl+1)-x(jl))
           DO k = jl+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF

      ENDIF
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE inter4(ng,xg,yg, n,x,y, FoldIn)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on a set of bins onto a different set of target     =*
*=  bins.                                                                    =*
*=  The input data are given on a set of bins (representing the integral     =*
*=  of the input quantity over the range of each bin) and are being matched  =*
*=  onto another set of bins (target grid).  A typical example would be an   =*
*=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
*=  vals, that has to be matched onto the working wavelength grid.           =*
*=  The resulting area in a given bin of the target grid is calculated by    =*
*=  simply adding all fractional areas of the input data that cover that     =*
*=  particular target bin.                                                   =*
*=  Some caution should be used near the endpoints of the grids.  If the     =*
*=  input data do not span the full range of the target grid, the area in    =*
*=  the "missing" bins will be assumed to be zero.  If the input data extend =*
*=  beyond the upper limit of the target grid, the user has the option to    =*
*=  integrate the "overhang" data and fold the remaining area back into the  =*
*=  last target bin.  Using this option is recommended when re-gridding      =*
*=  vertical profiles that directly affect the total optical depth of the    =*
*=  model atmosphere.                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
*=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
*=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
*=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
*=           y-value for bin i (i = 1..NG-1)                                 =*
*=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
*=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
*=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
*=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
*=           y-value for bin i (i = 1..N-1)                                  =*
*=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
*=           FoldIn = 0 -> No folding of "overhang" data                     =*
*=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
*=                         last target bin                                   =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  06/96  Added FoldIn switch                                               =*
C R. F. Esswein 020214 Change all REAL declarations to REAL*8
C                               Use generic names of intrinsic functions.
C $Id$
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      
* input:
      INTEGER n, ng
      REAL*8 xg(ng)
      REAL*8 x(n), y(n)

      INTEGER FoldIn

* output:
      REAL*8 yg(ng)

* local:
      REAL*8 a1, a2, sum
      REAL*8 tail
      INTEGER jstart, i, jl, k
*_______________________________________________________________________

* check whether flag given is legal
      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
         WRITE(0,*) '                        Must be 0 or 1'
         STOP
      ENDIF

* do interpolation

      jstart = 1

      DO 30, i = 1, ng - 1

         yg(i) = 0.
         sum = 0.
         jl = jstart

         IF (jl .LE. n-1) THEN

   20      CONTINUE

             IF (x(jl+1) .LT. xg(i)) THEN
                jstart = jl
                jl = jl+1
                IF (jl .LE. n-1) GO TO 20
             ENDIF               

   25      CONTINUE

           IF ((x(jl) .LE. xg(i+1)) .AND. (jl .LE. n-1)) THEN

              a1 = MAX(x(jl),xg(i))
              a2 = MIN(x(jl+1),xg(i+1))

              sum = sum + y(jl) * (a2-a1)        !inter3.f divides by width here

              jl = jl+1
              GO TO 25

           ENDIF

           yg(i) = sum /(xg(i+1)-xg(i))         !inter3.f does not divide here

        ENDIF

 30   CONTINUE


* if wanted, integrate data "overhang" and fold back into last bin

      IF (FoldIn .EQ. 1) THEN

         jl = jl-1
         a1 = xg(ng)     ! upper limit of last interpolated bin
         a2 = x(jl+1)     ! upper limit of last input bin considered

*        do folding only if grids don't match up and there is more input 
         IF ((a2 .GT. a1) .OR. (jl+1 .LT. n)) THEN
           tail = y(jl) * (a2-a1)/(x(jl+1)-x(jl))
           DO k = jl+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF

      ENDIF
*_______________________________________________________________________

      RETURN
      END
