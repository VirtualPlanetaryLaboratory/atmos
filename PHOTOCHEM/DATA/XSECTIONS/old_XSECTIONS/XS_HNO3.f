      SUBROUTINE XS_HNO3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
*=        HNO3 + hv -> OH + NO2                                              =*
*=  Cross section: Burkholder et al., 1993                                   =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
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

      IMPLICIT NONE
      INCLUDE 'photol_params.h'

* input

      INTEGER nw
      REAL*8 wl(kw), wc(kw)

      INTEGER nz

      REAL*8 tlev(kz)
      REAL*8 airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL*8 sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER n1, n2
      REAL*8 x1(kdata), x2(kdata)
      REAL*8 y1(kdata), y2(kdata)

* local

      REAL*8 yg1(kw), yg2(kw)
      INTEGER i, iw
      INTEGER ierr

**************** HNO3 photodissociation

       j = j + 1
       jlabel(j) = 'HNO3 -> OH + NO2'

* HNO3 cross section parameters from Burkholder et al. 1993

      OPEN(UNIT=kin,
     &  file=PHRT(1:INDEX(PHRT,' ')-1)//'DATAJ1/ABS/HNO3_burk.abs',
     &  STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      END DO
      n1 =  83
      n2 = n1
      DO i = 1, n1
         READ(kin,*) y1(i), y2(i)
         x1(i) = 184. + i*2.
         x2(i) = x1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,               zero,zero)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kdata,n1,            largest,zero)
      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
      CALL addpnt(x2,y2,kdata,n2,               zero,y2(1))
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
      CALL addpnt(x2,y2,kdata,n2,            largest,y2(n2))
      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield = 1
* correct for temperature dependence

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw) * 1.E-20
     $           * exp( yg2(iw)/1.e3*(tlev(i)-298.) )
         ENDDO
      ENDDO

      END
