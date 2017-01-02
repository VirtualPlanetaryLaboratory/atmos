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
C	R. F. Esswein	020214	Change all REAL declarations to REAL*8
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
C	R. F. Esswein	020214	Change all REAL declarations to REAL*8
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
