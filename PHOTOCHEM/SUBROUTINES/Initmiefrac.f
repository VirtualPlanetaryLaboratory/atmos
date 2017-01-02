      SUBROUTINE INITMIEFRAC(nw,wl,frak,monsize)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)

      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/MFBLOK.inc'


!local variables for hydrocarbon optical properties
      
 
!mie only has 108 lines, so putting in 118 breaks this for mie
!(but it works for fractals). 
     
      DIMENSION WAVLS(118), WAVUS(118), QEXTSTAND(118,34) !108->118
      dimension W0STAND(118,34), GSTAND(118,34)

      dimension W0STANDNEW(118,34)
      


      CHARACTER*31 root,filenames
      dimension filenames(34)


      REAL*8 deltax,biggest,zero
      PARAMETER (deltax = 1.E-4,biggest=1.E+36, zero=0.0)
* input
      INTEGER nw
      REAL*8 wl(nw+1)

* data arrays
      INTEGER n1, n2,n3 
      REAL*8 x1(kw), x2(kw), x3(kw)
      REAL*8 y1(kw), y2(kw), y3(kw)

* local
      REAL*8 yg(nw),yg1(nw),yg2(kw),yg3(kw)
      INTEGER i, iw, n
      INTEGER ierr
      REAL*8 monsize
      ierr = 0      

C-AP RADIUS of particles for which Mie calculations were run
C-AP **********************************************************
      DATA RSTAND/0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
     2  0.008, 0.009, 0.01, 0.03, 0.05, 0.07, 0.1, 0.13, 0.15,
     3  0.17, 0.2, 0.23, 0.27, 0.3, 0.33, 0.37, 0.4, 0.43, 0.47,
     4  0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1., 2./
C-AP ***********************************************************
!so 34 mie calcualtions were run.  I wonder how we will do this for sulfate which are 19 calculations from 0.01 to 1 micron


      !this was in the time-stepping loop of the original code, but that seems like a mistake
C-AP Since all model is in cm we should convert RSTAND
      if (frak.eq.1) then
         npts = 118
      else
         npts = 108
      endif



      DO k=1,34
       RSTAND(k) = RSTAND(k)/10000.
      ENDDO

      if (frak.eq.0) then  !use spherical MIE data
      root='PHOTOCHEM/DATA/MIE/fitmythol'
      filenames=['0001.DAT','0002.DAT','0003.DAT','0004.DAT','0005.DAT',
     $       '0006.DAT','0007.DAT','0008.DAT','0009.DAT','001.DAT ',
     $       '003.DAT ','005.DAT ','007.DAT ','01.DAT  ','013.DAT ',
     $       '015.DAT ','017.DAT ','02.DAT  ','023.DAT ','027.DAT ',
     $       '03.DAT  ','033.DAT ','037.DAT ','04.DAT  ','043.DAT ',
     $       '047.DAT ','05.DAT  ','055.DAT ','06.DAT  ','07.DAT  ',
     $       '08.DAT  ','09.DAT  ','1.DAT   ','2.DAT   ']
      else  !use fractal MIE data
         print *, 'using fractal MIE data'

       if (monsize.eq.0) then 
          root = 'PHOTOCHEM/DATA/MIE/f0/fractopts'
       endif
       if (monsize.eq.1) then 
          root = 'PHOTOCHEM/DATA/MIE/f1/fractopts'
       endif
       if (monsize.eq.2) then 
          root = 'PHOTOCHEM/DATA/MIE/f2/fractopts'
       endif
       if (monsize.eq.3) then 
          root = 'PHOTOCHEM/DATA/MIE/f3/fractopts'
       endif
       if (monsize.eq.4) then 
          root = 'PHOTOCHEM/DATA/MIE/f4/fractopts'
       endif

      filenames=['0.001um.txt','0.002um.txt','0.003um.txt','0.004um.txt'
     $          ,'0.005um.txt','0.006um.txt','0.007um.txt','0.008um.txt'
     $          ,'0.009um.txt','0.010um.txt','0.030um.txt','0.050um.txt'
     $          ,'0.070um.txt','0.100um.txt','0.130um.txt','0.150um.txt'
     $          ,'0.170um.txt','0.200um.txt','0.230um.txt','0.270um.txt'
     $          ,'0.300um.txt','0.330um.txt','0.370um.txt','0.400um.txt'
     $          ,'0.430um.txt','0.470um.txt','0.500um.txt','0.550um.txt'
     $          ,'0.600um.txt','0.700um.txt','0.800um.txt','0.900um.txt'
     $          ,'1.000um.txt','2.000um.txt']


      endif


      do j=1,34  !there are 34 particle sizes for hydrocarbons
       open(unit=125,file=trim(root)//trim(filenames(j)))  !open up MIE FILES

       
        do i=1,npts  !there are 118 wavlength bins in the FRACTAL data files 108->118
         READ(125,*) WAVLS(I),WAVUS(I),W0STAND(I,j),QEXTSTAND(I,j),
     2   GSTAND(I,j)
c         print *,J,I, GSTAND(I,J)
        enddo
      
      close (125)
      enddo   

c      stop

 998  FORMAT(1x,F6.1,1x,F6.1,2x,F6.4,2x,F6.4,2x,F6.4)

      do j=1,34
      n1=npts !108->118
      n2=n1
      n3=n1
      do i=1,n1
       x1(i)=wavls(i)
       x2(i)=x1(i)
       x3(i)=x1(i)
       y1(i)=W0STAND(i,j)
       y2(i)=QEXTSTAND(i,j)
       y3(i)=GSTAND(i,j)
      enddo   

!next, we interpolate these values to our wavelength grid, with an option to extend the optical properties to the UV

!assuming optical properties of the shortest wl bin of the data (~1750A) extend through Ly alpha
!per Jim Kasting
      if (frak.eq.0) then 
         extend=1               !change back to extend=1 if using MIE (not fractal) particles
      else
         extend=0
      endif

      if (extend.eq.1) then 
       CALL addpnt(x1,y1,kw,n1,wl(1),W0STAND(1,j))  
       CALL addpnt(x1,y1,kw,n1,wl(1)*(1.-deltax),zero)  
      else 
       CALL addpnt(x1,y1,kw,n1,x1(1)*(1.-deltax),zero)  
      endif
      CALL addpnt(x1,y1,kw,n1,                  zero,zero)
      CALL addpnt(x1,y1,kw,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kw,n1,                   biggest,zero)

      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0)   !inter2 is discrete grid points to bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
         STOP
      ENDIF

      if (extend.eq.1) then      
!extending hydrocarbon optical properties to the UV
      CALL addpnt(x2,y2,kw,n2,wl(1),QEXTSTAND(1,j))  
      CALL addpnt(x2,y2,kw,n2,wl(1)*(1.-deltax),zero)  
      else
      CALL addpnt(x2,y2,kw,n2,x2(1)*(1.-deltax),zero)  
      endif

      CALL addpnt(x2,y2,kw,n2,                  zero,zero)
      CALL addpnt(x2,y2,kw,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kw,n2,                   biggest,zero)
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,0) !inter2 is discrete grid points to bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
         STOP
      ENDIF

      if (extend.eq.1) then 
!extending hydrocarbon optical properties to the UV
      CALL addpnt(x3,y3,kw,n3,wl(1),GSTAND(1,j))  
      CALL addpnt(x3,y3,kw,n3,wl(1)*(1.-deltax),zero)  
      else
      CALL addpnt(x3,y3,kw,n3,x3(1)*(1.-deltax),zero)  
      endif

      CALL addpnt(x3,y3,kw,n3,x3(1)*(1.-deltax),zero)  
      CALL addpnt(x3,y3,kw,n3,                  zero,zero)
      CALL addpnt(x3,y3,kw,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kw,n3,                   biggest,zero)
      CALL inter2(nw+1,wl,yg3,n3,x3,y3,0)   !inter2 is discrete grid points to bins
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
         STOP
      ENDIF

c      print *, j, yg3

      do k=1,nw
      W0HC(k,j)=yg1(k)
      QEXTHC(k,j)=yg2(k)
      GHC(k,j)=yg3(k)

c       if (wl(k).eq. 2100) print *, wl(k),yg1(k),yg2(k),yg3(k)
      enddo


      enddo

      END



