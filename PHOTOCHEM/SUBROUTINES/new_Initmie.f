      SUBROUTINE INITMIE(nw,wl,frak,monsize)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)

      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/MBLOK.inc'


!local variables for hydrocarbon optical properties
      DIMENSION WAVLS(108), WAVUS(108), QEXTSTAND(108,34)
      dimension W0STAND(108,34), GSTAND(108,34)
      dimension W0STANDNEW(108,34)

!C-GA dimensions of fractal hydrocarbon arrays (they have more points)
!C-GA because the fractals were recalculated for a wide wavelength range
      DIMENSION WAVLSF(118), WAVUSF(118), QEXTSTANDF(118,34)
      dimension W0STANDF(118,34), GSTANDF(118,34)
      dimension W0STANDNEWF(118,34)

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
      ierr = 0      

C-AP RADIUS of particles for which Mie calculations were run
C-AP **********************************************************
      DATA RSTAND/0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
     2  0.008, 0.009, 0.01, 0.03, 0.05, 0.07, 0.1, 0.13, 0.15,
     3  0.17, 0.2, 0.23, 0.27, 0.3, 0.33, 0.37, 0.4, 0.43, 0.47,
     4  0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1., 2./
C-AP ***********************************************************

      print *, 'monsize is: ', monsize
C************SPHERICAL***************************************** 
      if (frak.eq.0) then  !using spherical sized arrays (MIE)
         DO k=1,34
            RSTAND(k) = RSTAND(k)/10000.
         ENDDO

            root='PHOTOCHEM/DATA/MIE/fitmythol'
            print *,'using spherical MIE data'
            filenames=[
     $           '0001.DAT','0002.DAT','0003.DAT','0004.DAT','0005.DAT',
     $           '0006.DAT','0007.DAT','0008.DAT','0009.DAT','001.DAT ',
     $           '003.DAT ','005.DAT ','007.DAT ','01.DAT  ','013.DAT ',
     $           '015.DAT ','017.DAT ','02.DAT  ','023.DAT ','027.DAT ',
     $           '03.DAT  ','033.DAT ','037.DAT ','04.DAT  ','043.DAT ',
     $           '047.DAT ','05.DAT  ','055.DAT ','06.DAT  ','07.DAT  ',
     $           '08.DAT  ','09.DAT  ','1.DAT   ','2.DAT   ']
 
        

         do j=1,34              !there are 34 particle sizes for hydrocarbons
            open(unit=125,file=trim(root)//trim(filenames(j))) !open up MIE FILES

            do i=1,108          !there are 108 wavlength bins in the spherical data files
             READ(125,*) WAVLS(I),WAVUS(I),W0STAND(I,j),QEXTSTAND(I,j),
     2              GSTAND(I,j)
            enddo
          close (125)
         enddo   

 998  FORMAT(1x,F6.1,1x,F6.1,2x,F6.4,2x,F6.4,2x,F6.4)

      do j=1,34 !for each particle size
         n1=108 !number of wavelength points
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
         extend=1

         if (extend.eq.1) then 
            CALL addpnt(x1,y1,kw,n1,wl(1),W0STAND(1,j))  
            CALL addpnt(x1,y1,kw,n1,wl(1)*(1.-deltax),zero)  
         else 
            CALL addpnt(x1,y1,kw,n1,x1(1)*(1.-deltax),zero)  
         endif
         CALL addpnt(x1,y1,kw,n1,                  zero,zero)
         CALL addpnt(x1,y1,kw,n1,x1(n1)*(1.+deltax),zero)
         CALL addpnt(x1,y1,kw,n1,                   biggest,zero)

         CALL inter2(nw+1,wl,yg1,n1,x1,y1,0) !inter2 is discrete grid points to bins

         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
            STOP
         ENDIF

         if (extend.eq.1) then      
!     extending hydrocarbon optical properties to the UV
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
!     extending hydrocarbon optical properties to the UV
            CALL addpnt(x3,y3,kw,n3,wl(1),GSTAND(1,j))  
            CALL addpnt(x3,y3,kw,n3,wl(1)*(1.-deltax),zero)  
         else
            CALL addpnt(x3,y3,kw,n3,x3(1)*(1.-deltax),zero)  
         endif
         
         CALL addpnt(x3,y3,kw,n3,x3(1)*(1.-deltax),zero)  
         CALL addpnt(x3,y3,kw,n3,                  zero,zero)
         CALL addpnt(x3,y3,kw,n3,x3(n3)*(1.+deltax),zero)
         CALL addpnt(x3,y3,kw,n3,                   biggest,zero)
         CALL inter2(nw+1,wl,yg3,n3,x3,y3,0) !inter2 is discrete grid points to bins
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
            STOP
         ENDIF

         do k=1,nw
            W0HC(k,j)=yg1(k)    !single scattering albedo for hcaer 
            QEXTHC(k,j)=yg2(k)  !Extinction coefficient for hcaer
            GHC(k,j)=yg3(k)     !asymmetry param for hcaer
            
         enddo
         

      enddo                     !end for each particle size

      endif  
***************END SPHERES*****************************
**************FRACTALS BELOW HERE**********************

      if (frak.eq.1) then       !use fractal  MIE data

         npts=118

         DO k=1,34
            RSTAND(k) = RSTAND(k)/10000.
         ENDDO

         print *, 'using fractal MIE data'

C-GA fractal files read in from different folders compared to
C-   spherical particles
         if (monsize.eq.0.) then 
            root = 'PHOTOCHEM/DATA/MIE/f0/fractopts'
            print *,'using monsize 0: 0.05 um particles'
       
         else if (monsize.eq.1.) then 
            root = 'PHOTOCHEM/DATA/MIE/f1/fractopts'
            print *,'using monsize 1: 0.01 um particles'
            
         else if (monsize.eq.2.) then 
            root = 'PHOTOCHEM/DATA/MIE/f2/fractopts'
            print *,'using monsize 2: 0.02 um particles'
            
         else if (monsize.eq.3.) then 
            root = 'PHOTOCHEM/DATA/MIE/f3/fractopts'
            print *,'using monsize 3: 0.07 um particles'
            
         else if (monsize.eq.4.) then 
            root = 'PHOTOCHEM/DATA/MIE/f4/fractopts'
            print *,'using monsize 4: 0.10 um particles'
         end if
         
        
      filenames=['0.001um.txt','0.002um.txt','0.003um.txt','0.004um.txt'
     $          ,'0.005um.txt','0.006um.txt','0.007um.txt','0.008um.txt'
     $          ,'0.009um.txt','0.010um.txt','0.030um.txt','0.050um.txt'
     $          ,'0.070um.txt','0.100um.txt','0.130um.txt','0.150um.txt'
     $          ,'0.170um.txt','0.200um.txt','0.230um.txt','0.270um.txt'
     $          ,'0.300um.txt','0.330um.txt','0.370um.txt','0.400um.txt'
     $          ,'0.430um.txt','0.470um.txt','0.500um.txt','0.550um.txt'
     $          ,'0.600um.txt','0.700um.txt','0.800um.txt','0.900um.txt'
     $          ,'1.000um.txt','2.000um.txt']


      


      do j=1,34                 !there are 34 particle sizes for hydrocarbons
         open(unit=125,file=trim(root)//trim(filenames(j))) !open up MIE FILES

       
         do i=1,npts            !there are 118 wavlength bins in the FRACTAL data files 108->118
          READ(125,*) WAVLSF(I),WAVUSF(I),W0STANDF(I,j),QEXTSTANDF(I,j),
     2        GSTANDF(I,j)

      enddo
      
      close (125)
      enddo   

 999  FORMAT(1x,F6.1,1x,F6.1,2x,F6.4,2x,F6.4,2x,F6.4)

      do j=1,34
         n1=npts                !108->118
         n2=n1
         n3=n1
      do i=1,n1
         x1(i)=wavlsf(i)
         x2(i)=x1(i)
         x3(i)=x1(i)
         y1(i)=W0STANDF(i,j)    !single scattering albedo
         y2(i)=QEXTSTANDF(i,j)  !extinction efficiency Qext
         y3(i)=GSTANDF(i,j)     !asymmetry parameter
      enddo   

C-GA we don't need to to the 'extend' thing for the fractals...
         extend=0
   
     
      CALL addpnt(x1,y1,kw,n1,x1(1)*(1.-deltax),zero)  
      CALL addpnt(x1,y1,kw,n1,                  zero,zero)
      CALL addpnt(x1,y1,kw,n1,x1(n1)*(1.+deltax),zero)
      CALL addpnt(x1,y1,kw,n1,                   biggest,zero)

      CALL inter2(nw+1,wl,yg1,n1,x1,y1,0) !inter2 is discrete grid points to bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
         STOP
      ENDIF

  
      CALL addpnt(x2,y2,kw,n2,x2(1)*(1.-deltax),zero)    
      CALL addpnt(x2,y2,kw,n2,                  zero,zero)
      CALL addpnt(x2,y2,kw,n2,x2(n2)*(1.+deltax),zero)
      CALL addpnt(x2,y2,kw,n2,                   biggest,zero)
      CALL inter2(nw+1,wl,yg2,n2,x2,y2,0) !inter2 is discrete grid points to bins

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
         STOP
      ENDIF

    
      
      CALL addpnt(x3,y3,kw,n3,x3(1)*(1.-deltax),zero)  
      CALL addpnt(x3,y3,kw,n3,x3(1)*(1.-deltax),zero)  
      CALL addpnt(x3,y3,kw,n3,                  zero,zero)
      CALL addpnt(x3,y3,kw,n3,x3(n3)*(1.+deltax),zero)
      CALL addpnt(x3,y3,kw,n3,                   biggest,zero)
      CALL inter2(nw+1,wl,yg3,n3,x3,y3,0)   !inter2 is discrete grid points to bins
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, ' ***Something wrong in Initmie***'
         STOP
      ENDIF

      do k=1,nw
         W0HC(k,j)=yg1(k)
         QEXTHC(k,j)=yg2(k)
         GHC(k,j)=yg3(k)

      enddo

      enddo
      
      endif


      END



