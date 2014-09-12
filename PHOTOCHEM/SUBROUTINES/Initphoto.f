      SUBROUTINE INITPHOTO(sq,columndepth,zy,nw,timega,IO2,msun)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      Integer nw,j,IO2
      real*8 mass
      REAL*8 dum
      REAL*8 f(kw)   !temp flux vector
      real*8 SQ(kj,nz,kw),o2col(NZ),columndepth(KJ,NZ) 
      character*8 PLANET,ISPEC
      character*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/QBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'

      CALL gridw(nw,wavl,wav,wavu,LGRID)    !this function defines the wavelength grid and returns
c nw, the number of points on the grid
c wavl, a vector of lower grid points  !note wavl(nw+1)=wavu(nw) i.e. final grid point for interpolation
c wavu, a vector of upper grid points (i.w. wavl+delta)
c wav, a vector of centered values (i.e. (wavl + wavu)/2 )


      CALL readflux(nw,wavl,flux,timega,msun)
c this subroutine returns the flux data interpolated to the wavelength grid along with 
c Claire et al. 2012 corrections for the specified time (timeGa set in PLANET.dat)

      print *, 'for a solar flux from ',timega,' Ga'

C *****  ***** ***** READ THE PHOTOLYSIS DATAFILE ***** ***** *****
C *****  ***** ***** READ THE PHOTOLYSIS DATAFILE ***** ***** *****
c eventually, this will die...

      READ(3,100)
 100  FORMAT(/)
      DO 41 L=1,108 
  41  READ(3,*)
 101  FORMAT(8X,F6.1,1X,F6.1,2X,E9.2,4(3X,E8.1))

      READ(3,102)
 102  FORMAT(////)
      DO 42 L=1,17
  42  READ(3,*)
 103  FORMAT(4X,I1,1X,9(1X,E13.6))
  
      READ(3,102)
      DO 43 L=1,17
  43  READ(3,*) 
 104  FORMAT(4X,I1,1X,5(1X,E13.6))
 
      READ(3,102)        !SO2HZ is still used. what is the best way to deal?
      DO 44 L=1,35
  44  READ(3,105) SO2HZ(L),dum,dum,dum,dum,dum,
     2  dum,dum,dum,dum,dum
 105  FORMAT(5X,11(E8.1,1X))
 
      READ(3,106)
 106  FORMAT(//)
      DO 45 L=1,68
  45  READ(3,*)
 107  FORMAT(6X,8E8.1,2(1X,F4.2))
 
      READ(3,106)
      DO 46 L=1,68
 46    READ(3,*)  

 108  FORMAT(6X,9E8.1)
 
!these lines might need to be kept...
!what I don't know is where the coefficients come from or how they would scale with a varying wavelength grid...

      READ(3,109)
 109  FORMAT(//////)
      DO 58 L=1,17
      READ(3,110) NK(L),(ALPHAP(L,K),K=1,4)
 110  FORMAT(6X,I1,1X,4F13.5)
  58  READ(3,111) (BETA(L,K),K=1,4)
 111  FORMAT(8X,4E13.4/)
 
      READ(3,106)
      DO 67 L=1,68
  67  READ(3,*)
 113  FORMAT(15X,7(E9.3,4X))
 
C *****  ***** ***** END READING PHOTOLYSIS DATAFILE ***** ***** *****
C *****  ***** ***** END READING PHOTOLYSIS DATAFILE ***** ***** *****

c at some point do away with all of the above, just taking what's relevant into single datafiles



C read in relevant cross sections and set the Jnumbers,photospec, and photolabel

      do i=1,kj
         photolabel(i)=''
      enddo   


      j=1
      do i=1,ks   !number of photolysis species
c         print *, photospec(i),ISPEC(INT(photospec(i))),j
       CALL XS(ISPEC(INT(photospec(i))),nw,wavl,wav,T,DEN,j,sq,
     $         columndepth,zy,IO2)         
       enddo   

c       J=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'CH4    ')

c       J=JS8L
c       print *, (wav(L),sq(J,20,L),L=1,NW)
c       print *,(wav(L),sq(J,20,L),sq(J+1,10,L),L=1,NW)
c       print *, (wav(L),sq(J,20,L),sq(J+1,20,L),sq(J+2,20,L)
c     $          ,L=1,NW)

c       print *, 'stopping in Initphoto'
c       stop

      RETURN
      END
