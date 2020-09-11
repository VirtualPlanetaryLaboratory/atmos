      SUBROUTINE INITPHOTO(sq,columndepth,zy,nw,timega,IO2,pstar,uvscale
     $       )
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      Integer nw,j,IO2
      real*8 mass
      REAL*8 dum
      character*8 pstar
      real*8 SQ(kj,nz,kw),columndepth(KJ,NZ)
      character*8 PLANET,ISPEC
      character*11 photolabel
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/QBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'



c this function defines the wavelength grid and returns
      CALL gridw(nw,wavl,wav,wavu,LGRID)
c nw, the number of points on the grid
c wavl, a vector of lower grid points  !note wavl(nw+1)=wavu(nw) i.e. final grid point for interpolation
c wavu, a vector of upper grid points (i.w. wavl+delta)
c wav, a vector of centered values (i.e. (wavl + wavu)/2 )
c this subroutine returns the flux data interpolated to the wavelength grid along with
c Claire et al. 2012 corrections for the specified time (timeGa set in PLANET.dat)
      CALL readflux(nw,wavl,timega,pstar,uvscale)

      print *, 'for a solar flux from ',timega,' Ga'

C *****  ***** ***** READ THE PHOTOLYSIS DATAFILE ***** ***** *****
c eventually, this will die...

      READ(3,100)
 100  FORMAT(/)
      DO L=1,108
      READ(3,*)
      END Do

      READ(3,102)
 102  FORMAT(////)
      DO L=1,17
      READ(3,*)
      END DO

      READ(3,102)
      DO L=1,17
      READ(3,*)
      END DO

c Below SO2HZ is still used. what is the best way to deal?
      READ(3,102)
      DO L=1,35
      READ(3,105) SO2HZ(L),dum,dum,dum,dum,dum,
     2  dum,dum,dum,dum,dum
      END DO
 105  FORMAT(5X,11(E8.1,1X))

      READ(3,106)
 106  FORMAT(//)
      DO L=1,68
      READ(3,*)
      END DO

      READ(3,106)
      DO L=1,68
       READ(3,*)
      END DO


!these lines might need to be kept...
!what I don't know is where the coefficients come from or how they would scale with a varying wavelength grid...

      READ(3,109)
 109  FORMAT(//////)
      DO L=1,17
      READ(3,110) NK(L),(ALPHAP(L,K),K=1,4)
 110  FORMAT(6X,I1,1X,4F13.5)
      READ(3,111) (BETA(L,K),K=1,4)
      END DO
 111  FORMAT(8X,4E13.4/)

      READ(3,106)
      DO L=1,68
      READ(3,*)
      END DO

C *****  ***** ***** END READING PHOTOLYSIS DATAFILE ***** ***** *****
c at some point do away with all of the above, just taking what's relevant into single datafiles


C read in relevant cross sections and set the Jnumbers,photospec,and photolabel

      do i=1,kj
         photolabel(i)=''
      enddo


      j=1
      do i=1,ks   !number of photolysis species
c         print *, photospec(i),ISPEC(INT(photospec(i))),j
       CALL newXS(ISPEC(INT(photospec(i))),nw,wavl,wav,T,DEN,j,sq,
     $         columndepth,zy,IO2,lgrid)
       enddo
c       stop
      RETURN
      END
