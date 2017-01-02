      SUBROUTINE CHOOSE_STAR(STARR,SOLINT)

c This subroutine reads the spectra from stars other than the SUN

      PARAMETER(NSOL=38)
      DIMENSION SOLINT(NSOL)
      CHARACTER :: STARR*3,DIRINOUT*2,DIRDATA*4
      COMMON/DIR/DIRINOUT,DIRDATA

      IF(STARR=="MdV") THEN
c   File with M star data provided by Martin cohen
       OPEN(UNIT=81,FILE= DIRINOUT//'/Model_T5100_surf.dat')
       READ(81,*)
       DO j=1,NSOL
         READ(981,*)xl,xf,x,x
         SOLINT(j) = xf
       ENDDO
       close(81)
c       IF(STARR=="GJ581") THEN
c        read(10,*) ! skips one line before reading 
c       DO I =1,38
c	  read(10,100)SOLINT(I)
c	  print 100, SOLINT(I)
c        ENDDO
      
c 100  format(6x,1pe11.4)
      ELSE
C This file contains normalized fluxes for K, G, and F type stars
c   NOTICE: G star is NOT the SUN. Data from Martin Cohen
       OPEN(UNIT=80,FILE= DIRDATA//'/fluxesKGF_surf.pdat')  
       READ(80,*)
        DO I=1,NSOL
         READ(80,*) ll,xf1,xf2,xf3,x,x
         IF(STARR=="K2V") SOLINT(I) = xf1
         IF(STARR=="G2V") SOLINT(I) = xf2    	
         IF(STARR=="F2V") SOLINT(I) = xf3
        ENDDO
       close(80)
      ENDIF

      RETURN
      END
 
