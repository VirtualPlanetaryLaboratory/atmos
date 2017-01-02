
        SUBROUTINE AERABSDATA(frak, monsize)        
        DIMENSION wavirst(55), wavsolst(38)
        CHARACTER :: DIRINOUT*8,DIRDATA*10
        COMMON/DIR/DIRINOUT,DIRDATA
        COMMON/HYDROCARB/Qextirst(73,55),w0irst(73,55),
     &  girst(73,55),Qextsolst(73,38),w0solst(73,38),gsolst(73,38),
     &  radstand(73)
C  radius of particles is given in cm
        DATA radstand/0.001E-5, 0.002E-5, 0.003E-5, 0.004E-5,
     &  0.005E-5, 0.006E-5, 0.007E-5, 0.008E-5, 0.009E-5, 
     &  0.001E-4, 0.002E-4, 0.003E-4, 0.004E-4,
     &  0.005E-4, 0.006E-4, 0.007E-4, 0.008E-4, 0.009E-4, 
     &  0.01E-4, 0.02E-4, 0.03E-4, 0.04E-4, 0.05E-4,
     &  0.06E-4, 0.07E-4,
     &  0.08E-4, 0.09E-4, 0.1E-4, 0.2E-4, 0.3E-4, 0.4E-4, 0.5E-4,
     &  0.6E-4, 0.7E-4, 0.8E-4, 0.9E-4,
     &  1.E-4, 2.E-4, 3.E-4, 4.E-4, 5.E-4, 6.E-4, 7.E-4, 8.E-4,
     &  9.E-4, 1.E-3/

        print *, 'monsize is: ', monsize
        IF (frak.eq.1) THEN
         if (monsize.eq.0) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract.DAT')
         print *, 'doing fractal particles monsize 0.05um'
         endif

         if (monsize.eq.1) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract1.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract1.DAT')
         print *, 'doing fractal particles monsize 0.01um'
         endif

         if (monsize.eq.2) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract2.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract2.DAT')
         print *, 'doing fractal particles monsize 0.02um'
         endif

         if (monsize.eq.3) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract3.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract3.DAT')
         print *, 'doing fractal particles monsize 0.07um'
         endif

         if (monsize.eq.4) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract4.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract4.DAT')
         print *, 'doing fractal particles monsize 0.1um'
         endif

        ELSE
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotal.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotal.DAT')
         print *,'doing spherical particles'
        ENDIF
C READING IR particle absorption data
         DO J=1,73
          READ(40,100)
         DO i=1,55
         READ(40,*) wavirst(i),Qextirst(J,i),w0irst(J,i),
     &   girst(J,i)
         ENDDO
         ENDDO
C READING SOLAR particle absorption data
         DO J=1,73
          READ(41,100)
         DO i=1,38
         READ(41,*) wavsolst(i),Qextsolst(J,i),w0solst(J,i),
     &   gsolst(J,i)
         ENDDO
         ENDDO
  100    format(/)
         CLOSE (40)
         CLOSE (41)
         RETURN
         END
