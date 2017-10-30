
        SUBROUTINE AERABSDATA(frak, ihztype)        
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
     &  0.08E-4, 0.09E-4, 
     &  0.1E-4,  0.13E-4, 0.15E-4, 0.17E-4,
     &  0.2E-4,  0.23E-4, 0.25E-4, 0.27E-4,
     &  0.3E-4,  0.33E-4, 0.35E-4, 0.37E-4,
     &  0.4E-4,  0.43E-4, 0.45E-4, 0.47E-4,
     &  0.5E-4,  0.53E-4, 0.55E-4, 0.57E-4,      
     &  0.6E-4,  0.63E-4, 0.65E-4, 0.67E-4,
     &  0.7E-4,  0.73E-4, 0.75E-4, 0.77E-4,
     &  0.8E-4,  0.83E-4, 0.85E-4, 0.87E-4,
     &  0.9E-4, 0.93E-4, 0.95E-4, 0.97E-4,
     &  1.E-4, 2.E-4, 3.E-4, 4.E-4, 5.E-4, 6.E-4, 7.E-4, 8.E-4,
     &  9.E-4, 1.E-3/


        print *, 'ihztype is: ', ihztype
        IF (frak.eq.1) THEN
         if (ihztype.eq.0) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract.DAT')
         print *, 'doing fractal particles ihztype 0.05um (Khare)'
         endif

         if (ihztype.eq.1) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract1.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract1.DAT')
         print *, 'doing fractal particles ihztype 0.01um (Khare)'
         endif

         if (ihztype.eq.2) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract2.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract2.DAT')
         print *, 'doing fractal particles ihztype 0.02um (Khare)'
         endif

         if (ihztype.eq.3) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract3.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract3.DAT')
         print *, 'doing fractal particles ihztype 0.07um (Khare)'
         endif

         if (ihztype.eq.4) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract4.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract4.DAT')
         print *, 'doing fractal particles monsize 0.1um (Khare)'
         endif

         if (ihztype.eq.5) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract_kmt.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract_kmt.DAT')
         print *, 'doing fractal 0.5um Khare+Mahjoub+Tran'
         endif

         if (ihztype.eq.6) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract_kh.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract_kh.DAT')
         print *, 'doing fractal particles Khare+Hasenkopf'
         endif

         if (ihztype.eq.7) then 
         OPEN(UNIT=40,FILE= DIRDATA//'/irtotalfract_gk.DAT')
         OPEN(UNIT=41,FILE= DIRDATA//'/soltotalfract_gk.DAT')
         print *, 'doing fractal particles Gavalin+Khare'
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
