C------------------------------------------------------------------------------
C-KK        Code added Sept 27 2001 for integration with the photochemical code
C-KK         now known as atm_chem.f (old primo3_g2.news.f) 
C------------------------------------------------------------------------------ 
      SUBROUTINE OUTPUT_PHOTO(T, FI, water, ALT)

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NZ=200)
      DIMENSION T(ND), FI(5,ND), alt_new(NZ), T_new(NZ)
      DIMENSION water(NZ), ALT(ND)
      CHARACTER :: DIRINOUT*2,DIRDATA*4
      COMMON/DIR/DIRINOUT,DIRDATA
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FH2

        alt_new(NZ) = 0.0
        do i = (NZ-1), 1, -1
C          interpolate temperatures at correct altitudes
          alt_new(i) = alt_new(i+1) + 0.5 !gna - used to be +1 but that's not how we are setting up photo's grid!
        end do

        ISTART = ND
        T_new(NZ) = T(ND)
        water(NZ) = FI(1,ND)
        
        DO J = (NZ-1), 1, -1
         DO i = ISTART,1,-1
           IS = i
           IF (ALT(IS) .GT. alt_new(J)) GOTO 350
         END DO
  350  CONTINUE
        IS1 = IS+1
        ISTART = IS
        DZI = ALT(IS) - ALT(IS1)
        DZJ = ALT(IS) - alt_new(J)
        FR = DZJ/DZI
C-KK                        begin temperature interpolation
        T1log=ALOG(T(IS))
        T2log=ALOG(T(IS1))
        T_temp = FR*T2log + (1-FR)*T1log
        T_new(J) = EXP(T_temp)
C-KK                        begin water interpolation
        FH1log = ALOG(FI(1,IS))
        FH2log = ALOG(FI(1,IS1))
        water_temp = FR*FH2log + (1-FR)*FH1log
        water(J) = EXP(water_temp)

        !print *, "FI(1,IS)", FI(1,IS)
        !print *, "water(J)", water(J)
        END DO


C-KK        This statement reverses the indexing (has to be done here and not
C-KK        prior because earlier reversal would confuse the issue during the
C-KK        interpolation procedure. The interpolated altitudes, temperatures,
C-KK        water profile and ozone profile are recorded to file TempH2Oout.dat

        JCOLD = 0

        DO J = NZ, 2, -1
           WRITE(116,*) alt_new(J), T_new(J), water(J)
           IF (JCOLD .EQ. 0) THEN
              IF (T_new(J) .LT. T_new(J-1)) JCOLD = NZ - J
           END IF
        END DO
        WRITE(116,*) alt_new(1), T_new(1), water(1)
        close (116)
        FCH4 = FI(3,ND)
        FCO2 = FI(2,ND)

c        CLOSE(14)
c        OPEN(unit=14,file= DIRINOUT//'/mixing_ratios.dat')        

c        WRITE(14,102) FAR, FCH4, FC2H6, FCO2, FN2, FO2,
c     &         FH2, FNO2, JCOLD
c  102        FORMAT(1PE10.3,10x,'!Argon'/,1PE10.3,10x,'!Methane'/
c     &  ,1PE10.3,10x,'!Ethane'/,1PE10.3,10x,'!Carbon Dioxide'/,
c     &  1PE10.3,10x,'!Nitrogen'/,1PE10.3,10x,'!Oxygen'/,
c     &  1PE10.3,10x,'!Hydrogen'/,1PE10.3,10x,'!Nitrogen Dioxide'/, 
c     &  I3,17x,'!Tropopause layer')

        END
