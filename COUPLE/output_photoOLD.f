C------------------------------------------------------------------------------
C-KK        Code added Sept 27 2001 for integration with the photochemical code
C-KK         atm_chem.f   
C------------------------------------------------------------------------------
  
      SUBROUTINE OUTPUT_PHOTO(T, FI, water, ALT,DELZkm,JCOLD)
      INCLUDE '../CLIMA/INCLUDE/header.inc'
      INCLUDE '../PHOTOCHEM/INPUTFILES/parameters.inc'
     
      DIMENSION T(ND), FI(4,ND), alt_new(NZ), T_new(NZ)
      DIMENSION water(NZ), ALT(ND)
         IFILL=0
        alt_new(NZ) = DELZkm/2.
        do i = (NZ-1), 1, -1
C          interpolate temperatures at correct altitudes
          alt_new(i) = alt_new(i+1) + DELZkm
        end do
        print*,ALT(1),alt_new(1)
        IF(ALT(1).LT.alt_new(1))IFILL =1
        ISTART = ND
      
        DO J = NZ, 1, -1
         DO i = ISTART,1,-1
            IF(IFILL.eq.1.and. i.eq.1) then
             NSAVE=J
             GO TO 320
            endif
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
        T1log=LOG(T(IS))
        T2log=LOG(T(IS1))
        T_temp = FR*T2log + (1-FR)*T1log
        T_new(J) = EXP(T_temp)
C-KK                        begin water interpolation
        FH1log = LOG(FI(1,IS))
        FH2log = LOG(FI(1,IS1))
        water_temp = FR*FH2log + (1-FR)*FH1log
        water(J) = EXP(water_temp)
        END DO
        GO TO 310

 320    do i=NSAVE,1,-1
         water(i)=FI(1,1)
         T_new(i)=T(1)
       enddo



C-KK        This statement reverses the indexing (has to be done here and not
C-KK        prior because earlier reversal would confuse the issue during the
C-KK        interpolation procedure. The interpolated altitudes, temperatures,
C-KK        water profile and ozone profile are recorded to file fromClima2Photo.dat

 310   JCOLD = 0

        DO J = NZ, 2, -1
           WRITE(54,*) alt_new(J), T_new(J), water(J)
            IF (JCOLD .EQ. 0) THEN
                IF (T_new(J) .LT. T_new(J-1)) JCOLD = NZ - J
c               IF (water(J) .LT. water(J-1)) JCOLD = NZ - J
           END IF
        END DO

        J=1
        WRITE(54,*) alt_new(J), T_new(J), water(J)
        CLOSE(54)

        CONTINUE
        END
