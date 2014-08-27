     
      SUBROUTINE INPUT_INTERP(NREAD,temp_alt, water, O3, Jcold, T, FI
     |   , FNO2_PHOTO, FNO2_NEW, FSO2_PHOTO, FSO2_NEW, CH4, CO2)
C-KK This subroutine interpolates the (alt),water, ozone file that is READ in 
C-KK from atm_chem.

      INCLUDE '../CLIMA/INCLUDE/header.inc'
      DIMENSION O3(NREAD),water(NREAD),temp_alt(NREAD)
      DIMENSION CH4(NREAD), CO2(NREAD), T(ND),FI(4,ND)
      DIMENSION FNO2_PHOTO(NREAD),FNO2_NEW(ND)
      DIMENSION FSO2_PHOTO(NREAD),FSO2_NEW(ND)
C
      CHARACTER :: DIRINOUT*8,DIRDATA*10
      COMMON/DIR/DIRINOUT,DIRDATA
      COMMON/ALTBLOK/DALT(ND-1),ALT(ND)

C-KK   both indices starting at ground level
       istart = 1
       j = ND
       print *, FI(1,j)
       print *, water(istart)
       FI(1,j) = water(istart)
       FI(2,j) = CO2(istart)
       FI(3,j) = CH4(istart)
       FI(4,j) = O3(istart)

       FNO2_NEW(j) = FNO2_PHOTO(1)
       FSO2_NEW(j) = FSO2_PHOTO(1)


       istart = istart+1
       top = temp_alt(NZ)           ! Was NREAD


C
      DO j = ND-1, 1, -1
        DO i = istart, NZ           ! Was NREAD
         IS = i
         if (temp_alt(i) .GT. ALT(j)) GOTO 350
        END DO
C
350   CONTINUE



C
C-TF  IT IS POSSIBLE THAT THE GRIDS IN THE CLIMATE MODEL IS
C-TF  MUCH SMALLER THAN THOSE IN THE PHOTOCHEMICAL MODEL.
C-TF  IN THIS CASE, IS=1 AND IS1=0. THERE IS NO VALUE IN
C-TF  THE ZERO GRID. AND WE CAN SIMPLY ASSUME A CONSTANT
C-TF  MIXING RATIO FOR THESE GRIDS (BELOW THE 1ST GRID IN 
C-TF  THE PHOTOCHEMICAL MODEL).
        IF (IS .EQ. 1) THEN
         FNO2_NEW(j) = FNO2_PHOTO(IS)
         FSO2_NEW(j) = FSO2_PHOTO(IS)
         FI(1,j) = water(IS)
         FI(2,j) = CO2(IS)
         FI(3,j) = CH4(IS)
         FI(4,j) = O3(IS)
C         GO TO 472
        ENDIF
C
        IS1 = IS-1
C        
        istart = IS1
        IF (IS .GT. NREAD) GOTO 471
C
        DZI = temp_alt(IS)-temp_alt(IS1)
        DZJ = ALT(j) - temp_alt(IS1)
        FR = DZJ/DZI
        IF (ALT(j) .LT. top) THEN
           O3log1 = LOG(O3(IS1))
           O3log2 = LOG(O3(IS))
           O3log = (FR*O3log2) + ((1-FR)*O3log1)
         ELSE O3log = FI(4,j+1)
        END IF
C



        FI(4,j) = EXP(O3log)

        F1log1 = LOG(water(IS1))                ! This is water layer F1
        F1log2 = LOG(water(IS))
        F11log = (FR*F1log2) + ((1-FR)*F1log1)
        FI(1,j) = EXP(F1log)

        F2log1 = LOG(CO2(IS1))                  ! This is CO2 layer F2
        F2log2 = LOG(CO2(IS))
        F2log = (FR*F2log2) + ((1-FR)*F2log1)
        FI(2,j) = EXP(F2log)

        F3log1 = LOG(CH4(IS1))                  ! This is CH4 layer F3
        F3log2 = LOG(CO2(IS))
        F3log = (FR*F3log2) + ((1-FR)*F3log1)
        FI(3,j) = EXP(F3log)
C
C-TF  INTERPOLATE NO2
        Flog1 = LOG(FNO2_PHOTO(IS1))
        Flog2 = LOG(FNO2_PHOTO(IS))
        flog = (FR*Flog2) + ((1-FR)*Flog1)
        FNO2_NEW(j) = EXP(flog)
C
C-TF  INTERPOLATE SO2
        Flog1 = LOG(FSO2_PHOTO(IS1))
        Flog2 = LOG(FSO2_PHOTO(IS))
        flog = (FR*Flog2) + ((1-FR)*Flog1)
        FSO2_NEW(j) = EXP(flog)
C

C
      END DO
471   CONTINUE

      DO k =1, ND
       if (ALT(k) .GE. temp_alt(JCOLD)) jcold_new = k
      END DO

      JCOLD = jcold_new

      RETURN
      END
