     
      SUBROUTINE INPUT_INTERP(temp_alt, water, O3, Jcold, T, FI
     |   , FNO2_PHOTO, FNO2_NEW)
C-KK This subroutine interpolates the (alt),water, ozone file that is READ in 
C-KK from atm_chem.

      INCLUDE '../INCLUDE/header.inc'
      PARAMETER (NZ=64)
      DIMENSION O3(NZ), water(NZ), temp_alt(NZ),T(ND),FI(4,ND)
     |   ,FNO2_PHOTO(NZ),FNO2_NEW(ND)
C
      CHARACTER :: DIRINOUT*2,DIRDATA*4
      COMMON/DIR/DIRINOUT,DIRDATA
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     & 	     ALT(ND)
      COMMON/PRESSURE/P(ND),PLOG(ND)

C-KK   both indices starting at ground level
       istart = 1
       j = ND
       FI(4,j) = O3(istart)
       FI(1,j) = water(istart)
       FNO2_NEW(j) = FNO2_PHOTO(1)
C
       istart = istart!+1
       top = temp_alt(NZ)
C
      DO j = ND-1, 1, -1
        DO i = istart, NZ
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
         FI(1,j) = water(IS)
         FI(4,j) = O3(IS)
         GO TO 472
        ENDIF
C
        IS1 = IS-1
C        
        istart = IS1
        IF (IS .GT. NZ) GOTO 471
C
        DZI = temp_alt(IS)-temp_alt(IS1)
        DZJ = ALT(j) - temp_alt(IS1)
        FR = DZJ/DZI
        IF (ALT(j) .LT. top) THEN
           O3log1 = ALOG(O3(IS1))
           O3log2 = ALOG(O3(IS))
           O3log = (FR*O3log2) + ((1-FR)*O3log1)
        ELSE O3log = FI(4,j+1)
        END IF
C
        FI(4,j) = EXP(O3log)
C
        Flog1 = ALOG(water(IS1))
        Flog2 = ALOG(water(IS))
        flog = (FR*Flog2) + ((1-FR)*Flog1)
        FI(1,j) = EXP(flog)
C
C-TF  INTERPOLATE NO2
        Flog1 = ALOG(FNO2_PHOTO(IS1))
        Flog2 = ALOG(FNO2_PHOTO(IS))
        flog = (FR*Flog2) + ((1-FR)*Flog1)
        FNO2_NEW(j) = EXP(flog)
C
472   CONTINUE
C      PRINT 99913, J, IS, ALT(J), FNO2_NEW(j), FR
99913 FORMAT("in input_inter.f:",2I5, 1P10E10.2)
C
      END DO
471   CONTINUE

      DO k =1, ND
       if (ALT(k) .GE. temp_alt(JCOLD)) jcold_new = k
      END DO

      JCOLD = jcold_new

      RETURN
      END
