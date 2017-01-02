     
      SUBROUTINE INPUT_INTERP(temp_alt,water,O3,CH4,CO2,ethane,Jcold,T,
     &           FI)
C-KK This subroutine interpolates the (alt),water, ozone file that is READ in 
C-KK from atm_chem.
!FI is matrix of mixing ratios
!strip down to be simpler so that it takes in only one gas concentration
!spits out profile of gas concentration
!FI_temp = one dimensional vector -> manually set the particular column of FI
!grid to the values in that vector
!CALL INPUT_INTERP(water, FI's column for water)

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NZ=200)
      DIMENSION O3(NZ),water(NZ),CH4(NZ),CO2(NZ),temp_alt(NZ),T(ND)
      DIMENSION ethane(NZ)
      DIMENSION FI(5,ND)
      CHARACTER :: DIRINOUT*8,DIRDATA*10
      COMMON/DIR/DIRINOUT,DIRDATA
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     &               ALT(ND)
      COMMON/PRESSURE/P(ND),PLOG(ND)
      print *,'read in defs'
C-KK   both indices starting at ground level
        T = T ! EWS - note unused dummy argument
        istart = 1
        j = ND
        FI(1,j) = water(istart)     ! Not defined in mixing_ratios.dat
        FI(2,j) = CO2(istart)
        FI(3,j) = CH4(istart)
        FI(4,j) = O3(istart)        ! Not defined in mixing_ratios.dat
        FI(5,j) = ethane(istart)
        
        istart = istart+1
        top = temp_alt(NZ)
        
        DO j = ND-1, 1, -1 !starting at ND-1, which is second to last layer, which in climate code is the surface 
           !clima and photochem alt codes are reversed from each other (counting by increments of -1)
         DO i = istart, NZ
           IS = i
           if (temp_alt(i) .GT. ALT(j)) GOTO 350
         END DO
  350  CONTINUE
         IS1 = IS-1
          istart = IS1
         IF (IS .GT. NZ) GOTO 471
         DZI = temp_alt(IS)-temp_alt(IS1)
         DZJ = ALT(j) - temp_alt(IS1)
         FR = DZJ/DZI
        

        IF (ALT(j) .LT. top) THEN !is this a bug that this layer logic is not built into the water? (top = top of photochem grid)
          F4log1 = ALOG(O3(IS1))
          F4log2 = ALOG(O3(IS))
          F4log = (FR*F4log2) + ((1-FR)*F4log1)        
         ELSE F4log = FI(4,j+1)
        END IF 
        FI(4,j) = EXP(F4log)
        F1log1 = ALOG(water(IS1))                     ! This is water layer F1
        F1log2 = ALOG(water(IS))
        F1log = (FR*F1log2) + ((1-FR)*F1log1)
        FI(1,j) = EXP(F1log)
        F2log1 = ALOG(CO2(IS1))                       ! This is CO2 layer F2
        F2log2 = ALOG(CO2(IS))
        F2log = (FR*F2log2) + ((1-FR)*F2log1)
        FI(2,j) = EXP(F2log)
        F3log1 = ALOG(CH4(IS1))                       ! This is CH4 layer F3
        F3log2 = ALOG(CH4(IS))
        F3log = (FR*F3log2) + ((1-FR)*F3log1)
        FI(3,j) = EXP(F3log)
        F5log1 = ALOG(ethane(IS1))                    ! This is the ethane layer F5 (gna)
        F5log2 = ALOG(ethane(IS))
        F5log = (FR*F5log2) + ((1-FR)*F5log1)
        FI(5,j) = EXP(F5log)
        
        END DO

  471  CONTINUE

        DO k =1, ND !k is index; looping over number of layers in climate code (ND = number of layers in clima)
           !goes from top of atmosphere since clima counts from top down
          if (ALT(k) .GE. temp_alt(JCOLD)) jcold_new = k !redefines cold trap until it matches the real cold trap
         END DO
        
        JCOLD = jcold_new

        !take out water statements and generalize the o3 statement to be for any gas.  
        !Will interpolate that specific gas to the climate code.  
        !We could probably use the same function to interpolate the number density of the particles.  
        !Then we can call this a bunch of times for each thing we're interested in.

        RETURN
        END
