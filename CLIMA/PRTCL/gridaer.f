
        SUBROUTINE GRIDAER(ICOUPLE, IHAZE)
        INCLUDE 'CLIMA/INCLUDE/header.inc'
       PARAMETER (NZ=100)
        DIMENSION WFALL(NZ),TAUSED(NZ),CONVER(NZ),TAUEDD(NZ),TAUC(NZ)
        REAL Z(NZ),NPART(NZ),RPART(NZ)
        COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     2    ALT(ND)

C     in this case NZ = 100 is appropriate since the hcaer files have 100 layers

      IF (IHAZE.eq.0) THEN
      print *, 'atmosphere has no haze'
      OPEN(unit=101,file="CLIMA/IO/hcaer_nohaze.out")
      ENDIF

      IF (ICOUPLE.eq.0.and.IHAZE.eq.1) THEN
      print *, 'reading in uncoupled haze file'
      OPEN(unit=101,file="CLIMA/IO/hcaer.out")
      ENDIF

      IF (ICOUPLE.eq.1.and.IHAZE.eq.1) THEN
      print *, 'reading in coupled haze file'
      OPEN(unit=101,file="COUPLE/hcaer.photoout.out")
      OPEN(unit=102,file="COUPLE/hcaer.climaout.out")
      ENDIF

      READ(101,*)
      READ(101,*)
      READ(101,*)
      READ(101,105) (Z(I),NPART(I),RPART(I),WFALL(I),TAUSED(I),
     2 TAUEDD(I),TAUC(I),CONVER(I),I=1,NZ)
c      print *, 'read it in'
      CLOSE(101)

      IF (ICOUPLE.eq.1) THEN
      WRITE(102,103)
      WRITE(102,104)
      WRITE(102,105) (Z(I),NPART(I),RPART(I),WFALL(I),
     2 TAUSED(I),TAUEDD(I),TAUC(I),CONVER(I),I=1,NZ)
      CLOSE(102)
      ENDIF

 103  FORMAT(1X,"HC AEROSOL PARAMETERS")
 104  FORMAT(/4X,"Z",8X,"AERSOL",5X,"RPART",6X,"WFALL",5X,"TAUSED",4X,
     2 "TAUEDD",4X,"TAUC",6X,"CONVER")
 105  FORMAT(1X,1P8E10.3)


C USE THIS PARTICLE DISTRIBUTION FOR FCO2=0.001 FCH4=0.001 
C        DATA NPART/7.7E-1, 2.05, 4.3, 14.5, 48.6, 47.6, 42.2, 40.3,
C     2  47.1, 72.1, 124., 216.6, 379.6, 664.8, 1163., 2026., 3448.,
C     3  5638./
C  Radius of particles in cm
C        DATA RPART/5.3E-5, 5.3E-5, 5.3E-5, 5.3E-5, 5.3E-5,
C     2  5.E-5, 4.5E-5, 3.8E-5, 3.1E-5, 2.3E-5, 1.7E-5, 
C     3  1.2E-5, 8.8E-6, 6.4E-6, 4.6E-6, 3.3E-6, 2.3E-6, 1.5E-6/
C
C USE THIS PARTICLE DISTRIBUTION FOR FCO2=0.002 FCH4=0.001 
C        DATA NPART/5.5E-1, 1.45, 3., 10.4, 37.2, 37.6, 33.9, 33.4,
C     2  41.4, 65.9, 114., 200.1, 350.7, 614.4, 1075., 1880., 3209.,
C     3  5167./
C  Radius of particles in cm
C        DATA RPART/4.6E-5, 4.6E-5, 4.6E-5, 4.6E-5, 4.6E-5,
C     2  4.4E-5, 3.9E-5, 3.3E-5, 2.6E-5, 1.95E-5, 1.4E-5, 
C     3  1.E-5, 7.4E-6, 5.4E-6, 3.9E-6, 2.8E-6, 1.96E-6, 1.2E-6/
C USE THIS PARTICLE DISTRIBUTION FOR FCO2=0.0025 FCH4=0.001 
C        DATA NPART/4.4E-1, 1.18, 2.5, 8.5, 31.5, 32.5, 29.6, 29.9,
C     2  38.4, 62.4, 108.5, 190., 333.1, 583.7, 1022., 1787., 3053.,
C     3  4851./
C  Radius of particles in cm
C        DATA RPART/4.3E-5, 4.3E-5, 4.3E-5, 4.3E-5, 4.3E-5,
C     2  4.1E-5, 3.6E-5, 3.E-5, 2.4E-5, 1.76E-5, 1.3E-5, 
C     3  9.E-6, 6.7E-6, 4.9E-6, 3.5E-6, 2.6E-6, 1.8E-6, 1.1E-6/


        DO I =1,NZ
         RPART(I) = RPART(I)
         Z(I) = Z(I)/1.e5 ! Convert from cm to km
        ENDDO 
        DO IL = 1,ND
          I = 1
          IF (ALT(IL).LE.Z(I)) THEN
            PARTICLES(IL) = 0.
            RAER(IL) = RPART(I)
            GOTO 4
          ENDIF
 20       IF (ALT(IL).GT.Z(I)) THEN
            I = I+1
            GOTO 20
          ENDIF
          PARTICLES(IL) = NPART(I-1)+((NPART(I)-NPART(I-1))/
     &    (Z(I)-Z(I-1)))*(ALT(IL)-Z(I-1))
          RAER(IL) = RPART(I-1)+((RPART(I)-RPART(I-1))/
     &    (Z(I)-Z(I-1)))*(ALT(IL)-Z(I-1))

    4   ENDDO
c        DO IL = 1,ND
c        PARTICLES(IL) = PARTICLES(IL)/100
c        ENDDO

        RETURN
        END
