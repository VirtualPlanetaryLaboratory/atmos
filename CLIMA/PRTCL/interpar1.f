     
        SUBROUTINE INTERPAR1(RAER) 
        INCLUDE 'CLIMA/INCLUDE/header.inc'
        PARAMETER (NF=55, NSOL=38, NGS=8)
        DIMENSION RAER(ND)
        COMMON/HYDROCARB/Qextirst(73,55),w0irst(73,55),
     &  girst(73,55),Qextsolst(73,38),w0solst(73,38),gsolst(73,38),
     &  radstand(73)

        COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1)

C        COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
C     &  ASYAIR(NF,ND-1),WEIGHT(8,3), xkappa(8,12,55,8,3),IO3,
C     &  QEXTIR(NF,ND-1)
        COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
     &  ASYAIR(NF,ND-1),IO3,QEXTIR(NF,ND-1)
     
****** IR fitting
        DO i = 1,55 !number of wavlengths
        DO j = 1,ND-1 !number of layers
        DO k = 1,72    !number of particle bins - 1
        IF ((RAER(j).GE.radstand(k)).and.(RAER(j).LT.radstand(k+1)))
     &   THEN
                     
         
        dra = RAER(j) - radstand(k)
        dr = radstand(k+1) - radstand(k)
        QEXTIR(i,j)=Qextirst(k,i)+((Qextirst(k+1,i)-Qextirst(k,i))
     &  /dr)*dra
        OMG0AIR(i,j)=w0irst(k,i)+((w0irst(k+1,i)-w0irst(k,i))
     &  /dr)*dra
        ASYAIR(i,j)=girst(k,i)+((girst(k+1,i)-girst(k,i))
     &  /dr)*dra
C        ELSE
C        QEXTIR(i,j) = Qextirst(1,i)
C        OMG0AIR(i,j) = w0irst(1,i)
C        ASYAIR(i,j) = girst(1,i)
        ENDIF
        ENDDO
        ENDDO
        ENDDO
****** SOLAR fitting 
        DO i = 1,38 !number of wavelengths 
        DO j = 1,ND-1 !number of layers
        DO k = 1,72   ! number of particle bins - 1
        IF ((RAER(j).GE.radstand(k)).and.(RAER(j).LT.radstand(k+1)))
     &   THEN
        dra = RAER(j) - radstand(k)
        dr = radstand(k+1) - radstand(k)
        QEXT(i,j)=Qextsolst(k,i)+((Qextsolst(k+1,i)-Qextsolst(k,i))
     &  /dr)*dra
        OMG0A(i,j)=w0solst(k,i)+((w0solst(k+1,i)-w0solst(k,i))
     &  /dr)*dra
        ASYA(i,j)=gsolst(k,i)+((gsolst(k+1,i)-gsolst(k,i))
     &  /dr)*dra
        ENDIF
C        IF (RAER(j).LT.radstand(1)) THEN
C        QEXT(i,j) = Qextsolst(1,i)
C        OMG0A(i,j) = w0solst(1,i)
C        ASYA(i,j) = gsolst(1,i)
C        ENDIF
        ENDDO
        ENDDO
        ENDDO
        RETURN
        END
