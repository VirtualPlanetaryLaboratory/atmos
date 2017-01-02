C------------------------------------------------------------------
C   Code added 6/15/01, Kara Krelove, for integration of Mlawer's RRTM
C------------------------------------------------------------------
        SUBROUTINE TRANSLATEM(G,FI,T,PF,ND1,DM,BKM)

C    This subroutine is designed to correctly translate between SurfTem
C    and RRTM, including flipping of layer indices, and READing 
C    molecular abundances into the appropriate grid. 

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NLAYERS = ND-1, NS=7)
c     PARAMETER (NLAYERS = 51, NS=7)

        REAL muH, muatm, Ncol, newalt

      DIMENSION T(ND), FI(5,ND), ALTF(ND), DALTF(ND-1), TF(ND), PF(ND)
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     &     ALT(ND)
      COMMON/PRESSURE/P(ND),PLOG(ND)

      COMMON/ MLAWERI/  layers, numspec, newalt(ND), TempT(0:NLAYERS), 
     &     Pres(0:NLAYERS), gasses(7, 0:NLAYERS), COLDEP(ND)

C        Assign mixing ratios to proper placement in array

        do i = 1, 7
          select case (i)          
           case (1) ! water
             do n = 0, NLAYERS
                gasses(1, n) = FI(1, ND-n)
             end do 
           case (2) ! CO2
             do n = 0, NLAYERS
                gasses(2, n) = FI(2, ND-n)
             end do
           case (3) ! ozone
             do n = 0, NLAYERS
                gasses(3, n) = FI(4, ND-n)
             end do
           case (4) 
              !OR READ from a different file somewhere; or a data block. N2O.
              do n = 0, NLAYERS
                gasses(4, n) = 1E-20
              end do
           case (5)  ! carbon monoxide
              !OR ditto from above. CO. 
              do n = 0, NLAYERS
                gasses(5, n) = 1E-20
              end do
           case (6)  ! methane
             do n = 0, NLAYERS
                gasses(6, n) = FI(3, ND-n)
             end do
           case (7)
             do n = 0, NLAYERS
                gasses(7, n) = FO2
             end do
        end select
        end do

C       Translate T grid:
C                T at Flux grid points
        do j = 1, ND-1
          TF(j) = (T(j)+T(j+1))/2.
        end do
        TF(ND) = T(ND)

        do i = 0, NLAYERS
          TempT(i) = TF(ND-i)
        end do

C        Translate P grid:
        do k = 0, NLAYERS
           Pres(k) = PF(ND-k)*1000.
        end do

C        Calculate ALT at flux points (currently at T points)

      ALTF(ND) = 0.
      DO JS=1,ND1
      J = ND - JS
      TAvg = 0.5*(TF(J) + TF(J+1))
      water = 0.5 * (FI(1,J) + FI(1,J+1))
      AM = 18.*water + DM*(1. - water)
      BMG = BKM/AM
      ALTF(J) = ALTF(J+1) + BMG*TAvg*ALOG(PF(J+1)/PF(J))*1.E-5
      DALTF(J) = ALTF(J) - ALTF(J+1)
      ENDDO

C       Calculate total atmospheric column depth
C        and altitude at flux points (midway between T points).        
        do ij = 0, NLAYERS
           newalt(ij) = ALTF(ND-ij)
        end do
        
        muH = 1.67E-24
c-jdh        muatm = ((FAR*(40.)*muH) + (FO2*(32.)*muH) + (FN2*(28.)*muH))
        muatm = 28.*muH
        muatm = muatm*1E-3
        
        do j = 1, ND
          Ncol = ((Pres(j-1)-Pres(j))*1.E2)/((G/100.)*muatm)
C                convert from m^-2 to cm^-2
c-jdh          COLDEP(j) = Ncol*1.E-4
          COLDEP(j) = Ncol*1.E-4*FN2
        end do

        layers = NLAYERS
        numspec = 7
        
        return
 
        END
