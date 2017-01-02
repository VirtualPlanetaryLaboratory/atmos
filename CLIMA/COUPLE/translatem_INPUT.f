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

      DIMENSION T(ND), FI(4,ND), ALTF(ND), DALTF(ND-1), TF(ND), PF(ND)
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FNO2
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     & 		     ALT(ND)
      COMMON/PRESSURE/P(ND),PLOG(ND)

      COMMON/ MLAWERI/  layers, numspec, newalt(ND), TempT(0:NLAYERS), 
     & 			Pres(0:NLAYERS), gasses(7, 0:NLAYERS), COLDEP(ND)


      DIMENSION ALTM(ND),PM(ND),TM(ND)
      DIMENSION FH2OM(ND),FCO2M(ND),FO3M(ND),FN2OM(ND),
     &          FCOM(ND),FCH4M(ND),FO2M(ND)
      DIMENSION FH2ON(ND),FCO2N(ND),FO3N(ND),FN2ON(ND),
     &          FCON(ND),FCH4N(ND),FO2N(ND)

      OPEN(UNIT=1, FILE='INPUT_RRTM.MLS')
   10 FORMAT(4x,F6.3,2x,F8.3,2x,F7.3)
   11 FORMAT(E10.4,E10.4,E10.4,E10.4,E10.4,E10.4,E10.4)
      DO J=ND,1,-1
        READ(1,10) ALTM(J),PM(J),TM(J)
        READ(1,11) FH2OM(J),FCO2M(J),FO3M(J),FN2OM(J),
     &             FCOM(J),FCH4M(J),FO2M(J)
      END DO


C        Assign mixing ratios to proper placement in array

        do i = 1, 7
	  select case (i)          
	   case (1) 		! water
             do n = 0, NLAYERS
        	gasses(1, n) = FH2OM(ND-n)/1.E6
             end do 
           case (2) 		! CO2
             do n = 0, NLAYERS
        	gasses(2, n) = FCO2
             end do
           case (3)		! ozone
             do n = 0, NLAYERS
        	gasses(3, n) = FO3M(ND-n)/1.E6
             end do
           case (4) 
              !OR READ from a different file somewhere; or a data block. N2O.
              do n = 0, NLAYERS
        	gasses(4, n) = FN2OM(ND-n)/1.E6
              end do
           case (5)
              !OR ditto from above. CO. 
              do n = 0, NLAYERS
        	gasses(5, n) = FCOM(ND-n)/1.E6
              end do
           case (6)		! methane
             do n = 0, NLAYERS
        	gasses(6, n) = FCH4
             end do
           case (7)
             do n = 0, NLAYERS
        	gasses(7, n) = FO2M(ND-n)/1.E6
             end do
        end select
        end do

C       Translate T grid:
C		T at Flux grid points
c	do j = 1, ND-1
c	  TF(j) = (T(j)+T(j+1))/2.
c	end do
c	TF(ND) = T(ND)

c-jdh Modified for comparision with Mlawer's MLS profile
        TF(1) = 272.32
	TF(ND) = 294.21
	do j = 2, ND-1
	  TF(j) = T(J-1) + 0.75*T(J) - 0.25*T(J+1) - 0.5*TF(J-1)
	end do


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
c     water = 0.5 * (FI(1,J) + FI(1,J+1))
      water = 0.5 * (FH2OM(J)/1.E6 + FH2OM(J+1)/1.E6)
      AM = 18.*water + DM*(1. - water)
      BMG = BKM/AM
      ALTF(J) = ALTF(J+1) + BMG*TAvg*ALOG(PF(J+1)/PF(J))*1.E-5
      DALTF(J) = ALTF(J) - ALTF(J+1)
      ENDDO

C       Calculate total atmospheric column depth
C	and altitude at flux points (midway between T points).	
	do ij = 0, NLAYERS
           newalt(ij) = ALTF(ND-ij)
	end do
        
	muH = 1.67E-24
c       muatm = ((FAR*(40.)*muH) + (FO2*(32.)*muH) + (FN2*(28.)*muH))
c       muatm = ((FN2/(FN2+FAR))*28. + (FAR/(FN2+FAR))*40.)*muH
        muatm = 28.*muH
        muatm = muatm*1E-3
 	
        do j = 1, ND
           Ncol = ((Pres(j-1)-Pres(j))*1.E2)/((G/100.)*muatm)
C		convert from m^-2 to cm^-2
c         COLDEP(j) = Ncol*1.E-4
          COLDEP(j) = Ncol*1.E-4*FN2
        end do

        layers = NLAYERS
	numspec = 7
	
	return
 
        END
