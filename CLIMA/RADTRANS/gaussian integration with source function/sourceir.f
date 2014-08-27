
      SUBROUTINE SOURCEIR(SRFALB,ASY,TAULAM,OMG0,
     &  FUP,FDN,BPLANCK,TAUTOP) 
!Written by Jim Kasting. 
! Source Function modification to two-stream hemispheric mean approximation by Ramses Ramirez and Ravi Kopparapu 5/10/2012

c
c jfk 6/25/08  Note: I added TAUTOP (the optical depth above the top
c     of the grid) to the call sequence    
C

      INCLUDE '../header.inc'


C

      PARAMETER(NZ=ND-1, NZ2=2*NZ)  ! NZ is useful because you can do NZ+1, for last layer at the surface
      REAL NORM,B1
      DIMENSION TAU(NZ), G(NZ), GAM1(NZ), GAM2(NZ),
     1  ALAM(NZ), CGAM(NZ), E1(NZ), E2(NZ), E3(NZ), E4(NZ),
     2  CP0(NZ), CPB(NZ), CM0(NZ), CMB(NZ), Y1(NZ), Y2(NZ),
     3  W0(NZ)
      DIMENSION A(NZ2), B(NZ2), D(NZ2), E(NZ2), Y(NZ2)
      DIMENSION ASY(NZ),TAULAM(NZ), OMG0(NZ),
     &     FUP(ND,3),FDN(ND,3),BPLANCK(ND)
      DIMENSION TERMUP1(NZ),TERMUP2(NZ),TERMUP3(NZ),TERMUP4(NZ),
     & TERMUP5(NZ),TERMDN1(NZ),TERMDN2(NZ),TERMDN3(NZ),TERMDN4(NZ),
     & TERMDN5(NZ), SIG1(NZ), SIG2(NZ),
     & ALPHA1(NZ), ALPHA2(NZ), AJ(NZ), AK(NZ), CG(NZ), H(NZ)

C
      REAL AMU
      DIMENSION AMU(3)
         DATA AMU/0.23861, 0.66120, 0.93246/  ! zenith angle weights for 3 gauss points
!          DATA AMU/0.46791, 0.36076, .17132/
      DATA C,HP,BK,SIGMA,PI,SM/3.E10, 6.63E-27, 1.38E-16, 5.67E-5,
     2  3.14159, 1.67E-24/
     


       DO I = 1,3  ! For all three gaussian weights


      U1 = AMU(I)

      
      ALB = SRFALB
      PI = 3.14159
      emissivity = 1.
      BCON = 2.*HP/C/C
      HK = HP/BK


!      U1 = 0.5   ! mu_1 , the zenith angle (due to diffuse scattering)
      U1M = 1./U1
      NZM1 = NZ - 1
      NZP1 = NZ + 1
      MZ2 = NZ2
      NORM = 2*U1*PI

      DO 1 N=1,NZ
      TAU(N) = TAULAM(N)
      W0(N) = OMG0(N)
      G(N) = ASY(N)
 1    CONTINUE

      DO 2 N= 1, NZ
      GAM1(N) = 2 - W0(N)*(1.+G(N)) 
      GAM2(N) = W0(N)*(1.-G(N))
      ALAM(N) = SQRT(GAM1(N)*GAM1(N) - GAM2(N)*GAM2(N))  ! wavelength expression
      CGAM(N) = (GAM1(N) - ALAM(N))/GAM2(N)
      EMLT = EXP(-ALAM(N)*TAU(N))
      E1(N) = 1. + CGAM(N)*EMLT
      E2(N) = 1. - CGAM(N)*EMLT
      E3(N) = CGAM(N) + EMLT
      E4(N) = CGAM(N) - EMLT

 2    CONTINUE



C   Top of atmosphere
      A(1) = 0.
      B(1) = E1(1)
      D(1) = -E2(1)
C   Odd coefficients
      DO 3 N=1,NZM1
      L = 2*N + 1
      A(L) = E2(N)*E3(N) - E4(N)*E1(N)
      B(L) = E1(N)*E1(N+1) - E3(N)*E3(N+1)
      D(L) = E3(N)*E4(N+1) - E1(N)*E2(N+1)

 3    CONTINUE
C
C   Even coefficients
      DO 4 N=1,NZM1
      L = 2*N
      A(L) = E2(N+1)*E1(N) - E3(N)*E4(N+1)
      B(L) = E2(N)*E2(N+1) - E4(N)*E4(N+1)
      D(L) = E1(N+1)*E4(N+1) - E2(N+1)*E3(N+1)

 4    CONTINUE

C
C   Bottom of atmosphere
      A(NZ2) = E1(NZ) - ALB*E3(NZ)
      B(NZ2) = E2(NZ) - ALB*E4(NZ)
      D(NZ2) = 0.
      B(1) = E1(1)


C   Now, set up the RHS of the equation:
C

      NORM = 2*U1*PI
      DO 6 N=1,NZ
C

      B0n = BPLANCK(N)
      B1 = BPLANCK(N+1)
      B1n = (B1-B0n)/TAU(N)
      CP0(N) = NORM*(B0n + B1n*(1/(GAM1(N)+GAM2(N))))
      CPB(N) = NORM*(B0n + B1n*(TAU(N)+1/(GAM1(N)+GAM2(N))))
      CM0(N) = NORM*(B0n + B1n*(-1/(GAM1(N)+GAM2(N))))
      CMB(N) = NORM*(B0n + B1n*(TAU(N)-1/(GAM1(N)+GAM2(N))))

    6 CONTINUE


      SSFC = emissivity*PI*BPLANCK(ND)

      F0M0 = 0.
      E(1) = - CM0(1) + F0M0

      DO 7 N=1,NZM1
      L = 2*N + 1
      E(L) = (CP0(N+1)-CPB(N))*E3(N) + (CMB(N)-CM0(N+1))*E1(N)
   
  7   CONTINUE

     
C   Even coefficients
      DO 8 N=1,NZM1
      L = 2*N
      E(L) = (CP0(N+1)-CPB(N))*E2(N+1) - (CM0(N+1)-CMB(N))*E4(N+1)

 8    CONTINUE

      E(NZ2) = SSFC - CPB(NZ) + ALB*CMB(NZ)
C     Call the tridiagonal solver.  Use Numerical Recipes for now.
                  
       CALL TRIDAG(A,B,D,E,Y,MZ2)

c----------------------------------------------------------- ALL THE NEW SOURCE FUNCTION STUFF c-rr 5/10/2012



      DO 10 N=1,NZ   ! Need all 100 layers of Ys and other parameters! 
      L = 2*N   ! for even Ys
      L1 = L-1  ! for odd Ys

!     Calculate Bs to put into below source function parameter arrays
      B0n = BPLANCK(N)
      B1 = BPLANCK(N+1)
      B1n = (B1-B0n)/TAU(N)

!     Calculating Source Function Parameters at every layer from layers 1 to 100
      Y1(N) = Y(L1)
      Y2(N) = Y(L)

      CG(N)=(Y1(N) + Y2(N))*((1./U1)-ALAM(N))
      H(N) = (Y1(N) - Y2(N))*CGAM(N)*(ALAM(N) + (1./U1))
      AJ(N) = (Y1(N) + Y2(N))*CGAM(N)*(ALAM(N) + (1./U1))
      AK(N) = (Y1(N) - Y2(N))*((1./U1)-(ALAM(N)))
      ALPHA1(N) = 2.*PI*B0n + 2*PI*B1n/(GAM1(N)+GAM2(N))
     &           -2.*PI*U1*B1n
      ALPHA2(N) = 2.*PI*B1n
      SIG1(N)=  2.*PI*B0n - 2*PI*B1n/(GAM1(N)+GAM2(N))-2.*PI*U1*B1n
      SIG2(N) = 2.*PI*B1n
! We changed B1n to B1, for now...
 

! MATLABABLE
!      print 36573,'Y1=',Y1(N),'Y2=',Y2(N),'G=',CG(N),'H=',H(N),
!     & 'J=',AJ(N),
!     & 'K=', AK(N), 'ALPHA1=', ALPHA1(N), 'ALPHA2=', ALPHA2(N), 
!     & 'SIG1=', SIG1(N), 'SIG2=', SIG2(N), 'B1n=', B1n,'B0n=', B0n,
!     & 'GAM1=', GAM1(N), 'GAM2=', GAM2(N), 'ALAM=', ALAM(N), 
!     &  'U1=', U1, 'CGAM=', CGAM(N), 'PI=', PI, 'TAU=', TAU(N), 
!     & 'N=',N
!      pause

36573    format(19(a9,1pe20.12,/),5x,a3,i2)

 10   CONTINUE 

!      print *, 'Y1=', Y1(102), 'L=', L, 'L1=', L1
!      pause


!-------------------BOUNDARY CONDITIONS AT TOP----------------------------------------------------------
       
      FDN(1,I) = 0.! At the top, the intensity and flux are zero.
! NOTE: OUR PLANCK FUNCTION IS IN FREQUENCY SPACE. SO INTENSITIES ARE IN UNITS OF ERGS/CM^2 (FLUX/WAVENGTH OR (ERGS/CM^2/SEC)/(CM^-^-1). 
! These intensities (in ergs/cm^2) are the FUPS and FDNS that are multiplied by the weight factors (W(I) in per sec) to get the wavelength integrated
! fluxes for that bin, in ergs/cm^2/sec

!------------------------------------------------------------------------------------------------------


! Calculating DOWNWARD intensities starting from the top of the atmosphere and moving downward
      DO 11 N=1,NZ   ! For all NZ layers, and get final value at NZ+1 

      TERMDN1(N) = FDN(N,I)*EXP(-TAU(N)/U1)  
      TERMDN2(N) = AJ(N)/(ALAM(N)*U1 + 1.)*(1. - 
     &   exp(-TAU(N)*(ALAM(N) + 1./U1)))
!      TERMDN2(N) = AJ(N)/(ALAM(N)*U1 + 1.)*(-1. + 
!     &   exp(TAU(N)*(ALAM(N) + 1./U1)))
      TERMDN3(N) = (AK(N)/((ALAM(N)*U1) -1.))*exp(-TAU(N)/U1)

     & -(AK(N)/((ALAM(N)*U1) -1.))*exp(-TAU(N)*ALAM(N))

      TERMDN4(N) = SIG1(N) - SIG1(N)*exp(-TAU(N)/U1)
      TERMDN5(N) = SIG2(N)*U1*EXP(-TAU(N)/U1) + SIG2(N)*TAU(N)
     &    -SIG2(N)*U1
      FDN(N+1,I) =(TERMDN1(N) + TERMDN2(N)+TERMDN3(N)
     &	  + TERMDN4(N) + TERMDN5(N))  ! Convert downward intensities to downward fluxes (I = F/2*pi*mu). STILL IN ERGS/CM^2


!         print 3333, 'TAU=', TAU(N), 'T1=', TERMDN1(N), 'T2=', 
!     &    TERMDN2(N),
!     &   'T3=', TERMDN3(N), 'T4=', TERMDN4(N), 'T5=', TERMDN5(N),
!     &    'FDN=', FDN(N+1)

!      print *, 'FDNA=', FDN(N)
!      pause

! 3333 	format(8(a7,1pe20.12,/))

 11   CONTINUE   



!----------------------------------BOUNDARY CONDITIONS AT TOP---------------------------------------------------------
        

! NOTE: OUR PLANCK FUNCTION IS IN FREQUENCY SPACE. SO INTENSITIES ARE IN UNITS OF ERGS/CM^2 (FLUX/WAVENGTH OR (ERGS/CM^2/SEC)/(CM^-^-1). 
! NOTE: OUR PLANCK FUNCTION IS IN FREQUENCY SPACE. SO INTENSITIES ARE IN UNITS OF ERGS/CM^2 (FLUX/WAVENGTH OR (ERGS/CM^2/SEC)/(CM^-^-1). 
! These intensities (in ergs/cm^2) are the FUPS and FDNS that are multiplied by the weight factors (W(I) in per sec) to get the wavelength integrated
! fluxes for that bin, in ergs/cm^2/sec


        FUP(ND,I) = BPLANCK(ND)*2*PI
           ! Defining the upward intensity and flux at the surface, and moving upward. 
!        print *, 'FUPND=',FUP(ND)
!        pause


!       Y values go from 1 - 100!.. NOT 101!

!---------------------------------------------------------------------------------------------------------------------------------

! Calculating UPWARD intensities starting from the bottom of the atmosphere going up

      DO 12 N = NZ, 1, -1  !For all NZ layers, and get final value at NZ + 1 = 1
!    We need TAU(ND).. If do NZ, not using the BC
!        print *, U1, N
!        pause

        TERMUP1(N) = FUP(N+1,I)*exp(-TAU(N)/U1)
        AA = CG(N)/(ALAM(N)*U1 - 1.)
        BB = (exp(-TAU(N)/U1)-exp(-TAU(N)*ALAM(N)))
        TERMUP2(N) = AA*BB
        TERMUP3(N) =  (H(N)/(ALAM(N)*U1 + 1.))*(1. - 
     &                EXP(-TAU(N)*ALAM(N)- TAU(N)/U1))
        TERMUP4(N) = ALPHA1(N) - ALPHA1(N)*EXP(-TAU(N)/U1) 
        TERMUP5(N) = ALPHA2(N)*(U1 - (TAU(N) + U1)*exp(-TAU(N)/U1))
        FUP(N,I)     = (TERMUP1(N) +TERMUP2(N)+ TERMUP3(N)+ 
     &                TERMUP4(N)+TERMUP5(N))
          ! Convert upward intensities to upward fluxes (I = F/2pi*mu). STILL IN ERGS/CM^2

!              print *, 'LAYER=', N
!			  pause
!        print 3456, 'TAU=', TAU(N), 'T1=', TERMUP1(N), 'T2=', 
!     &    TERMUP2(N),
!     &   'T3=', TERMUP3(N), 'T4=', TERMUP4(N), 'T5=', TERMUP5(N),
!     &    'FUP=', FUP(N), 
!     &    'G=', CG(N), 'LAM=', ALAM(N), 'U1=', U1, 
!     &    'ALPHA1=', ALPHA1(N), 'ALPHA2=', ALPHA2(N), 
!     &     'Y1=', Y1(N), 'Y2=', Y2(N)
!         print *,FUP(N+1),exp(-TAU(N)/U1),N
   
!          print *,TAU(N),N

!         pause
3456   format(15(a10,1pe20.12,/), (a10, i3),/)
       
 12     CONTINUE

!         print *, 'U1=', U1
!                 pause
         ENDDO ! Ends gaussian loop
       
!        print *, 'ALPHA1=',ALPHA1(5)

      RETURN
      END

