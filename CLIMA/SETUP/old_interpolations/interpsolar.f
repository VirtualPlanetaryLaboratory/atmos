		SUBROUTINE interpsolar(Temp,p1,j,Fwater)  ! Written by Ramses Ramirez and Ravi Kopparapu
c-rr    i is wavelength counter and j is the altitude counter used in solar, respectively.

	INCLUDE '../header.inc'
	PARAMETER(NSOL=38,IK=16, NGS=5, NS=3, NS1 = NS+1 )
	REAL TempR,tempg(8),press(8), FC(6), FW(9),p1,kmatrix_sol,Temp,
     &  Fwater, kappa_sol, AK(16), AN(16)
	INTEGER i,j
	DATA tempg/100., 150., 200., 250., 300.,350.,400., 600./  ! 8 temps 
	DATA press/.00001, .0001,.001,.01, 0.1, 1., 10., 100./  ! 8 pressures
        DATA FC/1.e-04, 1.e-03, 1.e-02, 1.e-01, 1., 3.3e-04/  ! 6 FCO2s but only first 5 are used for interpolation. Last is for T>400K
        DATA FW/1.e-08, 1.e-07, 1.e-06, 1.e-05, 1.e-04, 1.e-03,
     &         1.e-02, 1.e-01, 1./ ! 9 FH2Os
	COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL), 
     &  kappa_sol(NSOL,8,9,6,8,16) ! Added new kappa matrix for mixed CO2-H2O coefficients. 
!     Added O3 absorption coefficient vector (plus c and d vectors), first 14 terms. 3/19/2012 ! 38 layers, 7 temps, 9 FH2Os, 6 FCO2s, 9 pressures, 16 coefficients



      COMMON/DATASOLAR/weightco2_h2oSOL(16), weights(2,NSOL,IK)   ! new common block for weights and interpolated coefficients for mixed H20-CO2, methane. 3/27/2012
     & ,kmatrix_sol(NSOL,IK)

c-rr    Created new common block solar data and put in readsol.f and solar.f as well
!       COMMON/SOLARDATA/weightco2_h2oSOL(16), weights(2,NSOL,IK)
!     & ,kmatrix_sol(NSOL,IK) ! c-rr Added weightco2_h2O data array and weights array that combines co2_h2o and Ch4 weight arrays 3/19/2012  ! Added in 16 term weight array


      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND),
     &   FH22,FNC

c  	Assigning TempR to temp (due to variable name overuse)
	TempR=Temp

              

	   Do L=1,8 ! For 8 pressures
		LS=L 
		 
	     If (p1.le.press(L))then
               exit!If an individial p greater than curent p set Ls=6 
             else
	     LS = 9
             endif
	   enddo
	   

           Do M=1,8  ! For 8 temperature intervals 
	        MS=M
  
	      If (TempR.lt.tempg(M))then 
                exit
              else  
	      MS=9
              endif
	   enddo


         DO IW = 1,9  ! For 9 FH2os
                 IWS = IW

	      If (Fwater.lt.FW(IW))then
                exit  
              else
	      IWS=10
              endif 
	   enddo

	DO IC = 1,5  ! For 6 FCO2s
                 ICS = IC

	      If (FCO2.lt.FC(IC))then 
               exit
              else  
	      ICS=7
              endif 
	   enddo



        IF (LS.eq.1)THEN
        LS1=1
        LS = 2        
	PLOG=LOG10(press(LS))  ! higher grid point
	PLOG1=LOG10(press(LS1)) ! lower grid point
	P1L=LOG10(p1)
!	FY=(P1L-PLOG1)/(PLOG-PLOG1)
        ENDIF


        IF(LS.le.7)THEN
	LS1 = LS-1
        PLOG = LOG10(press(LS)) ! higher grid point
        PLOG1 = LOG10(press(LS1)) ! lower grid point
        P1L=LOG10(p1)
        ENDIF


        
        IF(LS.ge.8) THEN
        LS = 8		! P1 is at or above the highest grid pressure
	LS1 = 7
        PLOG = LOG10(press(LS))
        PLOG1 = LOG10(press(LS1))
        P1L = LOG10(p1)
	ENDIF

C

	IF(MS.eq.1)THEN
	MS1=MS
        MS = 2		!TempR is at or below the lowest grid Temperature
	TempR = tempg(MS)
        TempL = tempg(MS1)
        ! Change all tempg(MS1) to TempL

	ELSEIf(MS.ge.8) THEN
	MS=8	!TempR is at or above the highest grid Temperature
	MS1=MS-1		
        TempR = tempg(MS)
        TempL = tempg(MS1)
        ELSE
	MS1=MS-1  !This calculates MS1 for values between MS 2-7
        TempL = tempg(MS1) ! high grid point is already tempg(MS) which is done in DEN. TempR is just input Temp
	ENDIF



        IF ((IWS.eq.1).and.(Fwater.lt.1.e-08))THEN
        IWS=2
        IWS1 = 1
        WLOG = LOG10(FW(IWS))  ! higher grid point
        WLOG1 = LOG10(FW(IWS1)) ! lower grid point
        W1L = WLOG1

	
        ELSEIF(IWS.eq.1) THEN
        IWS = 2
        IWS1 = 1
        WLOG = LOG10(FW(IWS))  ! higher grid point
        WLOG1 = LOG10(Fw(IWS1)) ! lower grid point
        W1L = LOG10(Fwater)

c--------------------------Don't need this logic anymore (we think) because our FH2Os go down to 1e-06 at higher temperatures
!        ELSEIF ((MS.ge.6).and.(IWS.le.3))THEN
!        IWS = 9
!        IWS1 = 8  ! If Temp > 400 and MR is less than 1e-1 then the FH2O boundaries are at 1e-1 and 1, respectively.
!        WLOG = LOG10(FW(IWS))
!        WLOG1 = LOG10(FW(IWS1))
!        W1L = LOG10(Fwater)
c-------------------------------------


        ELSEIF(IWS.ge.9)THEN
        IWS = 9
        IWS1 = 8
        WLOG = LOG10(FW(IWS))
        WLOG1 = LOG10(FW(IWS1))
        W1L = LOG10(Fwater)
        ELSE  ! For IWS between 1 and 8
        IWS1 = IWS - 1
        WLOG = LOG10(FW(IWS))
        WLOG1 = LOG10(FW(IWS1))
        W1L = LOG10(Fwater)
        ENDIF



        IF ((ICS.eq.1).and.(FCO2.lt.1.e-04))THEN
        ICS1 = 1
        ICS = 2
        CCLOG = LOG10(FC(ICS))  ! higher grid point
        CCLOG1 = LOG10(FC(ICS1)) ! lower grid point
        CC1L = CCLOG1


        ELSEIF(ICS.eq.1)THEN
        ICS1 = 1
        ICS = 2
        CCLOG = LOG10(FC(ICS))  ! higher grid point
        CCLOG1 = LOG10(FC(ICS1)) ! lower grid point
        CC1L = LOG10(FCO2)

        ELSEIF(ICS.ge.5)THEN
        ICS = 5
        ICS1 = 4
        CCLOG = LOG10(FC(ICS))  ! higher grid point
        CCLOG1 = LOG10(FC(ICS1)) ! lower grid point
        CC1L = LOG10(FCO2)

!-------------------------
        ELSEIF(MS.ge.6)THEN   ! If Temp above 300 K, for inner edge FCO2 ALWAYS 3.3E-4. So FCO2 boundaries are now at 1e-4 and 3.3E-4
        ICS1 = 1  ! FCO2 lower boundary = 1.e-04
        ICS = 6   ! FCO2 upper boundary = 3.3e-04
        CCLOG = LOG10(FC(ICS))  ! higher grid point
        CCLOG1 = LOG10(FC(ICS1)) ! lower grid point
        CC1L = LOG10(FCO2)
!---------------------------------------

        ELSE  ! in between indices 2 and 4
        ICS1 = ICS - 1
        CCLOG = LOG10(FC(ICS))
        CCLOG1 = LOG10(FC(ICS1))
        CC1L = LOG10(FCO2)
        ENDIF

        DEN = (PLOG - PLOG1)*(tempg(MS) - TempL)*(WLOG - WLOG1)
     &  *(CCLOG - CCLOG1)  ! denominator

      IF (TempR.ge.tempg(6))DEN = (PLOG - PLOG1)*(tempg(MS) - TempL) ! If temp >= 400 then don't interpolate over CO2.
     & *(WLOG - WLOG1)



            ! initialize kmatrix_sol
           DO lam = 1,38
              DO IK1 = 1,16
                 kmatrix_sol(lam,IK1) = 1.e-60
              ENDDO
           ENDDO

      DO lam=1,38  ! interval loop

               

               

!        print 2223,PLOG,PLOG1,tempg(MS),TempL,WLOG,WLOG1,
!     &             CCLOG,CCLOG1,TempR,p1,Fwater,FCO2
!        pause

2223    FORMAT(1P17E14.5)
!FX=(TempR-TempL)/(tempg(MS)-TempL)
!FY=(P1L-PLOG1)/(PLOG-PLOG1)
!FZ = (W1L - WLOG1)/(WLOG - WLOG1)
!FA = (C1L - CLOG1)/(CCLOG - CCLOG1)


        IF (TempR.lt.tempg(6))THEN
	    DO IK1=1,16 !IK is a kcoefficient index.1 to 16 because 16 kcoefficient sums.
	    AK(1)=LOG10(kappa_sol(lam,MS1,IWS1, ICS, LS1,IK1))
	    AK(2)=LOG10(kappa_sol(lam,MS1, IWS1, ICS, LS,IK1))
	    AK(3)=LOG10(kappa_sol(lam,MS, IWS1, ICS, LS, IK1))
	    AK(4)=LOG10(kappa_sol(lam, MS, IWS, ICS, LS, IK1))
            AK(5)= LOG10(kappa_sol(lam, MS, IWS, ICS, LS1, IK1))
            AK(6)=LOG10(kappa_sol(lam, MS1, IWS, ICS, LS1, IK1))
            AK(7)=LOG10(kappa_sol(lam, MS1, IWS, ICS1, LS1, IK1))
            AK(8)=LOG10(kappa_sol(lam, MS1, IWS, ICS1, LS, IK1))
            AK(9)=LOG10(kappa_sol(lam, MS, IWS, ICS1, LS, IK1))
            AK(10)=LOG10(kappa_sol(lam, MS, IWS1, ICS1, LS, IK1))
            AK(11)=LOG10(kappa_sol(lam, MS1, IWS1, ICS1, LS, IK1))
            AK(12)=LOG10(kappa_sol(lam, MS1, IWS1, ICS1, LS1, IK1))
            AK(13)=LOG10(kappa_sol(lam, MS, IWS1, ICS1, LS1, IK1))
            AK(14)=LOG10(kappa_sol(lam, MS, IWS, ICS1, LS1, IK1))
            AK(15)=LOG10(kappa_sol(lam, MS, IWS1, ICS, LS1, IK1))
            AK(16)=LOG10(kappa_sol(lam, MS1, IWS, ICS, LS, IK1))

	    AN(1) = (PLOG - P1L)*(tempg(MS) - TempR)*(WLOG - W1L)
     &      *(CC1L - CCLOG1)
            AN(2) = (P1L - PLOG1)*(tempg(MS) - TempR)*(WLOG - W1L)
     &      *(CC1L - CCLOG1)
            AN(3) = (P1L - PLOG1)*(TempR - TempL)*(WLOG - W1L)
     &      *(CC1L - CCLOG1)
            AN(4) = (P1L - PLOG1)*(TempR - TempL)*(W1L - WLOG1)
     &      *(CC1L - CCLOG1)
            AN(5) = (PLOG - P1L)*(TempR - TempL)*(W1L - WLOG1)
     &      *(CC1L - CCLOG1)
            AN(6) =  (PLOG - P1L)*(tempg(MS) - TempR)*(W1L - WLOG1)
     &      *(CC1L - CCLOG1)
            AN(7) = (PLOG - P1L)*(tempg(MS) - TempR)*(W1L - WLOG1)
     &      *(CCLOG - CC1L)
            AN(8) = (P1L - PLOG1)*(tempg(MS) - TempR)*(W1L - WLOG1)
     &      *(CCLOG - CC1L)
            AN(9) = (P1L - PLOG1)*(TempR - TempL)*(W1L - WLOG1)
     &      *(CCLOG - CC1L)
            AN(10) = (P1L - PLOG1)*(TempR - TempL)*(WLOG - W1L)
     &      *(CCLOG - CC1L)
            AN(11) = (P1L - PLOG1)*(tempg(MS) - TempR)*(WLOG - W1L)
     &      *(CCLOG - CC1L)
            AN(12) = (PLOG - P1L)*(tempg(MS) - TempR)*(WLOG - W1L)
     &      *(CCLOG - CC1L)
            AN(13) = (PLOG - P1L)*(TempR - TempL)*(WLOG - W1L)
     &      *(CCLOG - CC1L)
            AN(14) = (PLOG - P1L)*(TempR - TempL)*(W1L - WLOG1)
     &      *(CCLOG - CC1L)
            AN(15) = (PLOG - P1L)*(TempR - TempL)*(WLOG - W1L)
     &      *(CC1L - CCLOG1)
            AN(16)= (P1L - PLOG1)*(tempg(MS) - TempR)*(W1L - WLOG1)
     &      *(CC1L - CCLOG1)

            AKL= 0.
            DO II = 1, 16
	    AKL= AKL + (AK(II)*AN(II)/DEN)
            ENDDO
       

          
	    kmatrix_sol(lam,IK1)=10**AKL  ! And then this needs to be interpolated over altitude (j)
!            if(lam .eq. 10)then
!            print 2223,kmatrix_sol(lam,IK1),real(lam),Temp,Fwater,FCO2
!     &            ,p1,real(IK1)
!            pause
!            endif
              if(AKL >100.0)then
                 print *, 'in interpsolar'                
!                 pause
              endif
  
	   	ENDDO ! ends IK1 loop



          ELSE   ! If TempR is 350 or greater
            	     DO IK1=1,16 !IK is a kcoefficient index.1 to 16 because 16 kcoefficient sums.
  	    AK(1) = LOG10(kappa_sol(lam,MS1,IWS,6,LS1, IK1))  ! ICS is 6 or FC(6) = 3.3E-04
            AK(2) = LOG10(kappa_sol(lam,MS,IWS,6,LS1, IK1))
            AK(3) = LOG10(kappa_sol(lam,MS1,IWS,6,LS, IK1))
            AK(4) = LOG10(kappa_sol(lam,MS,IWS,6,LS, IK1))
            AK(5) = LOG10(kappa_sol(lam,MS1,IWS1,6,LS1, IK1))
            AK(6) = LOG10(kappa_sol(lam,MS,IWS1,6,LS1, IK1))
            AK(7) = LOG10(kappa_sol(lam,MS1,IWS1,6,LS, IK1))
            AK(8) = LOG10(kappa_sol(lam,MS,IWS1,6,LS, IK1))
            AK(9)=0.
            AK(10)=0.
            AK(11)=0.
            AK(12)=0.
            AK(13)=0.
            AK(14)=0.
            AK(15)=0.
            AK(16) =0.

            AN(1) = (PLOG -P1L)*(tempg(MS) - tempR)*(W1L - WLOG1)
            AN(2) = (PLOG - P1L)*(TEMPR - TEMPL)*(W1L - WLOG1)
            AN(3) = (P1L - PLOG1)*(TEMPG(ms)- TEMPR)*(W1L - WLOG1)
            AN(4) = (P1L - PLOG1)*(TEMPR - TEMPL)*(W1L - WLOG1)
            AN(5) = (PLOG - P1L)*(TEMPG(MS) - TEMPR)*(WLOG - W1L)
            AN(6) = (PLOG - P1L)*(TEMPR- TEMPL)*(WLOG - W1L)
            AN(7) = (P1L - PLOG1)*(TEMPG(MS) - TEMPR)*(WLOG - W1L)
            AN(8) = (P1L - PLOG1)*(TEMPR - TEMPL)*(WLOG - W1L)
            AN(9) = 0.
            AN(10)=0.
            AN(11)=0.
            AN(12)=0.
            AN(13)=0.
            AN(14)=0.
            AN(15)=0.
            AN(16) =0.

            AKL= 0.
            DO II = 1, 16
	    AKL= AKL + (AK(II)*AN(II)/DEN)
            ENDDO ! Loop to calculate AKL

            kmatrix_sol(lam,IK1)=10**AKL  ! And then this needs to be interpolated over altitude (j)
            	  	ENDDO ! ends IK1 loop

            
  	    ENDIF

            ENDDO ! ends interval loop
       RETURN

       END ! Ends subroutine


