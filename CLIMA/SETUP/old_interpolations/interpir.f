
c     SUBROUTINE INTERP(Temp1,p1,xkappa,kappa) 
      SUBROUTINE INTERPIR(Temp1,p1,kappa,j,Fwater,kmatrix_ir) ! Written by Ramses Ramirez and Ravi Kopparapu


c-mm  This subroutine will READ in the exponential sum datafiles in the
c-mm  ir CH4 as well as the mixed CO2/H2O coefficients and output the 
!     appropriate value for each wavelength interval.  
c-mm  Indices for xkappa and kappa are as follows
c-mm  Index #1:  Temp   1=less than T   2=more than T
c-mm  Index #2:  Pres   1=less than p   2=more than p
c-mm  Index #3:  lam    wavelength bin
c-mm  Index #4:  i      gauss point

! RR AND RK removed the species index because xkappa, kappapp, kappapm, and kappa arrays only have CH4 now. 
c-rr    i is wavelength counter and j is the altitude counter used in solar, respectively.
C new common block, von Paris, 21/04/2006
	include '../header.inc'
	PARAMETER(NSOL=38,IK=16, NGS=5, NS=3, NS1 = NS+1 )
        INTEGER LINE,i,j,lam,Tindex,pindex,Tp1,pp1,Temp
        REAL kappapp(55,8),kappapm(55,8) ! methane kappas
        REAL kappa(55,6), AN(16), AK(16) ! methane kappas
        REAL TempR,tempg(8),press(8), FC(6), FW(9),p1,kappa_sol,
     &  Fwater,kmatrix_IR(55,IK), kappa_ir, p, Temp1
	DATA tempg/100., 150., 200., 250., 300., 350., 400., 600./ ! 8 temps 
	DATA press/.00001,.0001,.001,.01, 0.1, 1., 10., 100./  ! 8 pressures
        DATA FC/1.e-04, 1.e-03, 1.e-02, 1.e-01, 1., 3.3e-04/  ! First 5 FCO2s for interpolation for T less than 400K. For T greater than 400 K, FC = 3.3e-04
        DATA FW/1.e-08, 1.e-07, 1.e-06, 1.e-05, 1.e-04, 1.e-03, 1.e-02,
     &   1.e-01, 1./ ! 9 FH2Os
	COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL), 
     &  kappa_sol(NSOL,8,9,6,8,16) ! Added new kappa matrix for mixed CO2-H2O  
c-rr !3/23/11 put CIA matrix in IRDATA. !Added O3 absorption coefficient vector (plus c and d vectors), first 14 terms. 3/19/2012 ! 38 layers, 7 temps, 9 FH2Os, 6 FCO2s, 9 pressures, 16 coefficients


        COMMON/IRDATA/WEIGHTCH4(6),xkappa(3,12,55,8), ! xkappa is pre-loaded methane absorption coefficient array 
     &   CIA(7,55), CPRW(ND,55)  !c-rr !3/23/11 put CIA matrix in IRDATA. Added weight co2_h2O 3/20/2012

      COMMON/VARIR/kappa_ir(55, 8, 9, 6, 8, IK)! Added kappa matrix in IR for kpsectrum mixed CO2/H2O coefficients 3/20/2012  
      COMMON/weightsIR/ weightco2_h2oIR(IK)



c-rr    Created new common block solar data and put in readsol.f and solar.f as well
!       COMMON/SOLARDATA/weightco2_h2oSOL(16), weights(2,NSOL,IK)
!     & ,kmatrix_sol(NSOL,IK)  ! c-rr Added weightco2_h2O data array and weights array that combines co2_h2o and Ch4 weight arrays 3/19/2012  ! Added in 16 term weight array
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND),
     &   FH22,FNC
 

      
        
!       do i = 1,16
!         print *, 'weight=',weightco2_h2oIR(i)
!       enddo
!       read(*,*)
 
!  	Assigning TempR to temp (due to variable name overuse)
	TempR=Temp1

        DO nlam = 1,55
              DO IK1 = 1,16
                 kmatrix_ir(nlam,IK1) = 1.e-60
              ENDDO
           ENDDO

        

          

          DO I = 1,16   ! Initializing AK and AN
           AK(I) = 1.e-60
           AN(I) = 1.e-60
          ENDDO
       
      DO lam=1,55 ! interval loop


	   Do L=1,8  ! For 8 pressures
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
	      MS=8
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
        P1L  = LOG10(p1)     ! Added by Ravi 3/29/2012
        ENDIF

        	
        
        IF(LS.ge.8) THEN
        LS = 8		! P1 is at or above the highest grid pressure
	LS1 = 7
        PLOG = LOG10(press(LS))
        PLOG1 = LOG10(press(LS1))
        P1L=LOG10(p1)
	ENDIF

C----------------------------------------

	IF(MS.eq.1)THEN
  	MS1=MS
        MS = 2		!TempR is at or below the lowest grid Temperature
	TempR = tempg(MS)
        TempL = tempg(MS1)
        ! Change all tempg(MS1) to TempL

	ELSEIf(MS.ge.8) THEN
	MS=8		!TempR is at or above the highest grid Temperature
	MS1=MS-1		
        TempR = tempg(MS)
        TempL = tempg(MS1)
        ELSE
	MS1=MS-1  !This calculates MS1 for values within a temperature range 
        TempL = tempg(MS1)  ! high grid point is already tempg(MS) which is done in DEN. TempR is just input Temp
	ENDIF


!------------------------------------------------

!       Also, For intervals between 15 and 20, zero out our water kcoefficients and replace with 8-12 water continuum in ir.f 4/6/2012

	IF((lam.gt.14).and.(lam.lt.21))then
        IWS=2
        IWS1= 1
        WLOG = LOG10(FW(IWS))
        WLOG1 = LOG10(FW(IWS1))
        WIL = WLOG1
        ENDIF

        IF ((IWS.eq.1).and.(Fwater.lt.1.e-08))THEN
        IWS=2
        IWS1 = 1
        WLOG = LOG10(FW(IWS))  ! higher grid point
        WLOG1 = LOG10(FW(IWS1)) ! lower grid point
        W1L = WLOG1
     
!        print *,'p1 is ....',p1,P1L

	

        ELSEIF(IWS.eq.1) THEN
        IWS = 2
        IWS1 = 1
        WLOG = LOG10(FW(IWS))  ! higher grid point
        WLOG1 = LOG10(Fw(IWS1)) ! lower grid point
        W1L = LOG10(Fwater)

c--------------------------Don't need this logic anymore (we think) because our FH2Os go down to 1e-06 at higher temperatures
!        ELSEIF ((MS.ge.6).and.(IWS.le.8))THEN
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
        ICS1=1
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
        ELSEIF(MS.ge.6)THEN   ! If Temp between 300 and 400, for inner edge FCO2 ALWAYS 3.3E-4. So FCO2 boundaries are now at 1e-4 and 3.3E-4
        ICS1 = 1 ! FCO2 lower boundary = 1.e-04
        ICS = 6 ! FCO2 upper boundary = 3.3e-04
        CCLOG = LOG10(FC(ICS))  ! higher grid point
        CCLOG1 = LOG10(FC(ICS1)) ! lower grid point
        CC1L = LOG10(FCO2)
!---------------------------------------

        ELSE
        ICS1 = ICS - 1  ! in between indices 2 and 4
        CCLOG = LOG10(FC(ICS))
        CCLOG1 = LOG10(FC(ICS1))
        CC1L = LOG10(FCO2)
        ENDIF

        

        DEN = (PLOG - PLOG1)*(tempg(MS) - TempL)*(WLOG - WLOG1)
     &  *(CCLOG - CCLOG1)  ! denominator!

       IF (TempR.ge.tempg(6))DEN = (PLOG - PLOG1)*(tempg(MS) - TempL) ! If temp >= 400 then don't interpolate over CO2.
     & *(WLOG - WLOG1)
      

     
              

            ! initialize kmatrix_ir
!           DO nlam = 1,55
!              DO IK1 = 1,16
!                 kmatrix_ir(nlam,IK1) = 1.e-60
!              ENDDO
!           ENDDO

        

          

!          DO I = 1,16   ! Initializing AK and AN
!           AK(I) = 1.e-60
!           AN(I) = 1.e-60
!          ENDDO


!      DO lam=1,55 ! interval loop

!        if(NST==5)then
!         if ((lam.eq.15).and.(J.eq.69))then
!        print 2222,PLOG,PLOG1,tempg(MS),TempL,TempR,WLOG,
!     &             WLOG1,
!     &             CCLOG,CCLOG1,Temp1,p1,Fwater,FCO2,DEN,
!     &             real(j),
!     &             real(lam)
!        pause
        

!        print *, P1L, CC1L, W1L
!        endif

!FX=(TempR-TempL)/(tempg(MS)-TempL)
!FY=(P1L-PLOG1)/(PLOG-PLOG1)
!FZ = (W1L - WLOG1)/(WLOG - WLOG1)
!FA = (C1L - CLOG1)/(CCLOG - CCLOG1)

          
           IF (TempR.lt.tempg(6))THEN

	     DO IK1=1,16 !IK is a kcoefficient index.1 to 16 because 16 kcoefficient sums.
	    AK(1)=LOG10(kappa_ir(lam,MS1,IWS1, ICS, LS1,IK1))
           
	    AK(2)=LOG10(kappa_ir(lam,MS1, IWS1, ICS, LS,IK1))
	    AK(3)=LOG10(kappa_ir(lam,MS, IWS1, ICS, LS, IK1))
	    AK(4)=LOG10(kappa_ir(lam, MS, IWS, ICS, LS, IK1))
            AK(5)= LOG10(kappa_ir(lam, MS, IWS, ICS, LS1, IK1))
            AK(6)=LOG10(kappa_ir(lam, MS1, IWS, ICS, LS1, IK1))
            AK(7)=LOG10(kappa_ir(lam, MS1, IWS, ICS1, LS1, IK1))
            AK(8)=LOG10(kappa_ir(lam, MS1, IWS, ICS1, LS, IK1))
            AK(9)=LOG10(kappa_ir(lam, MS, IWS, ICS1, LS, IK1))
            AK(10)=LOG10(kappa_ir(lam, MS, IWS1, ICS1, LS, IK1))
            AK(11)=LOG10(kappa_ir(lam, MS1, IWS1, ICS1, LS, IK1))
            AK(12)=LOG10(kappa_ir(lam, MS1, IWS1, ICS1, LS1, IK1))
            AK(13)=LOG10(kappa_ir(lam, MS, IWS1, ICS1, LS1, IK1))
            AK(14)=LOG10(kappa_ir(lam, MS, IWS, ICS1, LS1, IK1))
            AK(15)=LOG10(kappa_ir(lam, MS, IWS1, ICS, LS1, IK1))
            AK(16)=LOG10(kappa_ir(lam, MS1, IWS, ICS, LS, IK1))

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
!             if ((lam.eq.15).and.(IK1.eq.16).and.(J.eq.1))then
!            print 2222, (AKL), AN(II),AK(II),DEN
!              endif
            
!            print 120,MS1, MS, LS1, LS,  
!     &                    ICS1, ICS, IWS1, IWS, tempg(MS),TempR,TempL,
!     &                    WLOG, W1L,WLOG1, P1L,PLOG, PLOG1, CCLOG,CC1L 
!     &                    ,CCLOG1,
!     &                    kmatrix_ir(lam,1)
!            pause
            ENDDO ! Loop to calculate AKL
           

              
   !          (kappa_ir(lam,MS1,IWS1, ICS, LS1,IK1))

            
            	    kmatrix_ir(lam,IK1)=10**AKL  ! And then this needs to be interpolated over altitude (j)
!                         if((NST.eq.5) .and. (j==7))then

!                        if((j .eq. 69).and.(lam.eq.15))then
!            	  	 print 2222, (kmatrix_ir(lam,II),II=1,16),tempR,p1,
!     &                              Fwater,FCO2
!                         print *, lam, MS1, IWS1, ICS, LS1, IK1
!                         pause
!                        endif



!               if ((lam.eq.15).and.(IK1.eq.16).and.(J.eq.1))then
!             print *, 'INTERPIR',kmatrix_ir(lam,IK1), Temp1, P1,Fwater,
!     &             FCO2
!             pause
!              print 2222, ((AK(II)), II=1,16), AKL, real(j),real(lam)
!             pause
!             endif
!                        
                     ENDDO ! ends IK1 loop




                
!                    print 120,MS1, MS, LS1, LS,  
!     &                    ICS1, ICS, IWS1, IWS, tempg(MS),TempR,TempL,
!     &                    WLOG, W1L,WLOG1, P1L,PLOG, PLOG1, CCLOG,CC1L 
!     &                    ,CCLOG1,
!     &                    kmatrix_ir(lam,1)
!                   pause
120 		format(8(2x,i3),1p17e14.5)

             

           ELSE   ! If TempR is 400 or greater
            	     DO IK1=1,16 !IK is a kcoefficient index.1 to 16 because 16 kcoefficient sums.
  	    AK(1) = LOG10(kappa_ir(lam,MS1,IWS,6,LS1, IK1))  ! ICS is 6 or FC(6) = 3.3E-04
            AK(2) = LOG10(kappa_ir(lam,MS,IWS,6,LS1, IK1))
            AK(3) = LOG10(kappa_ir(lam,MS1,IWS,6,LS, IK1))
            AK(4) = LOG10(kappa_ir(lam,MS,IWS,6,LS, IK1))
            AK(5) = LOG10(kappa_ir(lam,MS1,IWS1,6,LS1, IK1))
            AK(6) = LOG10(kappa_ir(lam,MS,IWS1,6,LS1, IK1))
            AK(7) = LOG10(kappa_ir(lam,MS1,IWS1,6,LS, IK1))
            AK(8) = LOG10(kappa_ir(lam,MS,IWS1,6,LS, IK1))
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
!             if((NST.eq.4)  )then

    
!           print 2222, ((AK(II)), II=1,16), real(j),real(lam)
!            pause
!            endif


            kmatrix_ir(lam,IK1)=10**AKL  ! And then this needs to be interpolated over altitude (j)
            	  	ENDDO ! ends IK1 loop



!             if(j .eq. 1)then
!            	  	 print 2222, (kmatrix_ir(lam,II),II=1,16),tempR
!                        endif
!            print 120, kmatrix_ir(lam,1), MS1, MS, LS1, LS, ICS1, ICS, IWS1, IWS
!            pause
!120 		format(1pe14.5, 8(2x,i3))


  	    ENDIF
       
             ENDDO ! ends interval loop             
!          if((NST.eq.5) .and. (j==7))then
!              pause
!          endif

!              print 2223, kappa_ir(lam,MS1,IWS1, ICS, LS1,IK1),
!     &                    AKL,kmatrix_ir(lam,IK1),Temp1,P1,
!     &                    lam,MS1,IWS1,
!     &                    ICS,LS1,IK1
!             print 2222,kmatrix_ir(lam,IK1),AKL,real(j),Temp1,p1,
!     &                  tempg(MS), TempL
!             print 2222, (AK(i), i=1,16), AKL,real(j)
!             PRINT 2222, (AN(i), i=1,16),DEN
              
!             print *,'------------------------------------'

!               if(AKL >100.0)then
!                 print *, 'in interpir'                
!                 pause
!              endif
!
!            
!             pause

2222    format(1p21e14.5)
         
!         pause


!2223    format(1p5e14.5,0p,6(2x,i3))

!--------------------------------------------------------------------

      do 9400 i=1,55
         do 9401 ij=1,8
               kappa(i,ij)=0.0
 9401    continue
 9400 continue
        Temp=int(Temp1)
        Tindex=Temp/50
        p = p1*1000
        pindex=int(log10(p*1.0e7))
        Tp1=Tindex+1
        pp1=pindex+1


c-mm  These two loops will linearly interpolate kappa values to solve
c-mm  for any T and p.  kappapm interpolates xkappa for the two values
c-mm  having a pressure lower than p.  kappapp interpolates xkappa for
c-mm  the two values having a pressure higher than p.  kappa then
c-mm  interpolates temperature.

 
        do 9607 lam=5,42  ! changed methane intervals from intervals 1-38 to 5-42 3/21/2012
	if (Temp.le.150) then  ! Also changed Tindices for xkappa from 3,6 and 8 to 1,2, and 3, respectively.
		do 9620 i=1,6
		kappapm(lam,i) = xkappa(1,pindex,lam,i)
		kappapp(lam,i) = xkappa(1,pp1,lam,i)
		kappa(lam,i) = (1.-(abs(log10(p)-int(log10(p)))))*
     2  kappapp(lam,i)+abs(log10(p)-int(log10(p)))*kappapm(lam,i)
 9620  continue
        elseif (Temp.le.300) then
                do 9609 i=1,6
                        kappapm(lam,i)=((mod(Temp,150))/150.)*
     2  xkappa(2,pindex,lam,i)+(1.-((mod(Temp,150))/150.))*
     3  xkappa(1,pindex,lam,i)
                        kappapp(lam,i)=((mod(Temp,150))/150.)*
     2  xkappa(2,pp1,lam,i)+(1.-((mod(Temp,150))/150.))*
     3  xkappa(1,pp1,lam,i)
                        kappa(lam,i)=(1.-(abs(log10(p)-
     2  int(log10(p)))))*kappapp(lam,i)+
     3  abs(log10(p)-int(log10(p)))*kappapm(lam,i)
 9609   continue
        else
                do 9608 i=1,6
                        kappapm(lam,i)=((mod(Temp,300))/300.)*
     2  xkappa(3,pindex,lam,i)+(1.-((mod(Temp,300))/300.))*
     3  xkappa(2,pindex,lam,i)
                        kappapp(lam,i)=((mod(Temp,300))/300.)*
     2  xkappa(3,pp1,lam,i)+(1.-((mod(Temp,300))/300.))*
     3  xkappa(2,pp1,lam,i)
                        kappa(lam,i)=(1-(abs(log10(p)-
     2  int(log10(p)))))*kappapp(lam,i)+
     3  abs(log10(p)-int(log10(p)))*kappapm(lam,i)
 9608           continue
        endif
 9607   continue

        RETURN
        END
