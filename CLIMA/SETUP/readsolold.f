      SUBROUTINE READSOL      
c-as This subroutine reads the incoming flux from the Sun and the
c-as parameters for the subroutine SOLAR and CONVEC

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NF=55)
      PARAMETER(NS=3, NS1=NS+1, NS4=NS+5)
      PARAMETER(NT=76, MT=36)
      PARAMETER(NSOL=38, NGS=5)
      PARAMETER(IK=8)
      REAL kappa, kmatrix, pres(6), temp(3), alch4
      DATA kappa/5472*0./ 
c-rr  The weights below here are those for CO2 in the near-IR with new coefficients 
c     from Richard R. Freedman's eb.txt 10/11/2010 
      DATA weight/1.6523105144e-1, 3.0976894856e-1, 3.0976894856e-1, 
     2 1.6523105144e-1, 8.6963711284e-3, 1.6303628872e-2, 
     3 1.6303628872e-2, 8.6963711284e-3/  
      DIMENSION ALPHAZ(4,2),BETAZ(4,2),NPROB(2),
     &  NG(2),SIGG(4,2,NSOL)

      COMMON/ABLOK/LTYPE(NF,3),XH2O(NF),YH2O(NF),XCO2(NF),YCO2(NF),
     2  AXH(NF),AYH(NF),BXH(NF),BYH(NF),AXC(NF),AYC(NF),BXC(NF),
     3  BYC(NF),PDOP(NF),CPR(NF),TPR(NF),PATH(NS4),PATHP(NS4),
     4  PATHT(NS4),P1,P2,T1,T2,TAU2,TAUP2,
     5  ALPHA(4),BETH2O(4,5,NF),
     6  BETCO2(4,5,NF),CA(19),CB(19),CC(19),CD(19),CE(19),CH(19),
     7  CK(19)
      COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1)
      COMMON/CH4BLOCK/ALPHACH4T188(4,17),BETACH4T188(4,17),
     & ALPHACH4T295(4,17),BETACH4T295(4,17),ALPHACH4Kark(4,21),
     & BETACH4Kark(4,21),GAMMAEXP188(17),GAMMAEXP295(17),
     & ALPHACH4NEW(6),BETACH4NEW(17,3,5,6),ALCH4(6,38)
      COMMON/FBLOK/TTAB(NT),PVAP(NT),DPVAP(NT),SVC(NT),DSV(NT),DSC(NT)
     2  ,RHOV(NT),DRHOV(NT),BETAM(70,75),TCP(75),PCP(70),DPDTL(70,75),
     3  DRDTL(70,75)
      COMMON/GBLOK/TCTAB(MT),PCVAP(MT),BETASC(MT),DPCVAP(MT),
     &  DRCVAP(MT),SVSC(MT),DSCC(MT),TKTAB(MT),TCC(25),PCC(36),
     &  BETAMC(25,36),CPC(25,36),DVDTC(25,36),DVDPC(25,36),
     &  DSVC(MT)
      COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL),kappa(NSOL,3,6,8)
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND)
c-rr Created new solardata common block to hold richard freedman's near-IR CO2 coefficients
      COMMON/SOLARDATA/kmatrix(NSOL,ND,IK), weight(8)

      DATA ALPHACH4NEW/0.08566225,0.1803808,0.23395695,0.23395695,
     1 0.1803808,0.08566225/

C READ THE STEAM TABLE
      do i=1,4
       READ(3,*)
      enddo
      READ(3,300) (TTAB(I),PVAP(I),DPVAP(I),SVC(I),DSV(I),DSC(I),
     2  RHOV(I),DRHOV(I),I=1,NT)
 300  FORMAT(1P8E12.5)
      DSV(NT) = 0.

C   READ THE HIGH TEMPERATURE (UNSATURATED) STEAM TABLE
      do i=1,2
       READ(3,*)
      enddo
      DO 58 I=1,7
      M1 = 10*(I-1) + 1
      M2 = M1 + 9
      READ(3,302) (PCP(M),M=M1,M2)
 302  FORMAT(5X,10F11.0)
C
      DO 59 N=1,75
      READ(3,303) TCP(N),(DPDTL(M,N),M=M1,M2)
      READ(3,304) (DRDTL(M,N),M=M1,M2)
      READ(3,304) (BETAM(M,N),M=M1,M2)
  59  READ(3,304)
  58  CONTINUE 
 303  FORMAT(1X,F6.0,1P10E11.4)
 304  FORMAT(7X,1P10E11.4)

c    READ THE SATURATED CO2 TABLE
      do ii=1,5
       READ(9,*)
      enddo
      READ(9,306) (TCTAB(I),PCVAP(I),BETASC(I),DPCVAP(I),DRCVAP(I),
     2  SVSC(I),DSVC(I),DSCC(I),I=1,MT)
 306  FORMAT(1X,F7.2,7E10.3)
C
c  READ THE UNSATURATED CO2 TABLE
      READ(9,177)
      DO 76 L=1,6
      IS = 1 + 6*(L-1)
      IF = IS + 5
      READ(9,307) (PCC(I),I=IS,IF)
 307  FORMAT(///10X,6(F4.1,7X))
      DO 68 K=1,25
      READ(9,308) TCC(K),(BETAMC(K,I),I=IS,IF)
 308  FORMAT(/F5.1,2X,6E11.4)
      READ(9,309) (CPC(K,I),I=IS,IF)
 309  FORMAT(7X,6E11.4)
      READ(9,309) (DVDTC(K,I),I=IS,IF)
  68  READ(9,309) (DVDPC(K,I),I=IS,IF)
  76  CONTINUE

***INPUT TO SOLAR-TWO-STREAM*** (DMW)
C
C   TAUAER = OPTICAl DEPTH OF AEROSOL (IF NO AEROSOL, TAUAER = 0)
C   SIGERT = EXTINCTION COEFFICIENT FOR AEROSOL (IF NO AEROSOL, SIGERT = 0)
C   FMA = "F" (SECOND MOMENT/5) FOR AEROSOL (IF NO AEROSOL, FMA = 0)
C
      DO 128 I=1,NSOL
         TAUAER(I) = 0.
         SIGERT(I) = 0.
         FMA(I) = 0.
 128  CONTINUE

C
C  READING of new data
C  I corresponds to interval 22-38
C  J corresponds to Temps 112,188,295 [Kelvin]
C  K corresponds to pressures 0.0001,0.001,0.01,0.1,1.0 [Bars]
	
	DO I=1,17
	READ(21,*)
	READ(21,*)
	DO J=1,3
	READ(21,*)
	READ(21,*)
	DO K=1,5
	READ(21,902) (BETACH4NEW(I,J,K,L),L = 1,6)
	END DO
	END DO
	END DO
	
902   FORMAT(15X,1PE11.5,2X,E11.5,2X,E11.5,2X,
     1 E11.5,2X,E11.5,2X,E11.5)

C
C  READING of Kathy Rages data (near IR CH4 exponential sums)
      do i=1,4
       READ(8,*)
      enddo
      DO 173 I = 1,17
       READ(8,182) GAMMAEXP188(I)
       READ(8,176) (ALPHACH4T188(K,I),K = 1,4)
       READ(8,179) (BETACH4T188(K,I),K = 1,4)
 173  CONTINUE
      do i=1,4
       READ(8,*)
      enddo 
      DO 174 I = 1,17
       READ(8,182) GAMMAEXP295(I)
       READ(8,176) (ALPHACH4T295(K,I),K = 1,4)
       READ(8,179) (BETACH4T295(K,I),K = 1,4)
 174  CONTINUE
     
c Reading Karkoshka data for CH4
      do i=1,5
        READ(8,*)
      enddo
      DO 172 I = 1,21
       READ(8,181) (ALPHACH4Kark(K,I),K = 1,4)
       READ(8,181) (BETACH4Kark(K,I),K = 1,4)
       READ(8,177)
 172  CONTINUE
     
 176  FORMAT(4X,F6.4,4x,F6.4,4x,F6.4,4x,F6.4)
 177  FORMAT(/)
 178  FORMAT(//)
 179  FORMAT(F9.5,2x,F9.4,1x,F9.3,1x,F12.2)
 181  FORMAT(2x,F8.5,2x,F8.5,2x,F8.5,2x,F8.5)
 182  FORMAT(4x,F8.5)

C   READ EXPONENTIAL SUM DATAFILES

**** H2O parameters for exponential sums
      DO 56 I=1,30
       READ(8,501)
  56   READ(8,500) ((BETH2O(K,L,I),K=1,4),L=1,4)
      DO 156 I=31,NF
       READ(8,501)
 156   READ(8,500) ((BETH2O(K,L,I),K=1,4),L=1,5)
   
C ***** TEMPORARY FILL FOR H2O EXP SUMS AT 10 BARS *****
      DO 57 K=1,4
      DO 57 I=1,30
  57  BETH2O(K,5,I) = BETH2O(K,4,I)

**** CO2 parameters for exponential sums
      READ(8,*)
      READ(8,*)
      DO 54 I=9,23
      READ(8,501)
  54  READ(8,500) ((BETCO2(K,L,I),K=1,4),L=1,5)
      DO 55 I=27,45
      READ(8,501)
  55  READ(8,500) ((BETCO2(K,L,I),K=1,4),L=1,5)
      READ(8,501)
      READ(8,500) ((BETCO2(K,L,48),K=1,4),L=1,5)

 500  FORMAT(12X,E14.8,5X,E14.8,5X,E14.8,5X,E14.8)
 501  FORMAT(////)
C
C   CONVERT TO CM2/GM
      DO 34 K=1,4
      DO 34 L=1,5
      DO 34 I=1,NF
      BETH2O(K,L,I) = BETH2O(K,L,I) * 6.023E23/18.
  34  BETCO2(K,L,I) = BETCO2(K,L,I) * 6.023E23/44.
C
C   FILL UP SOLAR ABSORPTION MATRICES
      DO 35 K=1,4
      DO 35 L=1,5
      DO 35 I=14,NSOL
      J = 69 - I
      BETIR1(K,L,I) = BETH2O(K,L,J)
  35  BETIR2(K,L,I) = BETCO2(K,L,J)
C
C   Read Solar Data - formerly Ackerman's routine SOLAR (7-95 DMW)
C   First find FCO2.
C
      DO 435 I=1,NSOL
      READ(4,*) IJUNK
      READ(4,*) ALAMBDA(I),(NPROB(L),L=1,2),SOLINT(I),(NG(L),L=1,2)
         DO 436 K=1,4
            READ(4,571) (ALPHAZ(K,L),L=1,2),(BETAZ(K,L),L=1,2)
 571        FORMAT(1X,2(F10.8,2X),2E14.7)
 436     CONTINUE
         DO 437 L=1,2
            NPR(L,I) = NPROB(L)
            NGAS(L,I) = NG(L)
            DO 438 K=1,4    
               SIGG(K,L,I) = BETAZ(K,L)
               WGHT(K,L,I) = ALPHAZ(K,L)
 438        CONTINUE
 437     CONTINUE  
 435  CONTINUE
      IF (FCO2 .LT. 0.1) THEN
         NPR(2,15) = 1
         NGAS(2,15) = 4
         WGHT(1,2,15) = 1.
         WGHT(2,2,15) = 0.
         SIGG(1,2,15) = 0.02
      END IF
      DO 440 K=1,4  !sums
         DO 441 L=1,5 !gas number
            DO 442 I=1,13 !wavelength
               BETIR1(K,L,I)=SIGG(K,1,I)
	       
 442        CONTINUE
            DO 443 I=1,20
               BETIR2(K,L,I)=SIGG(K,2,I)
	       
 443        CONTINUE
 441     CONTINUE
 440  CONTINUE




      
      
C-rr Temporarily replacing Ackerman's solar fluxes with Gliese581 solar fluxes 11/05/2010
c        read(10,*) ! skips one line before reading 
c       DO I =1,38
c	  read(10,100)SOLINT(I)
c	  print 100, SOLINT(I)
c        ENDDO
c 100  format(6x,1pe11.4)
       
c      DO I=1,38
c	print *, SOLINT(I)
c      ENDDO
      
C   RETURN
C        
c	DO I=1,13
c		DO K=1,4
c		print 500, (BETIR1(k,l,i), l=1,5)
c		ENDDO
c	ENDDO
C        
C-rr Test to print out BETIR1
C		print *, "BETIR1 INTERVALS"
C	DO 800 I=1,13
C		print *, "Interval=", I
C		print 700
C		DO 801 K=1,4
C		print 600, (BETIR1(K,L,I), L=1,5)
C 801		CONTINUE
C 800	CONTINUE
C-rr Test to print out BETIR2
C		print *, "BETIR2 INTERVALS"
C	DO 802 I=1,20
C		print *, "Interval=", I
C		print 700
C		DO 803 K=1,4
C		print 600, (BETIR2(K,L,I), L=1,5)
C 803		CONTINUE
C 802	CONTINUE
	
C 600    format(5(1pe11.4,1x))
C 700 	format(4x,'ksums',8x,'copy',8x,'copy',8x,'copy',8x,'copy')
     	

C-rr It turns out that above BETIR1 matrix is just the kcoefficients (betas) for gas 1(H20) and at least Intervals 1-13   C    are read from Tom Ackerman's solar38.pdat.
C    BETIR2 is just the betas for gas 2. Intervals 21-38 of BETIR2 are Kasting's CO2 kcoefficients from nearIRexpsums.pdat
C    (minus intervals 22 and 23 where there is no CO2 absorption) whereas Intervals 1-20 pertain to Ackerman's 
C   (solar_data_38.pdat). 
C    10/27/10

C-rr READ in Richard Freedman's CO2 KAPPA matrix for solar intervals 17-38 10/11/2010
        read(19,*)
	      do i=17,NSOL   
	        read(19,*) 
             do it=1,3 
	           read(19,*) 
			    do ip=1,6   
			      read(19,301),(kappa(i,it,ip,k), k=1,8) 
		      do k=1,8
                kappa(i,it,ip,k)=kappa(i,it,ip,k)
 
	          enddo
              do k=1,8
                 kappa(i,it,ip,k) = amax1(kappa(i,it,ip,k),1.e-60)
              enddo

			   enddo
		     enddo
         enddo
         
      i=38
      it=3
c      do ip=1,1   
c	    print 301,(kappa(i,it,ip,k), k=1,8) 
c      enddo

c      print *, 'BETCO2=',BETCO2(1,1,48)
c      print *, 'BETCO2=', BETIR2(1,1,31)
      
 301    format(7x,1pe11.4,2x,1pe11.4,6(1x, 1pe11.4))
 
c-rr COMBINING ALPHACH4Kark and ALPHACH4NEW into one ALCH4 matrix 10/23/2010
      DO I = 1,21
         DO K=1,4
         ALCH4(K,I) = ALPHACH4Kark(K,I)
         ENDDO   
      ENDDO
      DO I=22,38
	  DO K=1,6
          ALCH4(K,I) = ALPHACH4NEW(K)
          ENDDO
      ENDDO

c      DO I=1,38
c      print 600, (ALCH4(K,I), K=1,6)
c      ENDDO
c 600  format(f6.4,3x,5(f6.4,3x))
     
      END

