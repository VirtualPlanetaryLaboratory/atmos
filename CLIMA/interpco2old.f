	SUBROUTINE interpco2(Temp,p1,i,j)
c-rr    i is wavelength counter and j is the altitude counter used in solar, respectively.

	INCLUDE '../header.inc'
	PARAMETER(NSOL=38,IK=8)
	REAL TempR,tempg(3),press(6),p1,kappa,kmatrix,Temp
	INTEGER i,j
	DATA tempg/150., 225., 300./ 
	DATA press/.0001, .001, .01, 0.1, 1., 5./ 
	COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL), kappa(NSOL,3,6,8) 
c-rr    Created new common block solar data and put in readsol.f and solar.f as well
	COMMON/SOLARDATA/kmatrix(NSOL,ND,IK),weight(8)

c  	Assigning TempR to temp (due to variable name overuse)
	TempR=Temp


C	There are separate loops for pressures lower than 1e-4 bar,  
C	for pressures higher than 1 bar, and those in between. For 
C	the pressures in between:

	Do INTVS=17,38  ! THE MASTER LOOP FOR THE 22 intervals.
	   Do L=1,6
		LS=L 
		 
	     If (p1.le.press(L)) exit!If an individial p greater than curent p set Ls=6 
	     LS = 7
	   enddo
	   

              Do M=1,3  ! For three temperature intervals (150,225,300K)
		 MS=M  
		If (TempR.lt.tempg(M)) exit  
		MS=4 
	      enddo
	      


c       For pressures within the grid
	IF(LS.eq.1) GO TO 10
	IF(LS.ge.6) GO TO 11
	LS1 = LS-1
	PLOG=ALOG10(press(LS))
	PLOG1=ALOG10(press(LS1))
	P1L=ALOG10(p1)
	FY=(P1L-PLOG1)/(PLOG-PLOG1)
	GO TO 20
C
c	For pressures outside of grid:

  10	LS1 = LS	! P1 is at or below the lowest grid pressure
	FY = 1.
	GO TO 20
  11	LS = 6		! P1 is at or above the highest grid pressure
	LS1 = 6
	FY = 1.
  20	CONTINUE
C

	IF(MS.le.1) GO TO 12
	If(MS.ge.4) GO TO 13
	MS1=MS-1  !This calculates MS1 for values within a temperature range (MS between 2 and 3)
	FX=(TempR-tempg(MS1))/(tempg(MS)-tempg(MS1)) ! This calculates FX for values within a temperature range
	GO TO 33
c
c	For temperatures outside of the grid:

  13	MS=MS-1		!TempR is at or above the highest grid Temperature
	MS1=MS		
	FX = 1.
	GO TO 33
  12	MS1=MS		!TempR is at or below the lowest grid Temperature
	FX = 1.
	GO TO 33
  33	CONTINUE



 	  DO IK1=1,8  !IK is a kcoefficient index.1 to 8 because 8 kcoefficient sums.
	    AK1L=LOG10(kappa(INTVS,MS1,LS1,IK1))
	    AK2L=LOG10(kappa(INTVS,MS,LS1,IK1))
	    AK3L=LOG10(kappa(INTVS,MS,LS,IK1))
	    AK4L=LOG10(kappa(INTVS,MS1,LS,IK1))
	    AKL=FX*FY*AK3L+FX*(1.-FY)*AK2L+(1.-FX)*FY*AK4L+(1.-FX)*
     2	        (1.-FY)*AK1L
 200   FORMAT(1X,1P5E10.3)
 
          
	    kmatrix(INTVS,j,IK1)=10**AKL  ! And then this needs to be interpolated over altitude (j)

  
	  ENDDO
c            if ((J.eq.80).and.(INTVS.eq.17).and.(I.eq.20))then
c		  print *, 'T=',TempR
c		  print *, 'P=', p1
c		  print *, 'j=', j
c		  print *, 'KMATRIX=', KMATRIX(INTVS,j,8)
c		  endif
  	ENDDO  ! interval loop
	

        RETURN
        END
