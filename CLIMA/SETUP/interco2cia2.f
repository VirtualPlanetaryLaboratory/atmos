SUBROUTINE interpco2cia(Temp)
c	Written by Ramses Ramirez	
	PARAMETER(NF=55)	
	real tempgrid(7), FX, Temp, tempr, Tnum, Tden, IK, 
     &  CIA(7,NF), CPR(NF), CIAMS1L, CIAMSL, CPRL
	integer  M, MS, MS1
	COMMON/IRDATA/WEIGHT(8,3),XKAPPA(8,12,55,8,3), CIA(7,NF),  
     & CPR(NF)!c-rr 3/23/11 put CIA matrix in IRdata and Clima.f

	data tempgrid/100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0/


	tempr=Temp  ! Temperature around which CIAs are interpolated.

	DO M=1,7
		MS = M
		if (tempr.lt.tempgrid(M))then
			exit
		else
		MS = 8
		endif
	ENDDO

	if (MS.le.1) then   ! MS1 is left-most temperature grid point. MS is right-most temperature grid point
	   MS1=MS
	   FX=1   
	elseif (MS.ge.8) then
	   MS = MS - 1
	   MS1 = MS
	   FX =1 
	else
	   MS1 = MS-1
	
	   Tnum = tempr - tempgrid(MS1)  ! I am assuming temperature is not on a log scale? Otherwise Tnum and Tden are logged.
           Tden = tempgrid(MS) - tempgrid(MS1)
           FX = Tnum/Tden
	endif
	
	DO IK = 1,55   !do interpolation for all 17 coefficients at a given temperature tempr
			CIAMS1L = alog(CIA(MS1,IK))
			CIAMSL = alog(CIA(MS,IK))
			CPRL = CIAMS1L + FX*(CIAMSL-CIAMS1L)
			CPR(IK) = exp(CPRL)
	ENDDO

	END
