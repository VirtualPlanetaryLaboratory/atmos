      SUBROUTINE pickstar(STARR,SOLINT)
! Written by Ramses Ramirez 

      PARAMETER(NSOL=38)
      DIMENSION SOLINT(NSOL)
      CHARACTER(5) :: STARR*5
      real :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11, 
     & s12,s13,s14,s15,s16,s17,s18,s19,s20,  
     & s21, s22, s23, s24, s25, s26, s27,  
     & s28, s29, s30,TOTAL
      integer :: a,WL1, Wl2 



      read(10,*) 
        do I =1,NSOL
	  read(10,*)a,WL1,WL2,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 
     & ,s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24 
     & ,s25, s26, s27, s28, s29, s30 
          IF (STARR== "Sun") THEN
             SOLINT(I)=s1
!	     print *, SOLINT(I)
	  ELSE IF (STARR=="GJ581")THEN
             SOLINT(I)=s2
	     print *, SOLINT(I)
	  ELSE IF (STARR=="ADLEO")THEN
	     SOLINT(I)=s3
             print *, SOLINT(I)
	  ELSE IF (STARR =="GJ644") THEN
	     SOLINT(I)=s4
	     print *, SOLINT(I)
	  ELSE IF (STARR =="M2600") THEN
	     SOLINT(I)=s5
	     print *, SOLINT(I)
	  ELSE IF (STARR =="M2700") THEN
	     SOLINT(I)=s6
	     print *, SOLINT(I)
	  ELSE IF (STARR =="M2800") THEN
	     SOLINT(I)=s7
	     print *, SOLINT(I)
      ELSE IF (STARR =="M3000") THEN
	     SOLINT(I)=s8
	     print *, SOLINT(I)
      ELSE IF (STARR =="M3200") THEN
	     SOLINT(I)=s9
	     print *, SOLINT(I)
      ELSE IF (STARR =="M3400") THEN
	     SOLINT(I)=s10
	     print *, SOLINT(I)
      ELSE IF (STARR =="M3600") THEN
	     SOLINT(I)=s11
	     print *, SOLINT(I)
      ELSE IF (STARR =="M3800") THEN
	     SOLINT(I)=s12
	     print *, SOLINT(I)
	  ELSE IF (STARR =="K4000") THEN
	      SOLINT(I)=s13
	     print *, SOLINT(I)
      ELSE IF (STARR =="K4200") THEN
	     SOLINT(I)=s14
	     print *, SOLINT(I)
      ELSE IF (STARR =="K4400") THEN
	     SOLINT(I)=s15
	     print *, SOLINT(I)
      ELSE IF (STARR =="K4600") THEN
	     SOLINT(I)=s16
	     print *, SOLINT(I)
      ELSE IF (STARR =="K4800") THEN
	     SOLINT(I)=s17
	     print *, SOLINT(I)
      ELSE IF (STARR =="K5000") THEN
	     SOLINT(I)=s18
	     print *, SOLINT(I)		
      ELSE IF (STARR =="K5200") THEN
	     SOLINT(I)=s19
	     print *, SOLINT(I)		 
      ELSE IF (STARR =="G5400") THEN
	     SOLINT(I)=s20
	     print *, SOLINT(I)
      ELSE IF (STARR =="G5600") THEN
	     SOLINT(I)=s21
	     print *, SOLINT(I)
      ELSE IF (STARR =="G5800") THEN
	     SOLINT(I)=s22
	     print *, SOLINT(I)
      ELSE IF (STARR =="G6000") THEN
	     SOLINT(I)=s23
	     print *, SOLINT(I)
      ELSE IF (STARR =="F6200") THEN
	     SOLINT(I)=s24
	     print *, SOLINT(I)
      ELSE IF (STARR =="F6400") THEN
	     SOLINT(I)=s25
	     print *, SOLINT(I)
      ELSE IF (STARR =="F6600") THEN
	     SOLINT(I)=s26
	     print *, SOLINT(I)
      ELSE IF (STARR =="F6800") THEN
	     SOLINT(I)=s27
	     print *, SOLINT(I)
      ELSE IF (STARR =="F7000") THEN
	     SOLINT(I)=s28
	     print *, SOLINT(I)
      ELSE IF (STARR =="F7200") THEN
	     SOLINT(I)=s29
	     print *, SOLINT(I)
      ELSE IF (STARR =="F7400") THEN
	     SOLINT(I)=s30
	     print *, SOLINT(I)
          ENDIF
		  
        enddo
      close(10)
c	     DO I=1,38
c	      print *, SOLINT(I)
c         ENDDO
        TOTAL=0.
	DO I=1,NSOL
           TOTAL=SOLINT(I)+TOTAL
        ENDDO
        PRINT *,'TOTALFLUX=',TOTAL

	END
