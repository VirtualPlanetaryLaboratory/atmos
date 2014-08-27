
c     SUBROUTINE IREXPSUMS(weight,xkappa)
      SUBROUTINE IREXPSUMS

c-mm  This subroutine will READ in the exponential sum datafiles in the
c-mm  ir for both h2o and co2 and output the appropriate value for each
c-mm  wavelength interval.  This will crash if our values exceed the
c-mm  boundaries set forth in the data, but this shouldn't happen under
c-mm  normal Martian conditions.
c-mm  Indices for xkappa and kappa are as follows
c-mm  Index #1:  Temp   1=less than T   2=more than T
c-mm  Index #2:  Pres   1=less than p   2=more than p
c-mm  Index #3:  lam    wavelength bin
c-mm  Index #4:  i      gauss point
c-mm  Index #5:  Spec   1=co2  2=h2o  3=ch4

      include 'CLIMA/INCLUDE/header.inc'
c      PARAMETER (NF=55,NGS=5)
      INTEGER LINE,i,lam,junk,Tind,pind,Tinddum,j,k,l,m
      REAL zeroed, CIA1(7,55)
c     REAL xkappa(8,12,55,8,3)
c     DIMENSION WEIGHT(8,3)
      CHARACTER :: DIRINOUT*2,DIRDATA*4
      COMMON/DIR/DIRINOUT,DIRDATA
C new common block, von Paris, 21/04/2006
      COMMON/IRDATA/WEIGHT(8,3),xkappa(8,12,55,8,3), 
     & CIA(7,55), CPRW(ND,55)  
c-rr !3/23/11 put CIA matrix in IRDATA

      OPEN(unit=20,file= DIRDATA//'/ir_expsums.pdat')
      do ii=1,3
       read(20,*)
      enddo
      do 9400 i=1,8
         do 9401 j=1,12
            do 9402 k=1,55
               do 9403 l=1,8
                  do 9404 m=1,3
                     xkappa(i,j,k,l,m)=0.0
 9404             continue
 9403          continue
 9402       continue
 9401    continue
 9400 continue
      READ(20,9950)
        do 9951 i=1,8
          READ(20,9952) weight(i,1)
 9951   continue
        do 9700 i=1,6
           READ(20,*)
 9700   continue

        do 9500 Tind=1,7
           do 9502 pind=1,11
              do 9957 lam=9,48

c-mm  If zeroed=0.0, that implies there is only one set of numbers
c-mm  for that wavelength bin.  Otherwise, there are two and we
c-mm  READ from the second set.  xkappa is the particular kappa
c-mm  (1-8) for the given indices.  This will be later interpolated

          READ(20,9956) zeroed
          READ(20,*)
          do 9955 i=1,8
             if (zeroed.ne.(0.0)) then
                READ(20,9958) xkappa(Tind,pind,lam,i,1)
             else
                READ(20,9960) xkappa(Tind,pind,lam,i,1)
             endif
 9955     continue
        READ(20,*)
 9957   continue
        READ(20,*)
 9502   continue
 9500   continue

c-mm   H2O below this line----------------------
        READ(20,9950)
        do 9901 i=1,8
           READ(20,9952) weight(i,2)
 9901       continue
        do 9701 i=1,6
           READ(20,*)
 9701      continue
        do 9504 Tind=1,7
           do 9506 pind=1,11
              do 9508 lam=1,55
 
c-mm  If zeroed=0.0, that implies there is only one set of numbers
c-mm  for that wavelength bin.  Otherwise, there are two and we
c-mm  READ from the second set.  xkappa is the particular kappa
c-mm  (1-8) for the given indices.  This will be later interpolated
 
          READ(20,9956) zeroed
          READ(20,*)
          do 9510 i=1,8
             if (zeroed.ne.(0.0)) then
                READ(20,9958) xkappa(Tind,pind,lam,i,2)
             else
                READ(20,9960) xkappa(Tind,pind,lam,i,2)
             endif
 9510     continue
        READ(20,*)
 9508 continue
        READ(20,*)
 9506 continue
 9504 continue

c-mm  CH4 below this line--------------------------------
      weight(1,3)=0.08566225
      weight(2,3)=0.18038079
      weight(3,3)=0.23395697
      weight(4,3)=0.23395697
      weight(5,3)=0.18038079
      weight(6,3)=0.08566225
      do ii=1,3 
       READ(20,*)
      enddo
      do 9600 Tinddum=1,3
         Tind=Tinddum*3-Tinddum/3
         do 9601 pind=7,12
c            do 9602 lam=1,38
            do 9602 lam=5,42
               READ(20,*) junk,(xkappa(Tind,pind,lam,ii,3), ii=1,6)
 9602          continue
              READ(20,*)
              READ(20,*)
              READ(20,*)
 9601          continue
 9600          continue

        close(20)

        do 9570 i=7,12
c        do 9571 j=1,38
        do 9571 j=5,42
        do 9572 k=1,6
          xkappa(1,i,j,k,3)=xkappa(1,i,j,k,3)
c          xkappa(1,i,j,k,3)=xkappa(3,i,j,k,3)
 9572 continue
 9571 continue
 9570 continue

C-rr reads CO2 CIA matrix(3/22/2011)----------------------
	read(30,9800) !skips two lines
	do k=1,9
		read(30,9801) (CIA1(it,k),it=1,7)  ! Put CIA values in the right bins 
	enddo
        

	do k=20,27
		read(30,9801) (CIA1(it,k),it=1,7)  ! Put CIA values in the right bins
	enddo

         
        do it=1,7
		do k=1,55
			CIA(it,k)=amax1(CIA1(it,k),1.e-60) ! set all zero CIA values to 1.e-60
		enddo
	enddo

        
 9800   format(/)
 9801   format(19x, 1pe24.17,6(3x,1pe24.17))
 9802   format(i3, 3x,1pe24.17,6(3x,1pe24.17))

 
 9950   format(/)
 9952   format(27X,E16.10)
 9956   format(83X,E9.3)
 9958   format(56X,E12.6)
 9960   format(21X,E12.6)
 
        RETURN
        END
