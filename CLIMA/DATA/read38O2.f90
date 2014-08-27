Program read38O2

!DIMENSION :: BETAO3(14)
COMMON/BLAH/BETAO3(14), d(14), c(14)
	OPEN(unit=4, file='solar_data_38.pdat', status='old')


		DO I = 1, 14
   		READ(4,*) 
   		READ(4,*) 
   		READ(4,*) a, b, c(I), d(I)
	        READ(4,350)
                if (I.le.8)then ! For intervals 1-8, 3rd column is O3
                BETAO3(I) = c(I)
                print *, BETAO3(I)
                else  ! For intervals 9-14, 4th column is O3
                BETAO3(I) = d(I)
                print *, BETAO3(I)
                endif
                ENDDO
                
               read(4,360)
               read(4,*) a,b,f,BETAO2
               print *, BETAO2

 350 format(2(/))
 360 format(7(/))

END ! End program
