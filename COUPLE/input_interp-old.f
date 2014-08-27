    
      SUBROUTINE INPUT_INTERP(P,pphot, water, O3, Jcold, FI)

c This subroutine interpolates the water and  ozone profiles that are read
c  in the climate model and come from the photochemical model
      
      INCLUDE '../CLIMA/INCLUDECLIM/parND.inc'
      INCLUDE '../ATMCHEM/INCLUDECHEM/parNZ.inc'

      DIMENSION O3(NZ), water(NZ), pphot(NZ),FI(4,ND),P(ND)
      
      do i=1,NZ
      pphot(i) = pphot(i)/1.e6  !dyn/cm^2 to bar
      enddo
       
	istart = 1
	
	DO j = ND, 1, -1
	 DO i = istart, NZ-1
	   IS = i
	   if (pphot(i) .LT. P(j)) GOTO 350
	 END DO
  350  CONTINUE
	 IS1 = IS+1
 	 istart = IS1
	 DZI = pphot(IS1)- pphot(IS)
	 DZJ = P(j) - pphot(IS)
	 FR = DZJ/DZI
 	 O3log1 = ALOG(O3(IS))
	 O3log2 = ALOG(O3(IS1))
	 O3log = (FR*O3log2) + ((1-FR)*O3log1)
	 FI(4,j) = EXP(O3log)
 	 Flog1 = ALOG(water(IS))
	 Flog2 = ALOG(water(IS1))
	 flog = (FR*Flog2) + ((1-FR)*Flog1)
	 FI(1,j) = EXP(flog)
	END DO
        jcold_new = 1

	DO k = 1, ND
	  if (P(k) .LE. pphot(JCOLD)) jcold_new = k
 	END DO
       JCOLD = jcold_new
        
	RETURN
        END

