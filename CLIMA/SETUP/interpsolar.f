        SUBROUTINE interpsolar(Temp,p1,j,kmatrix_solco2,kmatrix_solh2o) 
c       by Ramses Ramirez 8/28/2012
c-rr    i is wavelength counter and j is the altitude counter used in solar, respectively.

        INCLUDE 'CLIMA/INCLUDE/header.inc'
        PARAMETER(NSOL=38,IK=8)
        REAL TempR,tempg(8),press(8),p1,kappa_solh2o, kappa_solco2, 
     &  kmatrix_solco2, kmatrix_solh2o,Temp
      DIMENSION KMATRIX_SOLH2O(NSOL,IK), KMATRIX_SOLCO2(NSOL,IK)
        INTEGER j !EWS - i not used
        DATA tempg/100., 150., 200., 250., 300., 350., 400., 600./ 
        DATA press/.00001, .0001, .001, .01, 0.1, 1., 10., 100./ 
       COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL),
     &  kappa_solh2o(NSOL,8,8,IK), kappa_solco2(NSOL,8,8,IK) ! Added new kappa matricies for each of CO2 and H2O coefficients. 8/26/2012 
      COMMON/DATASOLAR/weightco2_h2oSOL(8), weights(3,NSOL,IK)   ! new common block for weights and interpolated coefficients for CO2, H2O, and methane. 8/26/2012


c       Assigning TempR to temp (due to variable name overuse)
        TempR=Temp


C        There are separate loops for pressures lower than 1e-4 bar,  
C        for pressures higher than 1 bar, and those in between. For 
C        the pressures in between:

           j=j !EWS just to avoid compilation warnings
           Do L=1,8
                LS=L 
                 
             if (p1.le.press(L))then
                 exit!If an individial p greater than curent p set Ls=9 
                 else
             LS = 9
                 endif
           enddo
           

              Do M=1,8  ! For 8 temperatures
                 MS=M  
                if (TempR.lt.tempg(M))then
                exit  
                else
                MS=9
        endif                
              enddo
              


c       For pressures within the grid
        IF(LS.eq.1) THEN
          LS1 = LS        ! P1 is at or below the lowest grid pressure
        FY = 1.
       
        ELSEIF(LS.ge.8) THEN !        For pressures outside of grid:
         LS = 8                ! P1 is at or above the highest grid pressure
        LS1 = 8
        FY = 1.
        ELSE
        LS1 = LS-1
        PLOG=ALOG10(press(LS)) ! upper grid boundary in log units
        PLOG1=ALOG10(press(LS1)) ! lower grid boundary in log units
        P1L=ALOG10(p1) ! location of pressure within grid in log units
        FY=(P1L-PLOG1)/(PLOG-PLOG1) ! fractional distance of pressure(p1) within grid
        ENDIF
C



        IF(MS.le.1) THEN
          MS1=MS                !TempR is at or below the lowest grid Temperature
        FX = 1.
        ELSEIf(MS.ge.9)THEN ! For temperatures outside of the grid:
         MS=MS-1                !TempR is at or above the highest grid Temperature
        MS1=MS                
        FX = 1.
        ELSE
        MS1=MS-1  !This calculates MS1 for values within a temperature range (MS between 2 and 3)
        FX=(TempR-tempg(MS1))/(tempg(MS)-tempg(MS1)) ! This calculates FX for values within a temperature range
        ENDIF

          Do INTVS=1,NSOL  ! THE MASTER LOOP FOR THE 22 intervals.

           DO IK1=1,8  !IK is a kcoefficient index.1 to 8 because 8 kcoefficient sums.
            AK1L=LOG10(kappa_solh2o(INTVS,MS1,LS1,IK1))  !  These terms from H2O coefficient data table
            AK2L=LOG10(kappa_solh2o(INTVS,MS,LS1,IK1))
            AK3L=LOG10(kappa_solh2o(INTVS,MS,LS,IK1))
            AK4L=LOG10(kappa_solh2o(INTVS,MS1,LS,IK1))

            AK11L=LOG10(kappa_solco2(INTVS,MS1,LS1,IK1))  ! These terms from CO2 coefficient data table
            AK22L=LOG10(kappa_solco2(INTVS,MS,LS1,IK1))
            AK33L=LOG10(kappa_solco2(INTVS,MS,LS,IK1))
            AK44L=LOG10(kappa_solco2(INTVS,MS1,LS,IK1))

            AKL=FX*FY*AK3L+FX*(1.-FY)*AK2L+(1.-FX)*FY*AK4L+(1.-FX)*
     &                (1.-FY)*AK1L

            AKLL=FX*FY*AK33L+FX*(1.-FY)*AK22L+(1.-FX)*FY*AK44L+(1.-FX)*
     &                (1.-FY)*AK11L
         
            kmatrix_solh2o(INTVS,IK1)=10**AKL  ! H2O kmatrix.  This needs to be interpolated over altitude (j)
            kmatrix_solco2(INTVS,IK1)=10**AKLL  ! CO2 kmatrix. This needs to be interpolated over altitude (j)
            
!             if(INTVS .eq. 9)then
!                print *,kmatrix_solh2o(INTVS,IK1),Temp,p1
!                print *,IK1,MS,MS1,LS,LS1
!                pause
!             endif


  
!             print *,kmatrix_solco2(INTVS,IK1),IK1,INTVS
          ENDDO  ! ends AK loop


c            if ((J.eq.80).and.(INTVS.eq.17).and.(I.eq.20))then
c                  print *, 'T=',TempR
c                  print *, 'P=', p1
c                  print *, 'j=', j
c                  print *, 'KMATRIX=', KMATRIX(INTVS,j,8)
c                  endif

          ENDDO  ! interval loop
        

        RETURN
        END
