        SUBROUTINE interpir(TempR,p1,j, kappa,kmatrix_irco2,
     &  kmatrix_irh2o) 
c       by Ramses Ramirez 8/28/2012
c-rr    i is wavelength counter and j is the altitude counter used in ir, respectively.

        INCLUDE 'CLIMA/INCLUDE/header.inc'
        PARAMETER(NF=55,IK=8)
        INTEGER i,j,lam,Tindex,pindex,Tp1,pp1, Temp !EWS - LINE not used here
        REAL kappapp(55,8),kappapm(55,8) ! methane kappas
        REAL kappa(55,6)! methane kappas
        REAL TempR,tempg(8),press(8),p1,kappa_irh2o, kappa_irco2, 
     &  kmatrix_irco2, kmatrix_irh2o
        DIMENSION :: kmatrix_irco2(NF,IK), kmatrix_irh2o(NF,IK)
        DATA tempg/100., 150., 200., 250., 300., 350., 400., 600./ 
        DATA press/.00001, .0001, .001, .01, 0.1, 1., 10., 100./
        COMMON/IRDATA/WEIGHTCH4(6),xkappa(3,12,55,8), ! xkappa is pre-loaded methane absorption coefficient array 
     &   CIA(7,55), CPRW(ND,55)  !c-rr !3/23/11 put CIA matrix in IRDATA. Added weight co2_h2O 3/20/2012 
        COMMON/VARIR/kappa_irh2o(NF,8,8,IK), kappa_irco2(NF,8,8,IK)! Added kappa matrix in IR for kpsectrum Co2 and H2O coefficients 8/26/2012  
        COMMON/weightsIR/weightco2_h2oIR(IK)


C        There are separate loops for pressures lower than 1e-4 bar,  
C        for pressures higher than 1 bar, and those in between. For 
C        the pressures in between:
 
           j=j !EWS -avoid complaining         
           Do L=1,8
                LS=L 
                 
             if (p1.le.press(L)) then
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
        PLOG=ALOG10(press(LS))! upper grid boundary in log units
        PLOG1=ALOG10(press(LS1))! lower grid boundary in log units
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




            Do INTVS=1,NF  ! THE MASTER LOOP FOR THE 22 intervals.




           DO IK1=1,8  !IK is a kcoefficient index.1 to 8 because 8 kcoefficient sums.

!             if((INTVS .eq. 15) .and. (IK1 .eq. 8))then
!          print *, 'MS1=', MS1, 'MS=', MS, 'LS1=', LS1, 'LS=', LS
!          print *, 'Input T=', TempR, 'Input P=', p1, 'FX=', FX
!          pause
!            endif

            AK1L=LOG10(kappa_irh2o(INTVS,MS1,LS1,IK1))  !  These terms from H2O coefficient data table
            AK2L=LOG10(kappa_irh2o(INTVS,MS,LS1,IK1))
            AK3L=LOG10(kappa_irh2o(INTVS,MS,LS,IK1))
            AK4L=LOG10(kappa_irh2o(INTVS,MS1,LS,IK1))

            AK11L=LOG10(kappa_irco2(INTVS,MS1,LS1,IK1))  ! These terms from CO2 coefficient data table
            AK22L=LOG10(kappa_irco2(INTVS,MS,LS1,IK1))
            AK33L=LOG10(kappa_irco2(INTVS,MS,LS,IK1))
            AK44L=LOG10(kappa_irco2(INTVS,MS1,LS,IK1))

            AKL=FX*FY*AK3L+FX*(1.-FY)*AK2L+(1.-FX)*FY*AK4L+(1.-FX)*
     &                (1.-FY)*AK1L

            AKLL=FX*FY*AK33L+FX*(1.-FY)*AK22L+(1.-FX)*FY*AK44L+(1.-FX)*
     &                (1.-FY)*AK11L
         
            kmatrix_irh2o(INTVS,IK1)=10**AKL  ! H2O kmatrix.  This needs to be interpolated over altitude (j)
            kmatrix_irco2(INTVS,IK1)=10**AKLL  ! CO2 kmatrix. This needs to be interpolated over altitude (j)
 
!             print *,kmatrix_irh2o(INTVS,IK1),INTVS,IK1

             
!             if((INTVS .eq. 15) .and. (IK1 .eq. 8))then
!                 print *,kmatrix_irco2(INTVS,IK1),AKLL
!                  print *,AK11L,AK22L,AK33L,AK44L, 
!     &            'Press=', p1, 'Temp=', TempR
!                  pause
!            endif

          ENDDO  ! ends AK loop

!          print *,Temp,p1
!            if ((J.eq.80).and.(INTVS.eq.17).and.(I.eq.20))then
c                  print *, 'T=',TempR
c                  print *, 'P=', p1
c                  print *, 'j=', j
c                  print *, 'KMATRIX=', KMATRIX(INTVS,j,8)
c                  endif

          ENDDO  ! interval loop


!--------------------------------------------------------------------

      do 9400 i=1,NF
         do 9401 ij=1,6
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
