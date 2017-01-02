 
*  Created by Antigona Segura. April 4, 2005.
       
       program couple

*  This program runs the climate and photochemical model in a coupled mode
*  and generates an output file for diagnostics
       INCLUDE 'CLIMA/INCLUDE/header.inc'
       INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'

       CHARACTER :: STARR*3,dirIO*2

c Common blocks with CLIMATE model
      INCLUDE 'INCLUDE/comCLIM.inc'

c Common bloCKs with PHOTOCHEMICAL model
      INCLUDE 'INCLUDE/comPHOT.inc'
      INCLUDE 'INCLUDE/comH2BALANCE.inc'
      INCLUDE 'INCLUDE/comFLXLOW.inc'
      INCLUDE 'INCLUDE/comFLUXPHOTO.inc'
      INCLUDE 'INCLUDE/comDIAG.inc'
      INCLUDE 'INCLUDE/comSTR.inc'
      COMMON/SPECIES/FCH4SAVE,FC2H6SAVE,FCO2SAVE,
     &  FO2SAVE,FNO2SAVE,FO3SAVE,FNH3SAVE 

      dirIO = 'IO'
c INPUT files
c     OPEN(unit=1,file= dirIO//'/input_clima.dat')
c OUTPUT FILES
      open(unit=301, file= dirIO//'/output_couple.dat')
      OPEN(UNIT=198,FILE= dirIO//'/clima_allout.tab')
      OPEN(unit=64,file= dirIO//'/photchem_allout.dat') 
        OPEN(unit=291,file='IO/mixing_ratios.dat')
        READ(291,*) X                     !Argon
        READ(291,*) FCH4SAVE                 !Methane
        READ(291,*) FC2H6SAVE                !Ethane
        READ(291,*) FCO2SAVE                 !Carbon dioxide
        READ(291,*) FO2SAVE                  !Oxygen
        READ(291,*) FNO2SAVE                 !Nitrogen dioxide
        READ(291,*) FO3SAVE                  !Ozone
        READ(291,*) X                !Tropopause layer
        READ(291,*) FNH3SAVE                 !Ammonia
        close(291)
c READ PLANET PARAMETERS
      OPEN(unit=204,file='IO/PLANET.dat',status='old')
      READ(204,502) GRAV,X,PG1,SOLCON,X,X
 502  FORMAT(F9.2/,E9.3/,E9.3/,F9.2/,F9.3/,E9.2)
      PG1=PG1/1.e6
      CLOSE(204)

      IMODEL = 3
      niter = 5
      STARR='Sun'
      ICOUPLE = 1 
c Choose a star
      write(*,*)'Choose a star: Sun or  dMV'
      write(*,*)'Write it exactely as is shown'
      read(*,'(A3)') STARR
c Choose programs
      write(*,*)'Do you want to run the coupled mode?'
      write(*,*)'0- no, 1-yes'
      read(*,*) ICOUPLE
      
      if(ICOUPLE.eq.0) then
         write(*,*)'What model do you want to use?'
         write(*,*)'0-climate, 1-photochemical'
         read(*,*) IMODEL
      endif
      write(*,*)'How many model iterations?'
      read(*,*) niter
     
 
      if(IMODEL.ne.1) then
c Parameters to run the climate model
         write(*,*)'For the climate model:'
         write(*,*)'Number of steps (usually 100)'
         read(*,*)NSTEPSC
      endif
c Photochemical model 
      if(IMODEL.eq.1.or.ICOUPLE.eq.1) then
         write(*,*)'For the photochemical model:'
         write(*,*)'Number of steps (usually 500)'
         read(*,*)NSTEPSP
      endif
      print*,'end reading parameters'
c Printing the atmospheric parameters
      write(301,990)STARR
      write(301,993) FN2SAVE,FO2SAVE,FCO2SAVE,FCH4SAVE, 
     & FNO2SAVE,FNH3SAVE
      write(301,995), GRAV, PG1, SOLCON
 990  format(2x,'PARAMETERS FOR THE PLANET AROUND ', A4/)
 993  format (2x,'FN2= ',E10.3,3x, 'FO2= ',E10.3,3x, 'FCO2=',E10.3/,
     & 2x,'FCH4= ',E10.3,3x, 'FNO2= ',E10.3,2x,'FNH3= ',E10.3/)
 995  format (2x,'Gravity= ',F9.2,2x, 'Surface pressure= ',F6.2,
     & 2x,'SolCon= ',F6.3/)
      
c Starting loop for steady-state equilibrium
      do 1 ncouple=1,niter   !START steady-state loop
       print*,''
       print*,'loop number ',ncouple
       print*,''
c Running the climate code
       if(IMODEL.eq.0.or.ICOUPLE.eq.1) then
          DTC = 5.e3      !initial time step (s)
          dtmax = 1.e4      !maximum time step (s)
          PRINT *,"Calling CLIMA:"
          CALL CLIMA(ICOUPLE,DTC,NSTEPSC,dtmax,ncouple)
       print*,'end climate'
       endif
       print *,'icouple is',icouple
c Running the photochemical code
       if(IMODEL.eq.1.or.ICOUPLE.eq.1) then
          DTP = 1.e-6      !initial time step
          TSTOP = 1.e17
          CALL PHOTOCHEM(ICOUPLE,TSTOP,DTP,NSTEPSP,TIME,ncouple)
       endif
c Calculating diagnostic parameters for stady-state convergence
       if (ICOUPLE.eq.1) then
       call diagnostic(nconvdif)
       if(nconvdif.ne.0) then
          close(198)
          close(64)
          OPEN(UNIT=198,FILE= dirIO//'/clima_allout.tab')
          OPEN(unit=64,file= dirIO//'/photchem_allout.dat')
       endif
       if(ncouple.eq.niter.and.nconvdif.ne.0) then
          write(301,*)'** WARNING: STEADY STATE SOLUTION NOT REACHED **'
          call outdiag(ncouple)
          call outatmos
          print*,'*** WARNING: STEADY STATE SOLUTION NOT REACHED ***'
          STOP
       endif
       if(nconvdif.eq.0) go to 2
       endif     
  1    enddo                  !END steady-state loop 
  
c Printing diagnostic parameters for the steady solution
  2    call outdiag(ncouple)
       call outatmos
c      stop

  3    close(198)
       close(64)
       close(301)
       stop
       end          

***********************
      subroutine diagnostic(nconvdif)
c This subroutine indicates when the photochemical and
c climate model have reached convergence
       INCLUDE 'CLIMA/INCLUDE/header.inc'
       INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
c Common blocks with CLIMATE model
      INCLUDE 'INCLUDE/comCLIM.inc'

c Common blocks with PHOTOCHEMICAL model
      INCLUDE 'INCLUDE/comPHOT.inc'
      INCLUDE 'INCLUDE/comH2BALANCE.inc'
      INCLUDE 'INCLUDE/comFLXLOW.inc'
      INCLUDE 'INCLUDE/comDIAG.inc'
      COMMON/SPECIES/FCH4SAVE,FC2H6SAVE,FCO2SAVE,
     &  FO2SAVE,FNO2SAVE,FO3SAVE,FNH3SAVE 
      nconver = 0
      nconv1 = 0
      nconv2 = 0
      nconv3= 0

c equilibrium conditions for the photochemical model
c   *** for species with fixed surface flux
      numlb2tot = 0
      numlb2 =0
      do i=1,NQ
       if(LLBOUND(i).eq.2)then
        nconv1 = 1
        numlb2tot = numlb2tot + 1
        dflux(i) = XFLOW(i)/XSGFLUX(i)
        if(dflux(i).ge.0.99.and.dflux(i).le.1.01) 
     &   numlb2 = numlb2 + 1
       endif 
      enddo
      if(numlb2.gt.0) then
        perlb2= numlb2/numlb2tot
        if(perlb2.eq.1.)then
         nconver=nconver+1
         print*,'Eq. for fixed fluxes, reached'
        endif
      endif 
c   *** for species with fixed mixing ratio
      numlb1tot = 0
      numlb1 =0
      do i=4,NQ
       if(LLBOUND(i).eq.1) then
        nconv2 = 1
        numlb1tot = numlb1tot + 1
        dlflux(i) = XFLOW(i)/XTL(i)
        if(i.eq.2)then
          dmixo2 = XUSOL(i,1)/FO2SAVE
          if(dmixo2.ge.0.99.and.dmixo2.le.1.01) 
     &    numlb1 = numlb1 + 1
        endif
        if(i.eq.12)then
          dmixch4 = XUSOL(i,1)/FCH4SAVE
          if(dmixch4.ge.0.99.and.dmixch4.le.1.01) 
     &    numlb1 = numlb1 + 1
        endif
        if(i.eq.26)then
          dmixnh3 = XUSOL(i,1)/FNH3SAVE
          if(dmixnh3.ge.0.99.and.dmixnh3.le.1.01)
     &    numlb1 = numlb1 + 1
        endif
       endif
      enddo
      if(numlb1.gt.0)then
       perlb1= numlb1/numlb1tot
       if(perlb1.eq.1.) then
         nconver=nconver+1
         print*,'Eq. for fixed abundances reached'
       endif
      endif
c   *** for species with fixed velocity deposition
      numlb0tot = 0
      numlb0 =0
      do i=1,NQ
       if(LLBOUND(i).eq.0) then
        nconv3 = 1
        numlb0tot = numlb0tot + 1 
        dlospro(i) = XTL(i)/XTP(i)
         if(dlospro(i).ge.0.9.and.dlospro(i).le.1.1) 
     &   numlb0 = numlb0 + 1
       endif
      enddo      
      if(numlb0.gt.0) then 
 5       Aperlb0= (numlb0*1.)/(numlb0tot*1.)
        print*,'vel dep. ',perlb0
c       if(perlb0.ge.0.9) then
        if(perlb0.ge.0.7) then    !for AD Leo
        nconver=nconver+1
        print*,'Eq. for fixed deposition velocity, reached'
       endif
      endif

c For the H2 balance
      ratioH2= (H2SURF+H2VOLC)/(H2CHEM+H2ESC)
      dH2BAL= ABS(1.-ratioH2)
      print*,'H2 balance',dH2BAL
      if(dH2BAL.lt.1.e-1) then 
        nconver =nconver +1
       print*,'Eq. for H2 balance, reached'
      endif

C Equilibrium conditions for the climate code
c  *** Temperature at the surface, the top, the last convective
c   layer and the tropopause
c Ratio of fluxes at the top of the atmosphere (DIVF)
      do i=ND,1,-1
        dtn(i) = abs(Tstart(i)-T(i))
      enddo
      If(dtn(ND).le.1.5) nconver = nconver + 1
      if(dtn(1).le.3.) nconver = nconver + 1
      if(dtn(JCOLD).le.2.)nconver = nconver + 1
      if(dtn(JCONV).le.2.)nconver = nconver + 1
      if(divf(1).le.1e-3)nconver = nconver + 1
C Summing up the equilibrium conditions
C When all the conditions are reached nconvdif = 0   
      nconvtot = nconv1 + nconv2 + nconv3 + 6
      nconvdif = nconvtot - nconver
      write(*,110)nconver, nconvtot,nconvdif
 110  format(I2,' convergence criteria of ',I2,' have been reached. ',
     & I2,' more need to be reached')
      return
      end
*******************

      subroutine outdiag(ncouple)
       INCLUDE 'CLIMA/INCLUDE/header.inc'
       INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'

      INCLUDE 'INCLUDE/comCLIM.inc'

      INCLUDE 'INCLUDE/comPHOT.inc'
      INCLUDE 'INCLUDE/comH2BALANCE.inc'
      INCLUDE 'INCLUDE/comFLXLOW.inc'
      INCLUDE 'INCLUDE/comDIAG.inc'

      COMMON/SPECIES/FCH4SAVE,FC2H6SAVE,FCO2SAVE,
     &  FO2SAVE,FNO2SAVE,FO3SAVE,FNH3SAVE 

c Printing diagnostic parameters from the steady state solution 
      write(301,100)ncouple
 100  format(/'*** Output file for diagnostics after ',I2,
     &' coupled iterations ***')
      write(301,120)
 120  format(/'*** PHOTOCHEMICAL MODEL output'/)
      write(301,310)
 310  FORMAT(/1X,'H2 BUDGET BALANCE (PHOTO MODEL)'/)

      write(301, 311) H2SURF,H2VOLC,H2CHEM,H2ESC,H2BAL
 311  FORMAT(5X,'H2SURF =',1PE10.3,3X,'H2VOLC =',E10.3,3X,'H2CHEM =',
     2  E10.3,4X,'H2 ESCAPE =',E10.3,4X,'H2BAL = ',E10.3/)
      ratioH2= (H2SURF+H2VOLC)/(H2CHEM+H2ESC)
      dH2BAL= ABS(1.-ratioH2)

      write(301,312) dH2BAL
 312  format(2x,'Hydrogen balance = ',1E10.4/)

      write(301,121)
 121  format(/' Results sorted by Lower Boundary Condition'/) 
      
      ifix = 0
      do i=1,NQ
       if(LLBOUND(i).eq.2) then
         ifix=1
         goto 210
       endif
      enddo
 210  if(ifix.eq.1) then
      write(301,*)' Species with fixed surface FLUX'
      write(301,123)
 123  format(6x,'Specified flux',3x,'Calculated flux',4x,'calc/spec'
     & ,6x, 'TL') 
      do i=1,NQ
       if(LLBOUND(i).eq.2) write(301,130)IISPEC(i),XSGFLUX(i),
     & XFLOW(i),dflux(i),XTL(i)
      enddo
 130  format(1x,A8,2x,1PE9.2,4x,1PE9.2,9x,1PE9.2,3x,1PE9.2)
      endif
      imix = 0
      do i=4,NQ
       if(LLBOUND(i).eq.1) then
         imix=1
         goto 220
       endif
      enddo
 220  if(imix.eq.1) then
      write(301,*)
      write(301,*)' Species with fixed surface MIXING RATIO'
      write(301,124)
 124  format(7x,'Specified M R',3x,'Calculated M R',3x,'calc/spec',
     & 6x,'TL',9x,'TP',9x,'TP/TL') 
      do i=4,NQ
       if(LLBOUND(i).eq.1) then
        if(i.eq.2)write(301,131)IISPEC(i),FO2SAVE,XUSOL(i,1),dmixo2,
     & XFLOW(i),XTL(i), dlflux(i)
        if(i.eq.12)write(301,131)IISPEC(i),FCH4SAVE,XUSOL(i,1),dmixch4,
     & XFLOW(i),XTL(i), dlflux(i)
        if(i.eq.26)write(301,131)IISPEC(i),FNH3SAVE,XUSOL(i,1),dmixnh3,
     & XFLOW(i),XTL(i), dlflux(i)

       endif
      enddo
 131  format(1x,A8,2x,1PE9.2,4x,1PE9.2,6x,1PE9.2,3x,1PE9.2,
     & 3x,1PE9.2,3x,1PE9.2) 
      endif

      write(301,*)
      write(301,*)' Species with fixed VELOCITY DEPOSITION'
      write(301,125)
 125  format(15x,'TL',9x,'TP',9x,'TL/TP')
      do i=1,NQ
       if(LLBOUND(i).eq.0)write(301,132)IISPEC(i),XTL(i),
     &  XTP(i),dlospro(i)
      enddo      
 132  format(1x,A8,2x,1PE9.2,3x,1PE9.2,3x,1PE9.2) 

      write(301,150)
 150  format(/'*** CLIMATE MODEL output'/)
      write(301,*)'  Selected temperatures for diagnostic'
      write(301,151)
 151  format(1x,'Layer',4x,'T_Ncoup-1',7x,'T_Ncoup',7x,'(T_n-1)-(T_n)'
     & ,5x,'DIVF_Nstep')
      do i=ND,1,-1
       if(i.eq.ND.or.i.eq.ND-1)write(301,152)i,Tstart(i),T(i),dtn(i),
     & DIVF(i)
       if(i.eq.JCOLD)then
        write(301,*)' At the cold trap '
        write(301,152)i,Tstart(i),T(i),dtn(i),DIVF(i)
 152  format(1x,I3,4(4x,1PE12.4))
       endif
       if(i.eq.JCONV)then
        write(301,*)' Last convective layer '
        write(301,152)i,Tstart(i),T(i),dtn(i),DIVF(i)
       endif
       if(i.eq.1.or.i.eq.2)write(301,152)i,Tstart(i),T(i),dtn(i),
     & DIVF(i)
      enddo
 
      return
      end
**************************
      subroutine outatmos
c 
       INCLUDE 'CLIMA/INCLUDE/header.inc'
       INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'

      dimension z(nz),WAV(108),WAVEUV(10),PRESS(NZ)
      INCLUDE 'INCLUDE/comCLIM.inc'

      INCLUDE 'INCLUDE/comPHOT.inc'
      INCLUDE 'INCLUDE/comFLUXPHOTO.inc'
      INCLUDE 'INCLUDE/comFLXLOW.inc'

      DATA WAV/1762.,1778.,1794.,1810.,1827.,1844.,1861.,1878.,
     & 1896.,1914.,1933.,1952.,1971.,1990.,2010.,2031.,2052.,2073.,
     & 2094.,2117.,2140.,2163.,2187.,2211.,2235.,2260.,2286.,2313.,
     & 2340.,2367.,2396.,2425.,2454.,2485.,2516.,2548.,2581.,2615.,
     & 2650.,2685.,2722.,2759.,2798.,2837.,2878.,2920.,2963.,3008.,
     & 3054.,3101.,3150.,3200.,3250.,3300.,3350.,3400.,3438.,3500.,
     & 3600.,3700.,3800.,3900.,4000.,4100.,4200.,4300.,4400.,4500.,
     & 4600.,4700.,4800.,4900.,5000.,5100.,5200.,5300.,5400.,5500.,
     & 5600.,5700.,5800.,5900.,6000.,6100.,6200.,6300.,6400.,6500.,
     & 6600.,6700.,6800.,6900.,7000.,7100.,7200.,7300.,7400.,7500.,
     & 7600.,7700.,7800.,7900.,8000.,8100.,8200.,8300.,8400.,8500/
c 10 Far UV wavelengths at the lower limit of the interval. 
c The last wavelegth is for Lyman alpha
      data WAVEUV/1725.,1675.,1625.,1575.,1525.,1475.,1425.,
     &            1375.,1325.,1216./

      ifile = 301
      BK = 1.38E-16

      write(301,169) 
 169  format(/2x,'** PARAMETERS OF THE ATMOSPHERE AT STEADY STATE **')

 172  format(/3X,'*** Surface fluxes (cm^-2 s^-1) ***'/)
      write(ifile,173)IISPEC(8),IISPEC(12)
      write(ifile,174)XFLOW(8),XFLOW(12)
 173  FORMAT(5X,6(A6,5X))
 174  FORMAT(6(1X,1PE10.3))
 
      write(ifile,175)
 175  FORMAT(/3X,'*** Mixing ratios ***'/)
c printing mixing ratios of: H2O,H2,CH4
      write(ifile,176)IISPEC(2),IISPEC(3),IISPEC(8),IISPEC(12)
 176  FORMAT(1X,'Z (km)',1x,'P(dyn/cm^2)',3X,6(A8,2X))
      do i=1,NZ
c calculating pressure in dyn/cm2
        PRESS(i)=XDEN(i)*BK*XTCH(i)
        z(i) = i - 1. + 0.5
        write(ifile,177) Z(I),PRESS(i),XUSOL(2,i),XUSOL(3,i),
     &   XUSOL(8,i),XUSOL(12,i)
      enddo
 177  FORMAT(1x,f5.2,2(2X,1P7E10.3))
      

      write(ifile,178)
 178  format(/3X,'*** Number densities (cm^-3) ***'/)
      write(ifile,179) IISPEC(2)
 179  format(2x,'Z(km)',1x,'P(dyn/cm^2)',3x,'Total',2x,2(A8,4X))
      do i=1,NZ
        write(ifile,177) Z(I),PRESS(i),XDEN(i),XSL(2,i)
      enddo
      
      write(ifile,180)
 180  FORMAT(/4X,'ENERGY FLUXES IN W/m^2/nm and photons/m^2/s/nm
     & (NOT DIURNALLY AVERAGED)')       
      write(ifile,181)
 181  FORMAT(/8X,'WAV',8X,'EFLUX',8X,'GFLUX',8x,'PhEFLUX',
     & 6x,'PhGFLUX')
      do jj=10,1,-1
        feflux = ESFX(jj)
        fgflux = GSFX(jj)
        fpheflux = PhESFX(jj)
        fphgflux = PhGSFX(jj)
        write(ifile,182)WAVEUV(jj),feflux,fgflux,fpheflux,fphgflux
      enddo
      do jj=1,108
        feflux = EFLUX(jj)
        fgflux = GFLUX(jj)
        fpheflux = PhEFLUX(jj)
        fphgflux = PhGFLUX(jj)
        write(ifile,182)WAV(jj),feflux,fgflux,fpheflux,fphgflux
      enddo
 182  FORMAT(2X,1P7E13.5)
      
      write(ifile,190)
 190  format(//3x,'**** From the climate model ****'/)
      write(ifile,191)
 191  format(5x,'P (atm)',5x,'Alt (km)',9x,'T',9x,'DIVF')
      do i=1,ND
       write(ifile,192) XP(i),XALT(i),T(i),DIVF(i)
      enddo
 192  format(4(2x,1pe11.4))
      
      return
      end


**************************


