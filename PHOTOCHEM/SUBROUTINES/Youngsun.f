      SUBROUTINE youngsun(n,timega,grid,fluxmult)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
c  this subroutine calculates relative flux increments for each model wavelength
c based on the time before the present (timeGa)


!INPUTS
!n is number of elements in grid
!grid is wavelength grid (in A)
!timega is time in billion years ago (Ga)

!OUTPUTS
!fluxmult is the relative multiplier defined on the Thuillier et al grid.

      ! option 1 is using Kevin's original perscription 4X <1754 and 2X<2500
      ! option 2 uses Claire et al. 2012   Sept 30, 2012  Version 3.141

      REAL*8 grid(n),nm(n),fluxmult(n)
      REAL*8 chromofit(n),chromomult(n),relphoto(n)

      REAL*8 kuruczflux(43,1221),AGE(43),wave(1221),dummy(44,1221)
c     dummy is read from 'DATA/FLUX/kuruczflux.dat', which has 44 columns and 1221 rows. 44th column is wavelength, the first 43 columns are fluxes corresponding to the different solar ages for that row's wavelength. Final row is the age grid in increments of ~ 0.2 Gyr. Dummy is then used to assign the proper values to wave and kuruczflux in an embedded do-loop below in (1)

      REAL*8 photothen(1221), photonow(1221),relphotokuruczgrid(1221) !1221=number of wavelength elements that have corresponding photspheric fluxes on the Kurucz wavelength grid

      REAL*8 stronglines(16), slinebetas(16)

      REAL*8 RnormBahcall(41),LnormBahcall(41),BahcallAge(41)
      REAL*8 temp(41) !used by cubic spline
      REAL*8 Rnorm, Lnorm,now, LnormNOW, RnormNOW, lowerl, upperl !EWS - added some extra variable declarations

c      INTEGER diff !EWS - not currently used
      !used to calculate offset between wavl and chromomult/fit vectors, 
      !need to be an integer because it is used below in (9) 
      !to properly index those arrays

      DATA RnormBahcall/0.869,0.882,0.888,0.892,0.897,0.901,0.906,0.91,
     2  0.915,0.92,0.924,0.929,0.934,0.94,0.945,0.951,0.956,0.962,0.968,
     3  0.974,0.981,0.987,0.994,1.001,1.008,1.016,1.023,1.031,1.04,
     4  1.048,1.057,1.066,1.075,1.085,1.095,1.105,1.116,1.127,1.139,
     5  1.152,1.166/

      DATA LnormBahcall/0.677,0.721,0.733,0.744,0.754,0.764,0.775,0.786,
     2  0.797,0.808,0.82,0.831,0.844,0.856,0.869,0.882,0.896,0.91,0.924,
     3  0.939,0.954,0.97,0.986,1.003,1.02,1.037,1.055,1.073,1.092,1.112,
     4  1.132,1.152,1.172,1.193,1.214,1.235,1.256,1.278,1.304,1.332,
     5  1.363/
      
      DATA BahcallAge/0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,
     2  2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,
     3  5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0/

c      print *, 'testing in Youngsun.f'
      
      option=2
      
      if (option.eq.1) then
      !Zahnle et al. 2006
      do L=1,n
       if (grid(L).LE.2500.0 .AND. grid(L).GE.1754.) then 
         fluxmult(L)=2.0  
       else if (grid(L).LT.1754.) then 
        fluxmult(L)=4.0 
       else
        fluxmult(L)=1.0
       endif

       !testing reducing photospheric flux for faint young sun...
      relL=(1.0 / (1.0 + 0.0835 * timega))
      if (grid(L).gt.3000) then
c         fluxmult(L)=relL
      endif   

      enddo

      endif  !end option 1



      if (option.eq.2) then  !use Claire et al. 2012

      a=-1.21374
      b=0.000502690
      now=4.56  !age of the sun in Gyr

!initialize vectors (MAYBE NOT NEEDED)
      do L=1,n
         relphoto(L)=0.0
         chromofit(L)=0.0
         chromomult(L)=0.0
         nm(L)=grid(L)/10.  !grid converted to Angstroms for ease later...
      enddo

      do L=1,n
!first, initialize flux multipliers in shorwave continuum (will be updated between 1500/2000 later)
         if (grid(L).lt.2000) then
            power=a+b*grid(L)
            fluxmult(L)=((now-timega)**power)/(now**power)
         else
            fluxmult(L)=1.0
         endif

!then do the same for chromofit using 'diff' as the indexing offset 
         if (grid(L).ge.1500.0 .and. grid(L).le.2200.0) then   ! eqn 4 in paper
            chromofit(L)=-41.4287 + 0.665466*nm(L) - 0.00360537*
     $      nm(L)**2 + 6.54727E-6*nm(L)**3
            if (chromofit(L).gt.0.0) then  !check this...
               chromofit(L)=0.0
            endif
         endif
!now set chromospheric "multipliers"
         chromomult(L)=(((now-timega)**chromofit(L)) /
     $                   (now**chromofit(L)))-1.
      enddo

      time=now-timega

      call spline(BahcallAge,RnormBahcall,41,1.d40,1.d40,temp)
      call splint(BahcallAge,RnormBahcall,temp,41,time,Rnorm)
      call spline(BahcallAge,LnormBahcall,41,1.d40,1.d40,temp)
      call splint(BahcallAge,LnormBahcall,temp,41,time,Lnorm)

      RnormNOW=1.0
      LnormNOW=1.0

      Tnow=5777.  !photospheric temperature of the modern sun (kelvin)
      Tphoto=Tnow*(Lnorm/Rnorm**2)**0.25

c     OK, main code to put it all together starts here...

c     1)find the relative flux

c file from ~/youngsun/kurucz/fluxreadAGE.pro 
      open(65, file='PHOTOCHEM/DATA/FLUX/kuruczflux.dat',
     &         status='UNKNOWN')  
!note the .dat file has fluxes in W/m^2/ster converted from
!the frequency units output by kurucz models
      
 111  FORMAT(44(E13.6))
 112  FORMAT(43(F6.1,2X))

c     Read Kurucz fluxes into a dummy variable because .dat file is weird and has wavelengths at the end of each row of 43 fluxes assosicated with each age
      do i=1,1221
         READ(65,111)dummy(1:44,i)
      enddo
      READ(65,112)AGE

c     1) Assign wave and kuruczflux to approrpriate elements of dummy
      do j=1,1221
         do i=1,44
            if (i.eq.44) then
                wave(j)=dummy(i,j)
                !again, wavelengths are the 44th element of dummy
             else
                kuruczflux(i,j)=dummy(i,j) !every other element is a kuruczflux
             endif
          enddo
       enddo

c     Subroutine fails if time is outside the photospheric models computed by Claire et al. 2012
       if ((time.le.minval(AGE)).or.(time.ge.maxval(AGE))) then
          print * , 'Stopping in Youngsun.f'
          print *, 'NEED TO COMPUTE MORE KURUCZ FLUX MODELS'
          !the solar flux model is only applicable from >0.0 to 8.0 Gyrs.
          STOP
      endif

c     2) find where "time" is within the grid of Kurucz model ages
      joe1=minloc(AGE,1,AGE.ge.time)-1
      Alow=AGE(joe1)
      Ahigh=AGE(joe1+1)

c     3) compute weights for this age
      w2=(time-Alow)/(Ahigh-Alow)
      w1=1.0-w2

c     4) do the same for modern temperature
      joe2=minloc(AGE,1,AGE.ge.now)-1
      Alownow=AGE(joe2)
      Ahighnow=AGE(joe2+1)

      w2now=(now-Alownow)/(Ahighnow-Alownow)
      w1now=1.0-w2now

c     5) compute fluxes by linear interpolation between solar ages
      photonow=(kuruczflux(joe2,1:1221)*w1now)+(kuruczflux(joe2+1,1:1221
     $)*w2now)
      photothen=(kuruczflux(joe1,1:1221)*w1)+(kuruczflux(joe1+1,1:1221)*
     $w2)


c     6) These are fluxes at the solar surface. Now compute fluxes at earth
! technically, the fluxes are per steradian so should be multiplied by 4pi but we
! normalize soon so don't bother...

c     need to multiply by R_sun/R_earth
      Rsunnow=6.96e8       !solar radius in meters
      Rstar=Rnorm*Rsunnow  !solar radii at times
      Rearth=1.49598e11    !1AU - assuming earth always at 1AU 
      fluxfactor=(Rstar/Rearth)**2  !flux factor at times vector

      fluxfactornow=(Rsunnow/Rearth)**2

      photothen=photothen*fluxfactor
      photonow=photonow*fluxfactornow


c     7) compute relative flux at earth on kurucz grid
      relphotokuruczgrid(100:1221)=photothen(100:1221)/
     $photonow(100:1221)
c     wave[100]=69.5 nm, and photonow is 0 below this. In the final code,
c     we only consider the photospheric contribution above 150nm anway...


c     8) interpolate kurucz grid to Thuillier wavelength grid
      wave(1:1221)=wave(1:1221)*10.0 !convert kurucz wavelength grid to Angstroms
      call inter1(n,grid,relphoto,1221,wave,relphotokuruczgrid)
      !inter1 is grid to grid interpolation

c     9) add in chromospheric excess to photosphere
      do L=1,n
         if (grid(L).LE.2200.0 .AND. grid(L).GE.1500.) then
            relphoto(L)=relphoto(L)+chromomult(L) 
         endif
      enddo
c     so relphoto now contains the kurucz + chromosphere
   
      start=1500.0 !chromospheric weight interval #1 (lower bound in Angstroms)
      finish=1700.0 !chromspheric weight interval #2 (upper bound in Angstroms)

      del=finish-start
      

c     10)start the new chromoshere process, using the start/finish weights
c - for completeness, also add in the shortest wavelength Ribas beta
c    (although we usually don't go much shorter than Ly a....)
      do L=1,n

         if (grid(L).GE.start.AND.grid(L).LT.finish) then
            !compute weight
            w2=(grid(L)-start)/del
            w1=1.0-w2
            !scale flux with weights
            fluxmult(L)=(fluxmult(L)*w1)+(relphoto(L)*w2)
         endif 

         if (grid(L).GE.finish) then
            fluxmult(L)=relphoto(L)
         endif 

         if (grid(L).LE.20.) then
            fluxmult(L)= ((now-timega)**(-1.92))/(now**(-1.92)) 
!applying ribas 1-20 A beta value for 2 ATLAS points less than 25 A - revisit if higher energy solar flux is ever used.
         endif

      enddo

c     11)add in strong lines from Ribas et al. 2005

!some issue with 285., 305, 365,975,1025,1035
      stronglines=(/284.,304.,361.,977.,1026.,1032.,1038.,1176.,1206., 
     $             1216.,1304.,1335.,1400.,1550.,1640.,1657./)
c     slinebetas=strong line betas (beta means power law slope)
      slinebetas=(/ -1.79,-1.34,-1.86,-0.85,-1.24,-1.00,-1.02,-1.02,
     $              -0.94,-0.72,-0.78,-0.78,-0.97,-1.08,-1.28,-0.68/)

c     Determine width of line, then calculate relative flux in that bin
      do k=1,16 
         if (stronglines(k).LT.1200.) then
            width=5.0  !these are tied to the use of the ATLAS 1 grid...
         else
            width=1.0
            if (stronglines(k).EQ.1216.) width=8.0 !lyman alpha width of 16A
         endif
        lowerl=stronglines(k)-width
        upperl=stronglines(k)+width

         do j=1,n
            if ((grid(j).GE.lowerl).AND.(grid(j).LE.upperl)) then
               fluxmult(j)=((now-timega)**slinebetas(k)) / 
     $                     (now**slinebetas(k))
            endif
         enddo
      enddo

c      do j=1,n
c         print *, grid(j),fluxmult(j)
c      enddo
c      stop



      endif  
c     END OPTION 2
      
      
      

      RETURN
      END
