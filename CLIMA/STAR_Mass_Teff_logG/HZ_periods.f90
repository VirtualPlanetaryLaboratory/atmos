!**********************************************************************************
! This code calculates HZ distances and periods for stars in the range
! of 2600 K < Teff < 7200 K. It also gives the stellar fluxes incident
! on the top of the atmosphere of a planet at the HZ boundary limits.
! (Basically HZ fluxes).
!
!The Teff-Lum relationships for these stars were obtained in the following way:
!
! The files 'dartmouth.txt' and 'read_dartmouth.f90' were used to 
! generate  'star_lum.dat' which has temperatures and luminosities for 
! 3500 K < Teff < 6200 K.
!
! Dartmouth paper reference:
! http://adsabs.harvard.edu/abs/2008ApJS..178...89D
!
!
! The files 'baraffe.txt' and 'read_baraffe.f90' were used to generate
! 'star_lum_baraffe.dat' which has temperatures and luminosities for
! 2000 K < Teff < 3500 K.
!
! Baraffe reference:
! BARAFFE, CHABRIER, ALLARD, HAUSCHILDT, 1998, A&A 337, 403
!
!
! Both of these data sets were combined to make a file with which has
! luminosities in the range 2600 K < Teff < 6200 K. Then, I did a curve
! fit for Teff-Lum in this temperature range to derive a Teff-Lum
! relationship ('teff_logL.dat').
!
! A similar analysis was done for Teff-mass relationship ('teff_mass.dat')
!
! The files 'IHZ_distance.dat' and 'OHZ_distance.dat' were generated 
! combining the files 'IHZ_stars.dat',  'OHZ_stars.dat' and
! 'teff_logL.dat'.

!*********************************************************************************

implicit none

real *8 teff,mass,distance_ihz,distance_ohz,lum,period,junk
real *8 seff_ihz,seff_ohz,period_ihz,period_ohz
Real *8 svenus,smars,s2AU,disvenus,dismars,periodvenus,periodmars
real *8 dis2AU,period2AU
real *8 G,au,msun,day,pi

AU     = 1.49598*1.0d+11
G      = 6.67*1.0d-11
msun   = 2.0d0*1.0d+30
day    = 24.0d0*3600.0d0 
pi     = acos(-1.0d0)

open(9,file='IHZ_distance.dat')
open(10,file='OHZ_distance.dat')
open(11,file='total_teff_mass.dat')
open(12,file='HZ_periods.dat')

write(12,'(a100)')trim('# I fudged luminosities a *little* for stars >=6400 K to have a good looking HZ distance plot')
write(12,'(a100)')trim('# This will in no way affect the Seff fluxes as they are independent of the luminosities')
write(12,200)
200 format(5x,'TEFF (K)',7x,'Mass',10x,'Lum',8x,'IHZ(AU)',8x,'OHZ(AU)',6x,&
           'P_IHZ(days)',3x,'P_OHZ(days)',3x,'Flux_IHZ',6x,'FLUX_OHZ',2x,&
           'd_Venus(AU)',2x,'d_Mars',2x,'Pvenus',2x,'Pmars',2x,'flux_V',2x,'flux_M',&
           2x,'flux_2AU')

do while (eof(9)==.false.)
   read(9,*)teff,seff_ihz,junk,lum
   read(10,*)junk,seff_ohz
   read(11,*)junk,mass
   svenus = seff_ihz*1.777/1.015   ! (1.777/1.015) is the ratio of the recent Venus flux 
                                   ! to waterloss flux for the Sun

   smars  = seff_ohz*0.32d0/0.350   !(0.3/0.346) is the ratio of the Early Mars flux
                                    ! to maximum greenhouse flux for the Sun
   
   s2AU   = seff_ohz*0.25d0/0.350   ! (0.25/0.346) is the ratio of the flux at 2 AU
                                    ! to the maximum greenhouse flux for the Sun


   disvenus = sqrt(lum/svenus)
   dismars  = sqrt(lum/smars)
   dis2AU   = sqrt(lum/s2AU)
   distance_ihz = sqrt(lum/seff_ihz)
   distance_ohz = sqrt(lum/seff_ohz)
   period_ihz  = 2.0d0*pi*(sqrt((distance_ihz*AU)**3/(G*mass*msun)))/day
   period_ohz  = 2.0d0*pi*(sqrt((distance_ohz*AU)**3/(G*mass*msun)))/day
   
   periodvenus  = 2.0d0*pi*(sqrt((disvenus*AU)**3/(G*mass*msun)))/day
   periodmars   = 2.0d0*pi*(sqrt((dismars*AU)**3/(G*mass*msun)))/day
   period2AU    = 2.0d0*pi*(sqrt((dis2AU*AU)**3/(G*mass*msun)))/day
   write(12,100)teff,mass,lum,distance_ihz,distance_ohz,period_ihz,&
                period_ohz,seff_ihz,seff_ohz,disvenus,dismars,periodvenus,&
                periodmars,svenus,smars,s2AU

enddo



100 format(1p17e14.5)

close(9)
close(10)
close(11)
close(12)

end
