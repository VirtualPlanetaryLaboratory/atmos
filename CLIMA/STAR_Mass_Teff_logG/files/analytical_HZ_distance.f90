!************************************************************************************
! This code calculates 'continuous habitable zone' (CHZ) distances using stellar 
! evolutionary models from Baraffe et al.(1998). 
!
! Note that stars with Teff > 6000 K (Mass > 1.15 MSun) are not in the Baraffe data 
!file because while  calculating CHZ, a 5 Gyr time is assumed & high mass stars 
! do not survive that long.
!
! Ravi kumar Kopparapu Feb 25 2012
!************************************************************************************
implicit none
real *8 seff(7),seffsun(7),teff,a(7),b(7),c(7),d(7),tstar
real *8 mass,time,logg,logl,distance(7,26)
real *8 mgmg,venusmars
integer i,k,flag

!************************************************************************************
! 'HZ_distances.dat' is the output file
! 'baraffe_mass_time.dat' is an input file which contains stellar data needed for this
! program. See that file for description of the data.
open(9,file='HZ_distances.dat')
open(10,file='baraffe_mass_time_New.dat')


!************************************************************************************
! Coeffcients to be used in the analytical expression to calculate habitable zone flux 
! boundaries


seffsun(1) = 1.7763
seffsun(2) = 1.0385
seffsun(3) = 1.0146
seffsun(4) = 0.3507
seffsun(5) = 0.3207
seffsun(6) = 0.2484
seffsun(7) = 0.5408

a(1) = 1.4335e-4
a(2) = 1.2456e-4
a(3) = 8.1884e-5
a(4) = 5.9578e-5
a(5) = 5.4471e-5
a(6) = 4.2588e-5
a(7) = 4.4499e-5

b(1) = 3.3954e-9
b(2) = 1.4612e-8
b(3) = 1.9394e-9
b(4) = 1.6707e-9
b(5) = 1.5275e-9
b(6) = 1.1963e-9
b(7) = 1.4065e-10

c(1) = -7.6364e-12
c(2) = -7.6345e-12
c(3) = -4.3618e-12
c(4) = -3.0058e-12
c(5) = -2.7481e-12
c(6) = -2.1709e-12
c(7) = -2.2750e-12

d(1) = -1.1950e-15
d(2) = -1.7511E-15
d(3) = -6.8260e-16
d(4) = -5.1925e-16
d(5) = -4.7474e-16
d(6) = -3.8282e-16
d(7) = -3.3509e-16

print *,'What do you want to calculate?'
print *,'For Continuous HZ (0 Gyr to 5 Gyr), press 1'
print *,'For HZ at Zero Age Main-Sequence (ZAMS, at 0 Gyr), press 2'
read(*,*)flag



write(9,90)
90 format('#', 4x,'Mass',10x,'Teff(K)',8x,'Recent',7x,'Runaway',7x,'Moist',8x,'Maximum',7x,&
          'Early')
write(9,91)
91 format('#',33x,'Venus',8x,'Greenhouse',4x,'Greenhouse',3x,'Greenhouse',4x,'Mars')

!************************************************************************************
! Skip comments in the input data file
do i = 1,19
 read(10,*)
enddo

k=1

!************************************************************************************
! Read stellar parameters from 'baraffe_mass_time.dat' file

do while(eof(10)==.false.)
   read(10,*)mass,time,teff,logg,logl
   tstar = teff - 5780.0d0

!************************************************************************************
! Calculate HZ fluxes first & corresponding distances
   do i = 1,7
     seff(i) = seffsun(i) + a(i)*tstar + b(i)*tstar**2 + c(i)*tstar**3 + d(i)*tstar**4
     distance(i,k) = sqrt(10.**logl/seff(i))
   enddo

!************************************************************************************
! Identify Continuous Habitable Zone (CHZ) for various limits. The odd number values of
! 'k' are for 0.5 Gyr, the even number values of 'k' are for 5 Gyr.
!
! i = 1 --> Recent Venus
! i = 2 --> Runaway Greenhouse
! i = 3 --> Moist Greenhouse
! i = 4 --> Maximum Greenhouse
! i = 5 --> Early Mars
! i = 6 --> 2 AU cloud limit
! i = 7 --> First CO2 condensation
!


! Write the data only 'k' is even number. That means, at 5 Gyr to calculate 
! 'continuous HZs'
if(flag==1)then
   if(mod(k,2)==0)then
     write(9,100)mass,teff,distance(1,k),distance(2,k),distance(3,k),distance(4,k-1),distance(5,k)
     print *,k
   endif
endif

! Write the data only 'k' is odd number. That means, at 0 Gyr to calculate 
! 'Zero Age Main-Sequence (ZAMS) HZs'
if(flag==2)then
   if(mod(k,2) .ne. 0)then

     write(9,100)mass,teff,distance(1,k),distance(2,k),distance(3,k),distance(4,k),&
                 distance(5,k)
   endif
endif

   k = k + 1
enddo  ! ENDDO for reading 'baraffe_mass_time.dat' file


100 format(1p7e14.5)
200 format(1p7e14.5)
close(9)
close(10)
end
