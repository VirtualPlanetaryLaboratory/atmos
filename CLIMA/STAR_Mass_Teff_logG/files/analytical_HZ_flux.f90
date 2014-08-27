!************************************************************************************
! This code calculates habitable zone 'fluxes' using the expression given in the 
! Kopparapu et al. paper. The corresponding output file is 'HZ_fluxes.dat'. 
! It also generates a file 'HZ_coefficients.dat' that gives the coefficients for 
! the analytical expression.
!
! Ravi kumar Kopparapu Feb 25 2012
!************************************************************************************
implicit none
real *8 seff(7),seffsun(7),teff,a(7),b(7),c(7),d(7),tstar
integer i

!************************************************************************************
! Output files.

open(9,file='HZ_fluxes.dat')
open(10,file='HZ_coefficients.dat')

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

!************************************************************************************
! Writing coefficients into 'HZ_coefficients.dat' file

write(10,*)'# The coefficients are as follows. The columns, i, are arranged according to'
write(10,*)'# the HZ limits given in the paper.'
write(10,*)'#'
write(10,*)'# i = 1 --> Recent Venus'
write(10,*)'# i = 2 --> Runaway Greenhouse'
write(10,*)'# i = 3 --> Moist Greenhouse'
write(10,*)'# i = 4 --> Maximum Greenhouse'
write(10,*)'# i = 5 --> Early Mars'
write(10,*)'# i = 6 --> 2 AU cloud limit (Note: This limit is not given in the paper)'
write(10,*)'# i = 7 --> 1st CO2 condensation limit (Note: This limit is not given in the paper)'

write(10,*)'# First row: S_effSun(i)'
write(10,*)'# Second row: a(i)'
write(10,*)'# Third row:  b(i)'
write(10,*)'# Fourth row: c(i)'
write(10,*)'# Fifth row:  d(i)'
write(10,200)(seffsun(i),i=1,7)
write(10,200)(a(i),i=1,7)
write(10,200)(b(i),i=1,7)
write(10,200)(c(i),i=1,7)
write(10,200)(d(i),i=1,7)

!************************************************************************************
! Calculating HZ fluxes for stars with 2600 K < T_eff < 7200 K. The output file is
! 'HZ_fluxes.dat'

teff  = 2600.0d0
write(9,90)
90 format('#', 2x,'Teff(K)',8x,'Recent',8x,'Runaway',7x,'Moist',&
          9x,'Maximum',7x,'Early',10x,'2 AU',9x,'First CO2')
write(9,91)
91 format('#',17x,'Venus',9x,'Greenhouse',4x,'Greenhouse',4x,'Greenhouse'&
          ,4x,'Mars',8x,'Cloud limit',4x,'Condensation')

do while(teff <=7201.0d0)
  tstar = teff - 5780.0d0

  do i = 1,7
     seff(i) = seffsun(i) + a(i)*tstar + b(i)*tstar**2 + c(i)*tstar**3 + d(i)*tstar**4
     print *,seff(i),teff
  enddo
  write(9,100)teff,(seff(i),i=1,7)
  print *,''
  teff = teff + 200.0d0
enddo


100 format(1p8e14.5)
200 format(1p8e14.5)
close(9)
close(10)
end
