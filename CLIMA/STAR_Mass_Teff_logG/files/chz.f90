implicit none

real *8 time,teff,junk,logl
real *8 ihz,ohz
real *8 seff(6),seffsun(6),a(6),b(6),c(6),d(6),tstar


open(9,file='starLum_Time_0.3Msun.dat')
open(10,file='chz_0.3Msun.dat')

seffsun(1) = 1.7763
seffsun(2) = 1.0385
seffsun(3) = 1.0146
seffsun(4) = 0.3507
seffsun(5) = 0.2946
!seffsun(6) = 0.2484

a(1) = 1.4335e-4
a(2) = 1.2456e-4
a(3) = 8.1884e-5
a(4) = 5.9578e-5
a(5) = 4.9952e-5
!a(6) = 4.2588e-5

b(1) = 3.3954e-9
b(2) = 1.4612e-8
b(3) = 1.9394e-9
b(4) = 1.6707e-9
b(5) = 1.3893e-9
!b(6) = 1.1963e-9

c(1) = -7.6364e-12
c(2) = -7.6345e-12
c(3) = -4.3618e-12
c(4) = -3.0058e-12
c(5) = -2.5331e-12
!c(6) = -2.1709e-12

d(1) = -1.1950e-15
d(2) = -1.7511E-15
d(3) = -6.8260e-16
d(4) = -5.1925e-16
d(5) = -4.3896e-16
!d(6) = -3.8282e-16




do while(eof(9)==.false.)
  read(9,*)junk,time,teff,junk,logl
   tstar = teff - 5780.0d0
  seff(3) = seffsun(3) + a(3)*tstar + b(3)*tstar**2 + c(3)*tstar**3 + d(3)*tstar**4  
  seff(4) = seffsun(4) + a(4)*tstar + b(4)*tstar**2 + c(4)*tstar**3 + d(4)*tstar**4  
  ihz = sqrt(10.**logl/seff(3))
  ohz = sqrt(10.**logl/seff(4))
  write(10,100)time,ihz,ohz,logl
100 format(1p6e14.5)

enddo

end
