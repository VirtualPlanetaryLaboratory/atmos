implicit none
real *8 mass,logteff,logg,junk,logl
real *8 teff,radius,gg,Gconst,radsun,sigma,lum
integer i

!open(9,file='dartmouth.txt')
open(9,file='Dartmouth_Highmass.dat')
open(10,file='star_lum_highmass.dat')

Gconst = 6.67e-11
radsun = 695500000.0d0
sigma  = 5.672*1.0d-8

do i = 1,9
read(9,*)
enddo

do while(eof(9)==.false.)
   read(9,*)junk,mass,logteff,logg,logl
   teff = 10.0**logteff
   gg = 10.0**logg           ! cm/s^2
   gg = gg/10**(4.44)        ! Normalize to Sun's gravity of 4.44 cm/s^2
   radius = sqrt(mass/(gg))  ! Radius in Solar radius (mass is in solar mass)
   lum = (teff/5780.0)**4 * radius**2
   write(10,100)teff,log10(lum),logl,radius,mass

enddo
100 format(1p6e14.5)
close(9)
close(10)

end
