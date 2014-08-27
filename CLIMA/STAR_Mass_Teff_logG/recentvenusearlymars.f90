implicit none

integer i
real *8 teff,seffi,seffo,svenus,smars,junk,lum,distance
open(9,file='IHZ_distance.dat')
open(10,file='OHZ_distance.dat')
open(11,file='recentvenusearlymars.dat')

do i = 1,24
 read(9,*)teff,seffi,junk,lum
 read(10,*)junk,seffo
 svenus = (seffi*1.777/1.015)
 smars  = seffo*0.32d0/0.350
 distance = sqrt(lum/svenus)
 write(11,*)teff,svenus,smars
 print *,teff,svenus,distance
 
enddo



close(9)
close(10)
close(11)
end
