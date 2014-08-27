implicit none

real *8 rt,P0,t,Q,Mstar,msun,AU
open(9,file='tidal_lock.dat')

msun = 1.99e+33
Mstar = 0.1d0
P0 = 13.5*3600.0d0
t  = 1.0e+9*31536000.0d0 *2
Q = 100.0d0
AU = 1.4959787e+13

do while(Mstar <=1.4d0)
   rt = 0.027*(P0*t/Q)**(1.0d0/6.0d0)*(Mstar*msun)**(1.0d0/3.0d0)
   rt = rt/AU
   write(9,*)rt,Mstar
   Mstar = Mstar + 0.05d0
enddo

end
