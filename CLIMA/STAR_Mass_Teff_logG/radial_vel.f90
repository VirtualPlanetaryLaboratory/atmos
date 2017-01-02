implicit none
real *8 Kvel,mp,mstar,msun,mjup,period,incl,ecc,oneyear
real *8 pi,semi,AU,G,mearth,minmass,teff,x,eqn7505,eqn6005
real *8 lum,seff
open(9,file='radial_vel')

msun    = 2.0d0*1.0d+30
mjup    = 1.8987d0 *1.0d+27
mearth  = 5.9742*1.0d+24
!mp      = 0.0031464686d0 !*1.0d+27 Kg
!mp      = 1.0 *mearth
AU     = 1.49598*1.0d+11
G      = 6.67*1.0d-11

oneyear = 31536000.0d0
ecc     = 0.55d0
!period  = 1.0d0/365.0d0
!period  = 1353.6d0


semi = 1.45d0
pi      = acos(-1.0d0)
incl    = pi/2.0d0


minmass = 0.9d0
mstar   = 1.70d0
teff    = 3000.0d0 
!do while(teff <=5780.0d0)
!   x    = teff
!   mstar=eqn7505(x)
!   lum  = eqn6005(x)
!   lum  = 10.0**(lum)
    


!!        semi = 0.01d0

!        do while(semi <=1.0d0)
!           seff  = lum/semi**2
           period =  2.0d0*pi*sqrt((semi*AU)**3/(G*mstar*msun))
           Kvel = ((2.0d0*pi*G)**(1.0d0/3.0d0)/sqrt(1.0d0 - ecc**2))&
                 *(minmass*mjup) &
!                *(mp*sin(incl)) &
                 *(period)**(-1.0d0/3.0d0)*(mstar*msun)**(-2.0d0/3.0d0)
!          write(9,100)mstar,semi,Kvel
    
!          write(9,100)teff,seff,Kvel,period,semi,mstar
!          print *,teff,seff,Kvel,period/31536000.0,semi,mstar
          print *,Kvel,period/31536000.0,semi,mstar
!          semi = semi + 0.01d0
!        enddo

!        write(9,*)

!mstar = mstar + 0.1d0
!teff = teff + 100.0d0
!enddo
100 format(1p6e14.5) 
end


!-----------------------------------------------------------
REAL*8 FUNCTION eqn7505(x)
!-----------------------------------------------------------
! TableCurve C:\Users\winey-akka\Desktop\teff_mass.f90 Sep 19, 2012 6:59:12 PM 
! C:\Users\winey-akka\Desktop\combined_lum.dat 
! X= 2.00900E+03 
! Y= 7.50000E-02 
! Eqn# 7505  y=(a+cx^(0.5)+ex+gx^(1.5))/(1+bx^(0.5)+dx+fx^(1.5)) 
! r2=0.9994766832916201D0 
! r2adj=0.9994340877455891D0 
! StdErr=0.00753892057263051D0 
! Fval=27693.3865776872D0 
! a= -0.1900924095182077D0 
! b= -0.04290251362241933D0 
! c= 0.01239091923915334D0 
! d= 0.0005926755682815886D0 
! e= -0.0002609242006665851D0 
! f= -2.574529567515367D-06 
! g= 1.793738550364907D-06 
!-----------------------------------------------------------
REAL*8 x,y
  x=DSQRT(x)
  y=(-0.1900924095182077D0+x*(0.01239091923915334D0+&
   &x*(-0.0002609242006665851D0+x*(1.793738550364907D-06))))/&
   &(1.0+x*(-0.04290251362241933D0+x*(0.0005926755682815886D0+&
   &x*(-2.574529567515367D-06))))
  x= x**2
  eqn7505=y
  RETURN
END


!-----------------------------------------------------------
REAL*8 FUNCTION eqn6005(x)
!-----------------------------------------------------------
! TableCurve C:\Users\winey-akka\Desktop\teff_logL.f90 Sep 19, 2012 6:33:59 PM 
! C:\Users\winey-akka\Desktop\combined_lum.dat 
! X= 2.00900E+03 
! Y= -3.91273E+00 
! Eqn# 6005  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8 
! r2=0.9990477452211798D0 
! r2adj=0.998945717923449D0 
! StdErr=0.0330916163241971D0 
! Fval=11147.10320081128D0 
! a= -613.8590440214133D0 
! b= 1.265608220688683D0 
! c= -0.001111597221651978D0 
! d= 5.398863926836034D-07 
! e= -1.588171015817137D-10 
! f= 2.906250529922859D-14 
! g= -3.240702181938055D-18 
! h= 2.018811990555861D-22 
! i= -5.391899214613087D-27 
!-----------------------------------------------------------
  REAL*8 x,y
  y=-613.8590440214133D0+x*(1.265608220688683D0+&
   &x*(-0.001111597221651978D0+x*(5.398863926836034D-07+&
   &x*(-1.588171015817137D-10+x*(2.906250529922859D-14+&
   &x*(-3.240702181938055D-18+x*(2.018811990555861D-22+&
   &x*(-5.391899214613087D-27))))))))
  eqn6005=y
  RETURN
END

