!****************************************************************************
! Reference Selsis et al. (2007) A &A, 476, 1373
!****************************************************************************
implicit none
real *8 lin_sun,lout_sun,a_in,a_out,b_in,b_out,L_sun,pi
real *8 Tstar,Teff,lin,lout,Lum,logmass,radius,radsun,sigma
real *8 mass

!open(9,file='MSstardata.dat')
!open(10,file='inoutHZ.dat')
open(9,file='Selsis_50pct.dat')

pi     = acos(-1.0d0)
radsun = 695500000.0d0
sigma  = 5.672*1.0d-8
L_sun  = 3.845*1.0d+26
a_in  = 2.7619E-05
a_out = 1.3786E-04
b_in  = 3.8095E-09
b_out = 1.4286E-09

!*******************************
! 0% Cloud cover
!lin_sun  = 0.95
!lout_sun = 1.67

!*******************************
! 50% Cloud cover
lin_sun  = 0.76
lout_sun = 1.95

!*******************************
! 100% Cloud cover
!lin_sun  = 0.51
!lout_sun = 2.40


!do while(eof(9)==.false.)
!read(9,*)logmass,Lum,radius


!mass = 7.00d0
!do while(mass >=0.060d0)
!  logmass = log10(mass)
!************************************************************************
! These relations are for M-dwarf stars in the mass range 0.5 - 0.06 Msun
!  Lum = -1.01*logmass**3 -1.739*logmass**2 + 1.269*logmass-0.6154
!  radius = 0.5444*logmass**3 + 1.309*logmass**2 +1.765*logmass+0.196
!************************************************************************
!  Lum = -0.3827*logmass**4 -0.6862*logmass**3 + 0.8617*logmass**2&
!        +3.747*logmass -0.01675
!  radius = -0.03133*logmass**3 - 0.4221*logmass**2 +0.8194*logmass+0.009744
!  radius = 10**radius
!  Teff = ((10**Lum)*L_sun/(4.0d0*pi*(radius*radsun)**2*sigma ))**(0.25d0)

Teff = 3700.0d0
do while(Teff <=7201.0d0)
! radius = 0.211*radsun
!  Lum = sigma*Teff**4*4.0d0*pi*radius**2
!  Lum = Lum/L_sun
!   Lum = 0.5d0   
!  Lum = 10**Lum

if(Teff < 3700.0d0)then
   Teff = 3700.0d0
endif
Tstar = Teff-5700.0d0
lin = (lin_sun - a_in*Tstar - b_in*Tstar**2)*sqrt(Lum)
lout = (lout_sun - a_out*Tstar - b_out*Tstar**2)*sqrt(Lum)
write(9,*)Teff,1.0/(lin_sun - a_in*Tstar - b_in*Tstar**2)**2,&
        1.0/(lout_sun - a_out*Tstar - b_out*Tstar**2)**2
!print *,lin,lout
!print *,1.0/(lin_sun - a_in*Tstar - b_in*Tstar**2)**2,&
!        1.0/(lout_sun - a_out*Tstar - b_out*Tstar**2)**2
Teff = Teff + 200.0d0
enddo
!  mass = mass - 0.01d0
!if ( (Teff >3699.0d0) .and. (Teff <7201.0d0) )then
!write(10,100)mass,lin,lout,log10(Lum),Teff
!endif
!enddo
100 format(1p6e12.3)
close(9)
close(10)
end
