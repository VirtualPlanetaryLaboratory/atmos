!****************************************************************************
! Reference Selsis et al. (2007) A &A, 476, 1373
!****************************************************************************
implicit none
real *8 lin_sun,lout_sun,a_in,a_out,b_in,b_out,L_sun,pi
real *8 Tstar,Teff,lin,lout,Lum,logmass,radius,radsun,sigma
real *8 mass,eqn6005,tt,time,logg
integer i,k

!open(9,file='MSstardata.dat')
!open(9,file='total_teff_lum.dat')
open(9,file='files/baraffe_mass_time.dat')
open(10,file='SelsisHZ0cloud_distance.dat')
!open(11,file='SelsisHZ0cloud_seff.dat')

pi     = acos(-1.0d0)
radsun = 695500000.0d0
sigma  = 5.672*1.0d-8
L_sun  = 3.845*1.0d+26
a_in  = 2.7619E-05
a_out = 1.3786E-04
b_in  = 3.8095E-09
b_out = 1.4286E-09


!****************************************************************************
! For 0% cloud cover
lin_sun  = 0.95
lout_sun = 1.67

!****************************************************************************
! For Venus & Mars critical.
!lin_sun  = 0.72
!lout_sun = 1.77


do i = 1,19
read(9,*)
enddo

tt = 2600.0

k = 1
do while(eof(9)==.false.)
   read(9,*)mass,time,Teff,logg,Lum
   print *,mass
   if(Teff < 3700.0d0)then
      Teff = 3700.0d0
   endif
   Tstar = Teff-5700.0d0
   Lum = 10.**Lum


   lin = (lin_sun - a_in*Tstar - b_in*Tstar**2)*sqrt(Lum)
   lout = (lout_sun - a_out*Tstar - b_out*Tstar**2)*sqrt(Lum)
!  print *,lin,lout,Teff
   if(mod(k,2)==0)then
      print *,1.0/(lin_sun - a_in*Tstar - b_in*Tstar**2)**2,&
              1.0/(lout_sun - a_out*Tstar - b_out*Tstar**2)**2,tt
      write(10,100)lin,lout,mass,teff,lum
   endif
   k = k + 1
!write(11,100)1.0/(lin_sun - a_in*Tstar - b_in*Tstar**2)**2,&
!             1.0/(lout_sun - a_out*Tstar - b_out*Tstar**2)**2,&
!             tt

!  mass = mass - 0.01d0
!tt = tt + 200.0
enddo
100 format(1p6e12.3)
close(9)
close(10)
!close(11)
end

!-----------------------------------------------------------
REAL*8 FUNCTION eqn6005(x)
!-----------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! APPLICABLE ONLY FOR TEMPERATURES 2100 - 6200 K. But I am also using
! it until 7200 K
!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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





!-----------------------------------------------------------
!REAL*8 FUNCTION eqn6003(x)
!-----------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! APPLICABLE ONLY FOR TEMPERATURES BELOW 4000 K
!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TableCurve C:\Users\winey-akka\Desktop\eqn6003-rank76.f90 May 18, 2012 8:49:14 PM 
! C:\Users\winey-akka\Desktop\b3.    
! X=  
! Y=  
! Eqn# 6003  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6 
! r2=0.9998617913718321D0 
! r2adj=0.9998108724035597D0 
! StdErr=0.008419563419147614D0 
! Fval=24114.79429385633D0 
! a= -337.5466889727813D0 
! b= 0.1113123704147002D0 
! c= 0.0003140864810252149D0 
! d= -2.944835680841914D-07 
! e= 1.058468336681171D-10 
! f= -1.754824936838729D-14 
! g= 1.119603420214092D-18 
!-----------------------------------------------------------
!  REAL*8 x,y
!  y=-337.5466889727813D0+x*(0.1113123704147002D0+&
!   &x*(0.0003140864810252149D0+x*(-2.944835680841914D-07+&
!   &x*(1.058468336681171D-10+x*(-1.754824936838729D-14+&
!   &x*(1.119603420214092D-18))))))
!  eqn6003=y
!  RETURN
!END

