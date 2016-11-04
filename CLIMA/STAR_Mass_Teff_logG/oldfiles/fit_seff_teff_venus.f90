!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ihzflux  

 a= 0.5423625781192452D0 
 b= 0.000777762140255013D0 
 c= -2.971663454325353D-07 
 d= 5.041964558021062D-11 
 e= -2.762072807929421D-15 
 teff = 5780.0d0

 constant = a + b*teff + c*teff**2 + d*teff**3 + e*teff**4
 order1   = b + 2.0*c*teff + 3.0*d*teff**2 + 4.0*e*teff**3
 order2   = c + 3.0*d*teff + 6.0*e*teff**2
 order3   = d + 4.0*e*teff
 order4   = e

 print *,constant,order1,order2,order3,order4

  x = 2600.0
  do while(x <=7201.d0)
    y=eqn6001(x)
    ihzflux = constant + order1*(x-teff) + order2*(x-teff)**2 + &
               order3*(x-teff)**3 + order4*(x-teff)**4
     print *,x,ihzflux
    x = x + 200.0
  enddo

20 END
  
!-----------------------------------------------------------
REAL*8 FUNCTION eqn6001(x)
!-----------------------------------------------------------
! TableCurve C:\Users\winey-akka\Desktop\fit_seff_teff_venus.f90 Oct 23, 2012 2:03:56 PM 
! C:\Users\winey-akka\Desktop\HZ_periods.dat 
! X= 2.60000E+03 
! Y= 1.31726E+00 
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.999932657296582D0 
! r2adj=0.999913950990077D0 
! StdErr=0.002414890372811472D0 
! Fval=70529.98886425333D0 
! a= 0.5423625781192452D0 
! b= 0.000777762140255013D0 
! c= -2.971663454325353D-07 
! d= 5.041964558021062D-11 
! e= -2.762072807929421D-15 
!-----------------------------------------------------------
  REAL*8 x,y
  y=0.5423625781192452D0+x*(0.0007777621402550130D0+&
   &x*(-2.971663454325353D-07+x*(5.041964558021062D-11+&
   &x*(-2.762072807929421D-15)))) 
  eqn6001=y
  RETURN
END
  
