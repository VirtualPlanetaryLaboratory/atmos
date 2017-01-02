!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ihzflux  

 a= 0.3097398289189322D0 
 b= 0.0004442895257152834D0 
 c= -1.697501412889046D-07 
 d= 2.880070856747148D-11 
 e= -1.577741778149674D-15 
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
!    PRINT *,x,y
     ihzflux = constant + order1*(x-teff) + order2*(x-teff)**2 + &
               order3*(x-teff)**3 + order4*(x-teff)**4
     print *,x,ihzflux
    x = x + 200.0
  enddo
20 END
  
!-----------------------------------------------------------
REAL*8 FUNCTION eqn6001(x)
!-----------------------------------------------------------
  REAL*8 x,y
! TableCurve C:\Users\winey-akka\Desktop\fit_seff_teff_ihz.f90 Oct 22, 2012 8:40:19 PM 
! C:\Users\winey-akka\Desktop\HZ_periods.dat 
! X= 2.60000E+03 
! Y= 7.52400E-01 
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9999327238118487D0 
! r2adj=0.9999140359818067D0 
! StdErr=0.001378672163286898D0 
! Fval=70599.72582610537D0 
! a= 0.3097398289189322D0 
! b= 0.0004442895257152834D0 
! c= -1.697501412889046D-07 
! d= 2.880070856747148D-11 
! e= -1.577741778149674D-15 
!-----------------------------------------------------------

  y=0.3097398289189322D0+x*(0.0004442895257152834D0+&
   &x*(-1.697501412889046D-07+x*(2.880070856747148D-11+&
   &x*(-1.577741778149674D-15)))) 
  eqn6001=y
  RETURN
END
  
