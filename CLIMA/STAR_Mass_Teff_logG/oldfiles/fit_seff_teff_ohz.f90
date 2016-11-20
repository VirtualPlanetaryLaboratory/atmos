!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ohzflux

 a= 0.04726189202188948D0 
 b= 0.000147919392148732D0 
 c= -5.245184318275036D-08 
 d= 9.245313368028642D-12 
 e= -5.298370673711571D-16 
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
     ohzflux = constant + order1*(x-teff) + order2*(x-teff)**2 + &
               order3*(x-teff)**3 + order4*(x-teff)**4
     print *,x,ohzflux
    x = x + 200.0
  enddo
  


20 END
  
!-----------------------------------------------------------
REAL*8 FUNCTION eqn6001(x)
!-----------------------------------------------------------
! TableCurve C:\Users\winey-akka\Desktop\fit_seff_teff_ohz.f90 Oct 23, 2012 12:02:48 PM 
! C:\Users\winey-akka\Desktop\HZ_periods.dat 
! X= 2.60000E+03 
! Y= 2.15019E-01 
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9998500665654652D0 
! r2adj=0.9998084183892056D0 
! StdErr=0.0008944847424979146D0 
! Fval=31675.97561492581D0 
! a= 0.04726189202188948D0 
! b= 0.000147919392148732D0 
! c= -5.245184318275036D-08 
! d= 9.245313368028642D-12 
! e= -5.298370673711571D-16 
!-----------------------------------------------------------
  REAL*8 x,y
  y=0.04726189202188948D0+x*(0.0001479193921487320D0+&
   &x*(-5.245184318275036D-08+x*(9.245313368028642D-12+&
   &x*(-5.298370673711571D-16)))) 
  eqn6001=y
  RETURN
END
  
