!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ihzflux

 a= 0.2532950578807256D0 
 b= 0.0006074633221383245D0 
 c= -2.289434566541145D-07 
 d= 3.753451577274881D-11 
 e= -2.015477728793372D-15 
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
! TableCurve C:\Users\winey-akka\Desktop\fit_seff_teff_runaway.f90 Oct 23, 2012 2:32:46 PM 
! C:\Users\winey-akka\Desktop\runaway.dat 
! X=  
! Y=  
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9999277650566247D0 
! r2adj=0.999907699794576D0 
! StdErr=0.001529930806908031D0 
! Fval=65752.89828001567D0 
! a= 0.2532950578807256D0 
! b= 0.0006074633221383245D0 
! c= -2.289434566541145D-07 
! d= 3.753451577274881D-11 
! e= -2.015477728793372D-15 
!-----------------------------------------------------------
  REAL*8 x,y
  y=0.2532950578807256D0+x*(0.0006074633221383245D0+&
   &x*(-2.289434566541145D-07+x*(3.753451577274881D-11+&
   &x*(-2.015477728793372D-15)))) 
  eqn6001=y
  RETURN
END
  
