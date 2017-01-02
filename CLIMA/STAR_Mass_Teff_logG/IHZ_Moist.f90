!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ihzflux

 a= 0.6865332809420407D0 
 b= 0.0001495467969597905D0 
 c= -5.925533255493829D-08 
 d= 1.142005044405594D-11 
 e= -6.826073232324185D-16 
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
! TableCurve C:\Users\winey-akka\Desktop\HZ_files\IHZ_Moist.f90 Feb 25, 2013 8:16:46 PM 
! C:\Users\winey-akka\Desktop\HZ_files\IHZ_stars.dat 
! X= 2600.0 
! Y= 0.8441 
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9998756081295524D0 
! r2adj=0.9998410548322059D0 
! StdErr=0.001112621829092534D0 
! Fval=38181.02518698254D0 
! a= 0.6865332809420407D0 
! b= 0.0001495467969597905D0 
! c= -5.925533255493829D-08 
! d= 1.142005044405594D-11 
! e= -6.826073232324185D-16 
!-----------------------------------------------------------
  REAL*8 x,y
  y=0.6865332809420407D0+x*(0.0001495467969597905D0+&
   &x*(-5.925533255493829D-08+x*(1.142005044405594D-11+&
   &x*(-6.826073232324185D-16)))) 
  eqn6001=y
  RETURN
END
  
