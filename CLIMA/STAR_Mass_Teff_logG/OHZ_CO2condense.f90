!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ohzflux

 a= 0.3536743027838826D0 
 b= 7.369074296556717D-05 
 c= -2.75812393722029D-08 
 d= 5.472475628180711D-12 
 e= -3.350995109641247D-16 
 teff = 5780.0
  
 constant = a + b*teff + c*teff**2 + d*teff**3 + e*teff**4
 order1   = b + 2.0*c*teff + 3.0*d*teff**2 + 4.0*e*teff**3
 order2   = c + 3.0*d*teff + 6.0*e*teff**2
 order3   = d + 4.0*e*teff
 order4   = e

 print *,constant,order1,order2,order3,order4

 x = 2600.0

   do while(x <=7201.d0)
     y=eqn6001_CO2Condense(x)
     ohzflux = constant + order1*(x-teff) + order2*(x-teff)**2 + &
               order3*(x-teff)**3 + order4*(x-teff)**4
     print *,x,ohzflux
    x = x + 200.0

   enddo
20 END
  
!-----------------------------------------------------------
REAL*8 FUNCTION eqn6001_CO2Condense(x)
!-----------------------------------------------------------
! TableCurve C:\Users\winey-akka\Desktop\HZ_files\OHZ_CO2condense.f90 May 6, 2013 10:45:35 PM 
! C:\Users\winey-akka\Desktop\OHZ_stars_CO2condense.dat 
! X=  
! Y=  
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9998252933867361D0 
! r2adj=0.9997767637719406D0 
! StdErr=0.0007473345691503415D0 
! Fval=27183.68844122969D0 
! a= 0.3536743027838826D0 
! b= 7.369074296556717D-05 
! c= -2.75812393722029D-08 
! d= 5.472475628180711D-12 
! e= -3.350995109641247D-16 
!-----------------------------------------------------------
  REAL*8 x,y
  y=0.3536743027838826D0+x*(7.369074296556717D-05+&
   &x*(-2.758123937220290D-08+x*(5.472475628180711D-12+&
   &x*(-3.350995109641247D-16)))) 
  eqn6001_CO2Condense=y
  RETURN
END
  
