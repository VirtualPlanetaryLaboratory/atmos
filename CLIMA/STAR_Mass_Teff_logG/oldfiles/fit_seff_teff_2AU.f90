!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
 REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ohzflux

 a= 0.03415600290497983D0 
 b= 0.0001068727657748935D0 
 c= -3.789727936944653D-08 
 d= 6.679972705268996D-12 
 e= -3.828231890393183D-16 
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
! TableCurve C:\Users\winey-akka\Desktop\fit_seff_teff_2AU.f90 Nov 12, 2012 9:19:47 PM 
! C:\Users\winey-akka\Desktop\t3.    
! X=  
! Y=  
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9998501068251153D0 
! r2adj=0.9998084698320918D0 
! StdErr=0.0006462175551203522D0 
! Fval=31684.4847076821D0 
! a= 0.03415600290497983D0 
! b= 0.0001068727657748935D0 
! c= -3.789727936944653D-08 
! d= 6.679972705268996D-12 
! e= -3.828231890393183D-16 
!-----------------------------------------------------------
  REAL*8 x,y
  y=0.03415600290497983D0+x*(0.0001068727657748935D0+&
   &x*(-3.789727936944653D-08+x*(6.679972705268996D-12+&
   &x*(-3.828231890393183D-16)))) 
  eqn6001=y
  RETURN
END
  
