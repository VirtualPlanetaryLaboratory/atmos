!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
   real *8 order1,order2,order3,order4,constant
   real *8 ihzflux

 a= 1.201940515647336D0 
 b= 0.0002618173936222995D0 
 c= -1.037406152665933D-07 
 d= 1.999352646183429D-11 
 e= -1.195067188796266D-15
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
! TableCurve C:\Users\winey-akka\Desktop\HZ_files\IHZ_recentVenus.f90 Feb 25, 2013 8:22:25 PM 
! C:\Users\winey-akka\Desktop\HZ_files\recentvenusearlymars.dat 
! X= 2600.00000000000 
! Y= 1.47779869879829 
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9998756081295525D0 
! r2adj=0.999841054832206D0 
! StdErr=0.00194791030836003D0 
! Fval=38181.02518702024D0 
! a= 1.201940515647336D0 
! b= 0.0002618173936222995D0 
! c= -1.037406152665933D-07 
! d= 1.999352646183429D-11 
! e= -1.195067188796266D-15 
!-----------------------------------------------------------
  REAL*8 x,y
  y=1.201940515647336D0+x*(0.0002618173936222995D0+&
   &x*(-1.037406152665933D-07+x*(1.999352646183429D-11+&
   &x*(-1.195067188796266D-15)))) 
  eqn6001=y
  RETURN
END
  
