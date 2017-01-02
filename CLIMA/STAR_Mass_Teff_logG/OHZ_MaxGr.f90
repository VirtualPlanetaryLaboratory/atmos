!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,a,b,c,d,e,teff
  real *8 order1,order2,order3,order4,constant
  real *8 ohzflux
  
  a= 0.06310316806859836D0 
 b= 0.0001400768156763158D0 
 c= -5.029249936658971D-08 
 d= 8.999267479704257D-12 
 e= -5.192508276748419D-16
 teff=5780.0

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
! TableCurve C:\Users\winey-akka\Desktop\HZ_files\OHZ_MaxGr.f90 Feb 26, 2013 9:32:49 AM 
! # TEFF (K)       SEFF 
! X=  
! Y=  
! Eqn# 6001  y=a+bx+cx^2+dx^3+ex^4 
! r2=0.9998560592583773D0 
! r2adj=0.9998160757190376D0 
! StdErr=0.0008850356429469139D0 
! Fval=32994.94102875677D0 
! a= 0.06310316806859836D0 
! b= 0.0001400768156763158D0 
! c= -5.029249936658971D-08 
! d= 8.999267479704257D-12 
! e= -5.192508276748419D-16 
!-----------------------------------------------------------
  REAL*8 x,y
  y=0.06310316806859836D0+x*(0.0001400768156763158D0+&
   &x*(-5.029249936658971D-08+x*(8.999267479704257D-12+&
   &x*(-5.192508276748419D-16)))) 
  eqn6001=y
  RETURN
END
  
