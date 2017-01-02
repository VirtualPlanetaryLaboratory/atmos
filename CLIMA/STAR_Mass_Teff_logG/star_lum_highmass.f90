!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
   real *8 x,y,teff
  
  open(9,file='teff_lum_highmass.dat')
  teff=6400.0d0
  do while(teff <=7201.0)
    x = teff
    y=eqn1270(x)
    write(9,*)x,10.**y
    PRINT *,x,10.**y
    teff = teff + 200.0d0
  enddo
  close(9)
20 END
  
!-----------------------------------------------------------
REAL*8 FUNCTION eqn1270(x)
!-----------------------------------------------------------
! TableCurve C:\Users\winey-akka\Desktop\star.f90 Oct 22, 2012 11:02:25 AM 
! C:\Users\winey-akka\Desktop\combined_lum_highmass.dat 
! X= 6.45208E+03 
! Y= 3.48017E-01 
! Eqn# 1270  lny=a+bx^2+cx^3 
! r2=0.9999999996057281D0 
! r2adj=0.9999999988171842D0 
! StdErr=6.750846012248326D-06 
! Fval=1268160231.700293D0 
! a= -13.74665758038148D0 
! b= 6.91651378688345D-07 
! c= -5.994820106752163D-11 
!-----------------------------------------------------------
  REAL*8 x,y
  REAL*8 x1,x2 
  x1=x*x 
  x2=x*x*x 
  y=-13.74665758038148D0+6.916513786883450D-07*x1&
   &-5.994820106752163D-11*x2 
  y=DEXP(y) 
  eqn1270=y
  RETURN
END
  
