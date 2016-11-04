!-----------------------------------------------------------
!         TableCurve FORTRAN F90 Code Output
!    To modify generated output, edit FORT90.TCL
!-----------------------------------------------------------

PROGRAM main
  REAL*8 x,y,teff

  open(9,file='teff_mass_highmass.dat')
  teff=6400.0d0
  do while(teff <=7201.0)
    x = teff
    y=eqn1147(x)
    write(9,*)x,y
    PRINT *,x,y
    teff = teff + 200.0d0
  enddo
  close(9)

20 END
  
!-----------------------------------------------------------
REAL*8 FUNCTION eqn1147(x)
!-----------------------------------------------------------
! TableCurve C:\Users\winey-akka\Desktop\star_mass_highmass.f90 Oct 22, 2012 11:14:27 AM 
! C:\Users\winey-akka\Desktop\combined_lum_highmass.dat 
! X= 6.45208E+03 
! Y= 1.21627E+00 
! Eqn# 1147  y=a+b(lnx)^2+clnx 
! r2=0.9999903120009821D0 
! r2adj=0.9999709360029464D0 
! StdErr=0.0006639821440949857D0 
! Fval=51609.74470349887D0 
! a= 141.2598466600419D0 
! b= 2.065146622776594D0 
! c= -34.08034816557594D0 
!-----------------------------------------------------------
  REAL*8 x,y
  REAL*8 x1,x2 
  x1=DLOG(x)*DLOG(x) 
  x2=DLOG(x) 
  y=141.2598466600419D0+2.065146622776594D0*x1&
   &-34.08034816557594D0*x2 
  eqn1147=y
  RETURN
END
  
