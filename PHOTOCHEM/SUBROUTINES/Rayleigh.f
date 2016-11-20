      SUBROUTINE RAYLEIGH(wleff,ncomp,icomp,volmix,SIGR2)

cccccccccccccccccccccccccccccc  r a y l e i g h cccccccccccccccccccccccc
cc                                                                    cc
cc    p u r p o s e :                                                 cc
cc                                                                    cc
cc    this subroutine computes the rayleigh scattering cross section  cc
cc    per molecule (meters**-2)  for atmospheres composed of (1) air, cc
cc    (2) co2, (3) n2, (4) o2, (5) h2, (6) he, or any combination of  cc
cc    these gases.                                                    cc
cc    to find the rayleigh scattering optical depth at any level of   cc
cc    the atmosphere, these cross sections must be multiplied by      cc
cc    the pathlength-integrated number density, n(z)*dz.  in a        cc
cc    hydrostatic atmosphere, n(z)*dz = a0*dp/(u0*grav), where u0     cc
cc    is the mean molecular mass (kg/mole), and a0 is avagadro's      cc
cc    number(kg/kmole)                                                cc
cc                                                                    cc
cc    r e f e r e n c e s :                                           cc
cc                                                                    cc
cc    e.j. mccartney, optics of the atmosphere, wiley, p. 187-215,    cc
cc          1976.                                                     cc
cc    a.t. young, revised depolarization corrections for atmospheric  cc
cc          extinction, appl. opt. 19, 3427-3428, 1980.               cc
cc    c.w. allen, astrophysical quantities, athlone press, p. 87,     cc
cc         1964.                                                      cc
cc                                                                    cc
cc    i n p u t :                                                     cc
cc                                                                    cc
cc    ncomp = number of major atmosphere constituents                 cc
cc    icomp = atmosphere constituent index                            cc
cc            (1) air  (2) co2  (3) n2  (4) o2 (5) h2 (6) he          cc
cc   volmix = volume mixing ratio of each constituent                 cc
cc    wleff = effective wavelength (microns)                          cc
cc                                                                    cc
cc    o u t p u t :                                                   cc
cc                                                                    cc
cc    sigr2  - raleigh scattering cross section per molecule (cm**2)  cc
cc                                                                    cc
cccccccccccccccccccccccccccccc  r a y l e i g h  ccccccccccccccccccccccc

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 delta(6),a(6),b(6),wl2i,r,aniso,r2,sum,cnst
      dimension volmix(10,NZ),ncomp(NZ),icomp(10,NZ),SIGR2(NZ)

c       depolarization factors for air, co2, n2, o2 (young, 1980)
      data delta/0.02790e0,0.0780e0,0.0210e0,0.0580e0,0.0000,0.0000/

c***    wavelength dependence coefficients for the refractive index
c       (allen, 1964) (note: wavelengths must be in microns)

      data a/2.871e-04,4.39e-04,2.906e-04,2.663e-04,1.358e-4,3.48e-5/
      data b/5.67e-03,6.4e-03,7.7e-03,5.07e-03,7.52e-3,2.3e-3/

c****   Define the constant multiplier, cnst.  this 
c	constant is equal to 24*pi**3/(1.e-24*L**2), 
c	where L = loschmidt's number (mks units) ,
c	(L = 2.687e25 molecules/m3) and the factor 1.e-24 is 
c	needed to convert wl**4 from microns to meters.

      data cnst /1.031e-24/

      wl2i=1.0e0/(wleff*wleff)

      do j=1,nz
      sum = 0.
c      print *, ncomp(j)
        do  i=1,ncomp(j)
c          print *, icomp(i,j)
            r = 1.0e0 + a(icomp(i,j))*(1.0e0 + b(icomp(i,j))*wl2i)
            aniso = (6.0e0 + 3.0e0*delta(icomp(i,j)))/
     $              (6.0e0 - 7.0e0*delta(icomp(i,j)))
            r2 = r*r + 2.0e0
           sum = sum + volmix(i,j)*aniso*
     $            ((r*r - 1.)/r2)**2
            
        enddo
       !1d4 converts from m^2 to cm^2 for photohcemical model
       sigr2(j)=cnst*wl2i*wl2i*sum*1e4
      enddo

      return
      end
