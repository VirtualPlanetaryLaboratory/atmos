c****************************** plume **********************************
c
      subroutine plume ( ytflx, yqflx, yuflx, yvflx, yt, yq, yu, yv,
     *                   ylwc_c, ylwc_a, ylwc_s, 
     *                   fplumt, fplumq, fplumu, fplumv,
     *                   tracer, ypm, ypthic, ypmcap,
     *                   wplume0, wplume1, plumeprec, heimin, heimax,
     *                   dtplume)

c     Does convection (both pbl and free), by solving a subgrid
c     plume model. For pbl convection, the plume is initiated  
c     at the center of the lowest layer (ie, top of constant-flux layer)
c     with w,t,etc values scaled depending on the surface heat/momentum
c     fluxes based on constant-flux scaling, and is integrated 
c     vertically until it stops.

c     For free plumes, there is an option (commented by "csing" or 
c     "cmult") to do just one vertical integration from the bottom to
c     the top (csing) or to do multiplte integrations stating from
c     each layer in turn (cmult). In both cases, a free plume is 
c     initiated at the center of a model layer if it can rise to 
c     (or beyond) the center of the next layer. Plumes are integrated
c     vertically upwards from layer midpoint to midpoint until the 
c     plume stops. 

c     Condensation can occur within the plumes, forming wthin-plume
c     liquid cloud water, which can precipitate out if it exceeds a
c     maximum value.

c     The plume equations explicitly solve for plume w,theta,q,u,v where
c     all except w are perturbations from the ambient quantities. The
c     plume vertical finite-difference step is upwards from layer
c     midpoint to layer midpoint. These quantities, along with the plume
c     condensation amounts, then imply the large-scale fluxes at layer
c     midpoints and the large-scale precip and latent heat release.
c     These fluxes are accumulated over pbl and free
c     calculations, and then "shifted" downwards to the next lowest
c     layer interface to calculate large-scale convergences.
c     The large-scale flux convergences for heat, water vapor and 
c     momentum are actually put in reservoirs fplum* and released 
c     (see reserv, called from vdif) with a time scale of a few hours.
c     Liquid water cloud amounts and tracers are updated here with
c     no reservoir.

c     Lines beginning with "cloc" are to convert to standalone mode
c     for use with ./Plume/plume.f

c     supplied:
c     ytflx   = upward sensible heat flux from surface (W/m2)
c     yqflx   = upward h2o mass flux from surface (Kg/m2/s)
c     yuflx   = upward u-momentum flux from surface (N/m2)
c     yvflx   = upward v-momentum flux from surface (N/m2)
c     yt      = agcm temperatures (deg K)
c     yq      = agcm specific humidity (Kg/Kg)
c     yu      = agcm eastward  velocity (m/s)
c     yv      = agcm northward velocity  (m/s)
c     ypm     = agcm mid-layer pressures (N/m2)
c     ypthic  = agcm layer pressure-thicknesses (N/m2)
c     ypmcap  = (ypm/surface pressure)**cappa (N/m2)
c     dtplume = plume model timestep (s)

c     modified:
c     ylwc_c  = convective liquid water content (3-D)
c     ylwc_a  = anvil-cirrus liquid water content (3-D) 
c     ylwc_s  = stratiform   liquid water content (3-D)
c     fplum*  = reservoirs of "flux convergence*time" in each agcm layer
c     tracer  = agcm tracer fields. ylwc_[c,a,s] are equivalenced to
c               tracers #1-3, and are changed explicitly below.
c               fplum* are equivalenced to tracers #4-7. Only
c               tracers #8 to #ntrace are mixed as "tracers" below
c               (see local params ntraca,ntracb)
c     wplume0,1 = zonal mean plume vertical velocities (for diagnostics)
c     plumeprec = within-plume precipitation formation rate (3-D)
c                 (for reevap)
c     heimin    = daily min pbl plume heights (for history)
c     heimax    = daily max pbl plume heights (for history)

c ---------------------------------------------------------------------
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      INCLUDE 'complume'
c-----------------------------------------------------------------------

      dimension 
     *  ytflx(nlon,norec),        yqflx(nlon,norec),
     *  yuflx(nlon,norec),        yvflx(nlon,norec),

     *  yt(nlon,norec,nlev),      yq(nlon,norec,nlev),
     *  yu(nlon,norec,nlev),      yv(nlon,norec,nlev), 
     *  ylwc_c(nlon,norec,nlev),  ylwc_a(nlon,norec,nlev),
     *  ylwc_s(nlon,norec,nlev), 
     *  fplumt(nlon,norec,nlev),  fplumq(nlon,norec,nlev),
     *  fplumu(nlon,norec,nlev),  fplumv(nlon,norec,nlev), 
     *  tracer(nlon,norec,nlev,ntrace),

     *  ypm(nlon,norec,nlev),
     *  ypthic(nlon,norec,nlev),  ypmcap(nlon,norec,nlev),

     *  wplume0(norec,nlev),      wplume1(norec,nlev),
     *  plumeprec(nlon,norec,nlev),
     *  heimin(nlon,norec),       heimax(nlon,norec)

      dimension
     *  yh(nlon,nlev),         yrho(nlon,nlev),
     *  ylwc_t(nlon,nlev),

     *  dzu(nlon,nlev),        dzl(nlon,nlev),
     *  dztot(nlon,nlev),      dz(nlon,nlev),
     *  dhdz(nlon,nlev),       dqdz(nlon,nlev),
     *  dudz(nlon,nlev),       dvdz(nlon,nlev),
     *  dlwcdz(nlon,nlev),     dadz(nlon,nlev,ntrace),
     *  zelev(nlon,nlev)

      dimension
     *  zh(nlon,nlev),        zq(nlon,nlev), 
     *  zu(nlon,nlev),        zv(nlon,nlev),    
     *  zw(nlon,nlev),        zws(nlon,nlev),
     *  zt(nlon,nlev),        zrho(nlon,nlev),    
     *  zlw(nlon,nlev),       za(nlon,nlev,ntrace),
     *  zcond(nlon,nlev),     zlath(nlon,nlev),
     *  zprec(nlon,nlev),

     *  zconda(nlon,0:nlev),  zheata(nlon,0:nlev),
     *  zpreca(nlon,0:nlev),  zmfa(nlon,0:nlev),
     *  tplu(nlon,0:nlev),    qplu(nlon,0:nlev), 
     *  uplu(nlon,0:nlev),    vplu(nlon,0:nlev),
     *  cpluc(nlon,0:nlev),   cplua(nlon,0:nlev),
     *  cplus(nlon,0:nlev),   aplu(nlon,0:nlev,ntrace)

c     tracers 1-3 are used for liquid water contents (mixed by plume),
c     and tracers 4-7 are used for flux reservoirs (not mixed by plume).
      equivalence 
     *  (dlwcdz,  dadz(1,1,1)), 
     *  (zlw,     za(1,1,1)),
     *  (cpluc,   aplu(1,0,1)),
     *  (cplua,   aplu(1,0,2)),
     *  (cplus,   aplu(1,0,3))
      parameter (ntraca=ntrace, ntracb=ntrace) !EWS - ntraca was hardcoded to '8', but could never go up that high due
                                               ! to dimensions of declared variabiables, so I reset it to 'ntrace'. 
                                               ! Hopefully didn't screw something up, need to test - 9/4/2015

      dimension
     *  zarea(nlon,nlev),      zrad(nlon,nlev),
     *  zent(nlon),            zdrag(nlon),

     *  zcfl(nlon),            zmf(nlon,nlev),
     *  ztotc1(nlon),          ztotc2(nlon),
     *  ztota1(nlon),          ztota2(nlon),
     *  ztots1(nlon),          ztots2(nlon),
     *  zlwcold(nlon,nlev)

      dimension
     *  heipbl(nlon),          
     *  laypbl(nlon),           numpbl(nlev),
     *  laytop(nlon),           numfret(nlev),
     *  laybot(nlon),           numfreb(nlev),
     *  wfre1(nlon,nlev),       wfre2(nlon,nlev)

      logical ifpbl,            ifanv(nlon)
      character*3 chaplu(nlev)

c     Diagnostic only:
      dimension qrat(nlev), nqrat(nlev), qneg(nlev), nqneg(nlev)
c     parameter (nprinz=1)
      parameter (nprinz=2)                                         !cloc
      dimension iuprinz(nprinz)

      parameter (pi = 3.14159265358979)

c ---------------------------------------------------------------------
c     statement function for plume radius vs. elevation (not used):
c     plumerad (z) = plumerad1 + min (0.5*z, plumerad2)
c-----------------------------------------------------------------------
c#include "comsatplume"
      INCLUDE 'comsatplume'
c-----------------------------------------------------------------------

      cpair  = 1004.64
      gravit = 9.80616
      latvap = 2.5104e06
      latice = 3.336e5
      latsub = latvap + latice
      rair   = 287.04
      rh2o   = 4.61e2
      rhoh2o = 1.e3
      ch2o = 4.218e3
      cice = 2.106e3                  
      cappa  = rair/cpair
      zvir   = rh2o/rair -1.
      cpwv   = 1.81e3
      cpvir  = cpwv/cpair - 1.

      nstep = 0
      nrstrt = nstep
      twodt = 2.*dtime
      dtbud = dtime
      facbud = 1.

      alw   = .95
      blw_c = 1.
      blw_a = 1.
      blw_s = 1.
      clw_c = .00010
      clw_a = .00003
      clw_s = .00003
      clw_p = .0010
      dlw_l = 20.
      dlw_i = .5
      elw   = .3
      flw_c = 0.2
      flw_a = 0.5
      flw_s = 1.0
      glw = 0.3 ! (Smith, 1990, eq.2.29, c_a)
      tliqa = -15.
      tliqb =  10.
      siganvil = .45








c        Initialize local diagnostics

      zptot = 0.
      zetot = 0.
      zqtot = 0.
      do 10 jk=1,nlev
        numpbl(jk) = 0
        numfret(jk) = 0
        numfreb(jk) = 0
   10 continue

c        Overall loop over latitude

c***********************
      do 1000 jj=1,norec
c***********************

c         Calculate ambient pot.temp, density and total cloud liquid
c         water at layer midpoints

      do 100 jk=1,nlev
        do 102 ji=1,nlon
          yh(ji,jk)     = yt(ji,jj,jk) / ypmcap(ji,jj,jk)
          yrho(ji,jk)   = ypm(ji,jj,jk)
     *                  / (rair*(1.+zvir*yq(ji,jj,jk))*yt(ji,jj,jk))
          ylwc_t(ji,jk) = ylwc_c(ji,jj,jk) + ylwc_a(ji,jj,jk) 
     *                  + ylwc_s(ji,jj,jk)
  102   continue
  100 continue

c        Compute ambient dz, d*/dz for full layers. For bottom layer,
c        just use top half of layer.

      do 120 jk=1,nlev
        do 122 ji=1,nlon
          dztot(ji,jk) = ypthic(ji,jj,jk) / (gravit*yrho(ji,jk))
          dzu(ji,jk) = dztot(ji,jk) * (sig(jk)-sigkmh(jk))/dsigma(jk)
          dzl(ji,jk) = dztot(ji,jk) - dzu(ji,jk)
  122   continue
  120 continue

      do 124 jk=1,nlev-1
        do 126 ji=1,nlon
          dz(ji,jk) = dzl(ji,jk) + dzu(ji,jk+1)
          dhdz(ji,jk)   = (yh(ji,jk)     - yh(ji,jk+1)    ) / dz(ji,jk)
          dqdz(ji,jk)   = (yq(ji,jj,jk)  - yq(ji,jj,jk+1) ) / dz(ji,jk)
          dudz(ji,jk)   = (yu(ji,jj,jk)  - yu(ji,jj,jk+1) ) / dz(ji,jk)
          dvdz(ji,jk)   = (yv(ji,jj,jk)  - yv(ji,jj,jk+1) ) / dz(ji,jk)
          dlwcdz(ji,jk) = (ylwc_t(ji,jk) - ylwc_t(ji,jk+1)) / dz(ji,jk)
  126   continue

        if (ntraca.le.ntracb) then
          do 1260 n=ntraca,ntracb
            do 1262 ji=1,nlon
             dadz(ji,jk,n) = (tracer(ji,jj,jk,n)-tracer(ji,jj,jk+1,n))
     *                       /dz(ji,jk)
 1262       continue
 1260     continue
        endif
  124 continue

c        Compute elevations of layer midpoints (for plume radii calc)

      do 130 ji=1,nlon
        zelev(ji,nlev) = dzl(ji,nlev)
  130 continue
      do 132 jk=nlev-1,1,-1
        do 134 ji=1,nlon
          zelev(ji,jk) = zelev(ji,jk+1) + dz(ji,jk)
  134   continue
  132 continue

c        Zero plume fluxes accumulated over loop 500

      call zero (tplu,  nlon*nlevp)
      call zero (qplu,  nlon*nlevp)
      call zero (uplu,  nlon*nlevp)
      call zero (vplu,  nlon*nlevp)
      call zero (cpluc, nlon*nlevp)
      call zero (cplua, nlon*nlevp)
      call zero (cplus, nlon*nlevp)
      if (ntraca.le.ntracb)
     *  call zero (aplu(1,0,ntraca), nlon*nlevp*(ntracb-ntraca+1))
      call zero (zconda, nlon*nlevp)
      call zero (zheata, nlon*nlevp)
      call zero (zpreca, nlon*nlevp)
      call zero (zmfa,   nlon*nlevp)


c        Loop twice, for "pbl" plumes (ifpbl=.true.) and for "free"
c        plumes (ifpbl=.false., start at second lowest layer)

c*********************
c     do 500 itype=0,1
      do 500 itype=0,0   ! for pbl  only (diagnostics ok only for one)
c     do 500 itype=1,1   ! for free only (or the other)            !cloc
c*********************

      if (itype.eq.0) then
        ifpbl = .true.
        jkinita = nlev
        jkinitb = nlev
      else
        ifpbl = .false.
        jkinita = nlev  
        jkinitb = nlev                                            !csing
cmult   jkinitb = 2                                               !cmult
      endif

c         Initialize diagnostics for setting of plume velocities
c         (wfre*,wplume*) and indices of bottom,top layers (lay*, num*)

      call zero (wfre1, nlon*nlev)
      call zero (wfre2, nlon*nlev)

      do 140 ji=1,nlon
        heipbl(ji) = dzl(ji,nlev)
        laypbl(ji) = nlevp
        laytop(ji) = nlevp
        laybot(ji) = 0
  140 continue

c************************************
      do 510 jkinit = jkinita,jkinitb
c************************************

c        Assume no new plumes initiated above 150 mb

      if (sig(jkinit).lt.0.150) goto 510

c        Initialize all plume variables to zero or ambient, and
c        initialize diagnostic plume top/bottom locations 

      call zero (zw,   nlon*nlev)
      call zero (zws,  nlon*nlev)
      call zero (zt,   nlon*nlev)
      call zero (zh,   nlon*nlev)
      call zero (zq,   nlon*nlev)
      call zero (zu,   nlon*nlev)
      call zero (zv,   nlon*nlev)
      call zero (zlw,  nlon*nlev)
      call scopy (nlon*nlev, yrho, 1, zrho, 1)
      if (ntraca.le.ntracb)
     *  call zero (za(1,1,ntraca), nlon*nlev*(ntracb-ntraca+1))
      call zero (zarea,nlon*nlev)
      call zero (zcond, nlon*nlev)
      call zero (zlath, nlon*nlev)
      call zero (zprec, nlon*nlev)

c       For "pbl" plumes, initialize plume conditions
c       at midpoint of bottom layer based on surface fluxes from 
c       constant-flux layer. If net buoyancy flux 
c       (ytflx/cp+.622*yqflx) is < 0, leave initial plume z* = zero.

c--------------------
      if (ifpbl) then
c--------------------

        do 150 ji=1,nlon

c         zarea(ji,nlev) = plumearea1  ! .005
c         zrad(ji,nlev) = plumerad1    ! 100.
          zarea(ji,nlev) = xare1                                   !cloc
          zrad (ji,nlev) = xrad1                                   !cloc

          if ( ytflx(ji,jj)/cpair+.622*yqflx(ji,jj) .gt. 0. ) then

            ztau = sqrt (yuflx(ji,jj)**2 + yvflx(ji,jj)**2)
            ztauv = max (.01, sqrt(ztau/yrho(ji,nlev)) )
            zconv = ( (max(ytflx(ji,jj),0.)/(yrho(ji,nlev)*cpair))
     *                * rair * (1.-sig(nlev))
     *              ) ** (1./3.)
c-----
c old:
c           zwt = ytflx(ji,jj) / (cpair*yrho(ji,nlev)*zarea(ji,nlev))
c           zw(ji,nlev) = max (ztauv, (abs(zwt))**(1./3.))
c           zt(ji,nlev) = zwt / zw(ji,nlev)
c           if (abs(zt(ji,nlev)).gt.10.) then
c             zt(ji,nlev) = min (10., max(-10., zt(ji,nlev)))
c             zw(ji,nlev) = zwt / zt(ji,nlev)
c           endif
c           zmfbot = zarea(ji,nlev) * yrho(ji,nlev) * zw(ji,nlev)
c           zq(ji,nlev) = yqflx(ji,jj) / zmfbot
c           zu(ji,nlev) = yuflx(ji,jj) / zmfbot
c           zv(ji,nlev) = yvflx(ji,jj) / zmfbot
c-----
c new:
            zfac = 1. / zarea(ji,nlev)
            zwt = ytflx(ji,jj) / (yrho(ji,nlev)*cpair)
            zw(ji,nlev) = sqrt(zfac) * max(ztauv,zconv) 
            zt(ji,nlev) = zfac*zwt / zw(ji,nlev)
            if (abs(zt(ji,nlev)).gt.10.) then
              zt(ji,nlev) = min (10., max(-10., zt(ji,nlev)))
              zw(ji,nlev) = zfac*zwt / zt(ji,nlev)
            endif
            zq(ji,nlev)= zfac*yqflx(ji,jj)/(yrho(ji,nlev)*zw(ji,nlev))
            zu(ji,nlev)= zfac*yuflx(ji,jj)/(yrho(ji,nlev)*zw(ji,nlev))
            zv(ji,nlev)= zfac*yvflx(ji,jj)/(yrho(ji,nlev)*zw(ji,nlev))

            zh(ji,nlev) = zt(ji,nlev) / ypmcap(ji,jj,nlev)
            zws(ji,nlev)= zw(ji,nlev)**2
            zrho(ji,nlev)= ypm(ji,jj,nlev)
     *                  / ( rair*(1.+zvir*(zq(ji,nlev)+yq(ji,jj,nlev)))
     *                      *(zt(ji,nlev)+yt(ji,jj,nlev)) )

          endif

  150   continue

c---------
      else
c---------

c          For "free" plumes, set initial plume area and radius. 
c          Nb: plume radius affects the physics only via the 
c          entrainment and drag coefficients zent,zdrag. 
c          Plume area does not affect within-plume physics, only the
c          net flux convergences for changing large-scale quantities.
       
        jk = jkinit
        do 160 ji=1,nlon
          zarea(ji,jk) = xare2                                     !cloc
          zrad (ji,jk) = xrad2                                     !cloc
  160   continue

c       Free plume radius dep on column precipitable water vapor:
c       call zero (zpwv, nlon)
c       do 160 jk=1,nlev
c         do 162 ji=1,nlon
c           zpwv(ji) = zpwv(ji) + yq(ji,jj,jk)*ypthic(ji,jj,jk)
c 162     continue
c 160   continue
c       do 164 ji=1,nlon
c         zarea(ji,jkinit) = plumearea2
c         zrad(ji,jkinit) = max (plumerad1, plumerad2*(zpwv(ji)/gravit))
c 164   continue

c       do 160 ji=1,nlon
c         zarea(ji,jkinit) = plumearea2
c         zrad(ji,jkinit) = plumerad2
c 160   continue

c----------
      endif
c----------

c        Integrate plumes upwards from initiating layer to top

c================================
      do 200 jk = jkinit-1, 1, -1
c================================

        jkp = jk+1

c          Step plume from midpoint of layer jkp to midpoint of layer jk

        do 210 ji=1,nlon

c            Force for free plume, or non-zero pbl plume from below

          if( (.not.ifpbl)                   .or.zw(ji,jkp).gt.0.)!csing
cmult     if(((.not.ifpbl).and.jkp.eq.jkinit).or.zw(ji,jkp).gt.0.)!cmult
     *      then

            zent(ji)  = 0.20 / zrad(ji,jkp)
            zdrag(ji) = 0.5  / (pi*zrad(ji,jkp))

            zexp = exp(-zent(ji)*dz(ji,jk))
            zh(ji,jk) = zh(ji,jkp)*zexp
     *                  - (dhdz(ji,jk)/zent(ji))  *(1.-zexp)
            zq(ji,jk) = zq(ji,jkp)*zexp
     *                  - (dqdz(ji,jk)/zent(ji))  *(1.-zexp)
            zu(ji,jk) = zu(ji,jkp)*zexp
     *                  - (dudz(ji,jk)/zent(ji))  *(1.-zexp)
            zv(ji,jk) = zv(ji,jkp)*zexp
     *                  - (dvdz(ji,jk)/zent(ji))  *(1.-zexp)
            zlw(ji,jk)= zlw(ji,jkp)*zexp
     *                  - (dlwcdz(ji,jk)/zent(ji))*(1.-zexp)
            zt(ji,jk) = zh(ji,jk)*ypmcap(ji,jj,jk)

c              Test for parcel saturation, and condense amount that 
c              results in exactly saturated air (zcond, zlath).
c              Add to plume liquid water content, converting excess 
c              over maximum (namelist param clw_p) to within-plume
c              precip (zprec).
c              zlath is latent heat of condensation, determined by
c              ambient level temperature (ie, rain or snow), plus
c              correction for specific heats to be consistent with lsx
c              (imagining all phase changes occur at tmelt).

            ztabs = zt(ji,jk) + yt(ji,jj,jk)
            zqabs = zq(ji,jk) + yq(ji,jj,jk) 
            zqsat = qstblf (ztabs, ypm(ji,jj,jk))
            if (zqabs.gt.zqsat) then
c             zlath(ji,jk) = cvmgt ( latvap + cpwv*(ztabs-tmelt)
c    *                               - ch2o*(yt(ji,jj,jk)-tmelt),
c    *                               latsub + cpwv*(ztabs-tmelt)
c    *                               - cice*(yt(ji,jj,jk)-tmelt),
c    *                               yt(ji,jj,jk).ge.tmelt )
c             For now, all cloud water is liquid with no specific heat
              zlath(ji,jk) = latvap + cpwv*(ztabs-tmelt)
              dqz = (zlath(ji,jk)/rh2o) * zqsat / (ztabs**2)
              cpz = cpair*(1.+cpvir*zqabs)
              delq = (zqabs-zqsat) / (1.+dqz*zlath(ji,jk)/cpz)
              zqnew = zqabs - delq
              delt = delq*zlath(ji,jk) / (cpair*(1.+cpvir*zqnew))
              ztnew = ztabs + delt
              zq(ji,jk) = zqnew - yq(ji,jj,jk)
              zt(ji,jk) = ztnew - yt(ji,jj,jk)
              zh(ji,jk) = zt(ji,jk)/ypmcap(ji,jj,jk)

              zlwnew = zlw(ji,jk) + ylwc_t(ji,jk) + delq
c             zclw_p = clw_p                               ! (kg/kg dep)
              zclw_p = clw_p / yrho(ji,jk)                 ! (kg/m3 dep)
              delw = max (0., zlwnew-zclw_p)
              zlw(ji,jk) = min (zlwnew, zclw_p) - ylwc_t(ji,jk)
            else
              delq = 0.
              delt = 0.
              delw = 0.
            endif

            zrho(ji,jk) = ypm(ji,jj,jk)
     *                  / ( rair*(1.+zvir*(zq(ji,jk)+yq(ji,jj,jk)))
     *                      *(zt(ji,jk)+yt(ji,jj,jk)) )

            drhog0= gravit * (yrho(ji,jkp)-zrho(ji,jkp))/zrho(ji,jkp)
            drhog1= gravit * (yrho(ji,jk) -zrho(ji,jk) )/zrho(ji,jk) 
            a_rho =  (drhog0 + (drhog1-drhog0)/(1.-zexp)) / zent(ji)
            b_rho = -(         (drhog1-drhog0)/(1.-zexp)) / zent(ji)

            zexp2 = zexp*zexp
c           or: include vertical drag
c           zexp2 = exp(-2.*(zent(ji)+zdrag(ji))*dz(ji,jk))
c           a_rho = a_rho* zent(ji) / (zent(ji)+zdrag(ji))
c           b_rho = b_rho* zent(ji) / (2.*(zent(ji)+zdrag(ji))-zent(ji))

            zws(ji,jk) = zws(ji,jkp)*zexp2
     *                 + a_rho*(1.-zexp2) + 2.*b_rho*(zexp-zexp2)
            zw(ji,jk) = sqrt ( max(zws(ji,jk),0.) )

c             Test if detrainment >=0 is violated. If so, do corrections

            if (zw(ji,jk).gt.0.) then
              drho0 = yrho(ji,jkp)-zrho(ji,jkp)
              drho1 = yrho(ji,jk) -zrho(ji,jk)
              zwsmax = (   zrho(ji,jkp)*zws(ji,jkp)
     *                   + gravit*0.5*(drho0+drho1)*dz(ji,jk)
     *                 ) / zrho(ji,jk)
              if (zwsmax.gt.0. .and. zws(ji,jk).gt.zwsmax) then
                zws(ji,jk) = zwsmax
                zw(ji,jk) = sqrt(zws(ji,jk))
                zma = zrho(ji,jkp)*zw(ji,jkp)
                zmb = zrho(ji,jk) *zw(ji,jk)
                zent(ji) = (zmb-zma)/(0.5*(zma+zmb)*dz(ji,jk))
                zexp = exp(-zent(ji)*dz(ji,jk))
                zh(ji,jk) = zh(ji,jkp)*zexp 
     *                      - (dhdz(ji,jk)/zent(ji))  *(1.-zexp)
     *                      + delt/ypmcap(ji,jj,jk)
                zq(ji,jk) = zq(ji,jkp)*zexp
     *                      - (dqdz(ji,jk)/zent(ji))  *(1.-zexp)
     *                      - delq
                zu(ji,jk) = zu(ji,jkp)*zexp
     *                      - (dudz(ji,jk)/zent(ji))  *(1.-zexp)
                zv(ji,jk) = zv(ji,jkp)*zexp
     *                      - (dvdz(ji,jk)/zent(ji))  *(1.-zexp)
                zlw(ji,jk)= zlw(ji,jkp)*zexp
     *                      - (dlwcdz(ji,jk)/zent(ji))*(1.-zexp)
     *                      + delq - delw
                zt(ji,jk) = zh(ji,jk)*ypmcap(ji,jj,jk)
              endif
            endif

c              If vertical velocity = 0, end of plume;
c              else, re-compute zu and zv with (non-linear) form drag

            if (zw(ji,jk).eq.0.) then

              zh(ji,jk) = 0.
              zq(ji,jk) = 0.
              zu(ji,jk) = 0.
              zv(ji,jk) = 0.
              zlw(ji,jk) = 0.
              zt(ji,jk) = 0.
              zrho(ji,jk) = yrho(ji,jk)
              zcond(ji,jk) = 0.
              zlath(ji,jk) = 0.
              zprec(ji,jk) = 0.

c             Either:
c             heipbl(ji)= heipbl(ji) + dz(ji,jk) * zws(ji,jkp)
c    *                / (max(zws(ji,jkp),1.e-10) - min(zws(ji,jk),0.))
c             Or:
c             Solve quadratic eqn for exp(-zent*z),ie,for z where w=0,
c             for accurate increment to pbl height (diagnostic only).
c             Only valid for no vertical drag(if commented out above).
c             cvmgt(zz1) mod is fudge to avoid quadratic singularity.
              if (ifpbl) then
                zz1 = zws(ji,jkp) - a_rho - 2.*b_rho
                zz1 = cvmgt (zz1, zz1+1.e-9, abs(zz1).gt.1.e-10)
                zz2 = (-b_rho + sqrt(max(0.,b_rho**2-zz1*a_rho))) / zz1
                heipbl(ji) = heipbl(ji) - alog(zz2)/zent(ji)
                laypbl(ji) = jkp
              else
                if (zw(ji,jkp).gt.0.) laytop(ji) = min (laytop(ji),jkp)
              endif

              zws(ji,jk) = 0.

            else

              zenu= zent(ji) + zdrag(ji)*abs(0.5*(zu(ji,jkp)+zu(ji,jk)))
     *                         / (0.5*(zw(ji,jkp)+zw(ji,jk)))
              zenv= zent(ji) + zdrag(ji)*abs(0.5*(zv(ji,jkp)+zv(ji,jk)))
     *                         / (0.5*(zw(ji,jkp)+zw(ji,jk)))
              zexu  = exp(-zenu*dz(ji,jk))
              zexv  = exp(-zenv*dz(ji,jk))
              zu(ji,jk) = zu(ji,jkp)*zexu - (dudz(ji,jk)/zenu)*(1.-zexu)
              zv(ji,jk) = zv(ji,jkp)*zexv - (dvdz(ji,jk)/zenv)*(1.-zexv)

              zcond(ji,jk) = zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delq
              zprec(ji,jk) = zarea(ji,jkp)*zrho(ji,jk)*zw(ji,jk)*delw

              if (ifpbl) then
                heipbl(ji) = heipbl(ji) + dz(ji,jk)
              else
                laybot(ji) = max (laybot(ji),jkp)
              endif

            endif

          endif

  210   continue

c          In same way, integrate plume tracers through curr layer

        if (ntraca.le.ntracb) then
          do 2100 n=ntraca,ntracb
            do 2102 ji=1,nlon
              if (zw(ji,jk).gt.0.) then
                zexp = exp(-zent(ji)*dz(ji,jk))
                za(ji,jk,n) = za(ji,jkp,n)*zexp
     *                      - (dadz(ji,jk,n)/zent(ji))*(1.-zexp)
              else
                za(ji,jk,n) = 0.
              endif
 2102       continue
 2100     continue
        endif

c          Set new plume radius and area.
c          Also count number of plumes, skip out of vertical 
c          integration for pbl or multiple-initiation plumes
c          if none left for this strip.

        nplu = 0
        do 230 ji=1,nlon
          zrad(ji,jk)  = zrad(ji,jkp)
          zarea(ji,jk) = zarea(ji,jkp)
c         zw(ji,jk)  = zw(ji,jk) * zarea(ji,jkp)/zarea(ji,jk)
c         zws(ji,jk) = zw(ji,jk)**2
          if (zw(ji,jk).gt.0.) nplu = nplu + 1
  230   continue
        if (ifpbl .and.  nplu.eq.0) goto 202                      !csing
cmult   if (             nplu.eq.0) goto 202                      !cmult

c=============
  200 continue
  202 continue
c=============
  
c        Accumulate large-scale fluxes at layer midpoints. 
c        Also accumulate mass flux zmfa for cfl check below

      do 300 jk=1,nlev

        do 302 ji=1,nlon
          zmf(ji,jk)  = zarea(ji,jk)*zrho(ji,jk)*zw(ji,jk)
          zmfa(ji,jk) = zmfa(ji,jk) + zmf(ji,jk)
  302   continue

        if (ifpbl .and. jk.eq.nlev) then
          do 304 ji=1,nlon
            tplu(ji,jk) = tplu(ji,jk) + ytflx(ji,jj)
            qplu(ji,jk) = qplu(ji,jk) + yqflx(ji,jj)
            uplu(ji,jk) = uplu(ji,jk) + yuflx(ji,jj)
            vplu(ji,jk) = vplu(ji,jk) + yvflx(ji,jj)
  304     continue
        else
          do 306 ji=1,nlon
            tplu(ji,jk) = tplu(ji,jk) + zt(ji,jk)*zmf(ji,jk)*cpair
            qplu(ji,jk) = qplu(ji,jk) + zq(ji,jk)*zmf(ji,jk)
c           comment out next two lines to turn off "cumulus friction"
            uplu(ji,jk) = uplu(ji,jk) + zu(ji,jk)*zmf(ji,jk)
            vplu(ji,jk) = vplu(ji,jk) + zv(ji,jk)*zmf(ji,jk)

            cpluc(ji,jk) = cpluc(ji,jk)
     *                   + ( zlw(ji,jk) + ylwc_t(ji,jk)
     *                                  - ylwc_c(ji,jj,jk) )*zmf(ji,jk)
            cplua(ji,jk) = cplua(ji,jk) - ylwc_a(ji,jj,jk)  *zmf(ji,jk)
            cplus(ji,jk) = cplus(ji,jk) - ylwc_s(ji,jj,jk)  *zmf(ji,jk)

            zconda(ji,jk)= zconda(ji,jk) + zcond(ji,jk)
            zheata(ji,jk)= zheata(ji,jk) + zcond(ji,jk)*zlath(ji,jk)
            zpreca(ji,jk)= zpreca(ji,jk) + zprec(ji,jk)
  306     continue
        endif

        if (ntraca.le.ntracb) then
          do 3060 n=ntraca,ntracb
            do 3062 ji=1,nlon
              aplu(ji,jk,n) = aplu(ji,jk,n) + za(ji,jk,n)*zmf(ji,jk)
 3062       continue
 3060     continue
        endif

  300 continue

c        For "pbl" plumes,set history daily min/max pbl heights,
c        accumulate diagnostic number of pbl-tops in each layer(numpbl),
c        and save diagnostic zonal mean plume vertical velocities in
c        wplume0. For "free" plumes, accumulate wfre1,2 for later
c        setting (after all initiating layers done) of wplume1.

c--------------------
      if (ifpbl) then
c--------------------

        do 400 ji=1,nlon
          heimin(ji,jj) = min (heipbl(ji), heimin(ji,jj))
          heimax(ji,jj) = max (heipbl(ji), heimax(ji,jj))
          if (laypbl(ji).ne.nlevp)
     *      numpbl(laypbl(ji)) = numpbl(laypbl(ji)) + 1
  400   continue

        do 402 jk=1,nlev
          wplume0(jj,jk) = 0.
          nzw = 0
          do 404 ji=1,nlon
            if (zw(ji,jk).gt.0.) then
              wplume0(jj,jk) = wplume0(jj,jk) + zw(ji,jk)
              nzw = nzw + 1
            endif
  404     continue
          if (nzw.gt.0) wplume0(jj,jk) = wplume0(jj,jk)/nzw
  402   continue

c---------
      else
c---------

        do 412 jk=1,nlev
          do 414 ji=1,nlon
c           wfre1(ji,jk) = wfre1(ji,jk) + zarea(ji,jk)*zrho(ji,jk)
c    *                                    *zw(ji,jk)
c           wfre2(ji,jk) = wfre2(ji,jk) + zarea(ji,jk)*zrho(ji,jk)
            wfre1(ji,jk) = max (wfre1(ji,jk),zw(ji,jk))
            wfre2(ji,jk) = 1.
  414     continue
  412   continue

c----------
      endif
c----------

c        End of loop over initiating layers

c*************
  510 continue
c*************

c        Set diagnostic zonal "free" plume bottom and top locations
c        and vertical velocity wplume1

      if (.not.ifpbl) then
  
        do 520 ji=1,nlon
          if (laytop(ji).ne.nlevp)
     *      numfret(laytop(ji)) = numfret(laytop(ji)) + 1
          if (laybot(ji).ne.0) 
     *      numfreb(laybot(ji)) = numfreb(laybot(ji)) + 1
  520   continue

        do 530 jk=1,nlev
          wplume1(jj,jk) = 0.
          nzw = 0
          do 532 ji=1,nlon
            if (wfre1(ji,jk).gt.0. .and. wfre2(ji,jk).gt.0.) then
              nzw = nzw + 1
              wplume1(jj,jk)= wplume1(jj,jk) + wfre1(ji,jk)/wfre2(ji,jk)
            endif
  532     continue
          if (nzw.gt.0) wplume1(jj,jk) = wplume1(jj,jk)/nzw
  530   continue

      endif

c        End of "pbl" or "free" loop

c*************
  500 continue
c*************

c        Check CFL mass-flux  criterion, find worst violation for each
c        column, and reduce fluxes uniformly for a column to make CFL
c        ok. (But not for bottom layer, otherwise wouldn't conserve
c        surface fluxes).
 
      call zero (zcfl, nlon)

      do 600 jk=1,nlev-1
        do 602 ji=1,nlon
          zcfl(ji) = max ( zcfl(ji), 
     *                     zmfa(ji,jk)*dtplume*gravit
     *                     / min (ypthic(ji,jj,jk),ypthic(ji,jj,jk+1))
     *                   ) 
  602   continue
  600 continue
 
c     cflmax = 1.e6
      cflmax = 0.75

      do 610 jk=1,nlev-1
        do 612 ji=1,nlon
          zfac = cflmax / max (cflmax, zcfl(ji))
          zmfa(ji,jk)  = zmfa(ji,jk)  * zfac
          tplu(ji,jk)  = tplu(ji,jk)  * zfac
          qplu(ji,jk)  = qplu(ji,jk)  * zfac
          uplu(ji,jk)  = uplu(ji,jk)  * zfac
          vplu(ji,jk)  = vplu(ji,jk)  * zfac
          cpluc(ji,jk) = cpluc(ji,jk) * zfac
          cplua(ji,jk) = cplua(ji,jk) * zfac
          cplus(ji,jk) = cplus(ji,jk) * zfac
          zconda(ji,jk)= zconda(ji,jk)* zfac
          zheata(ji,jk)= zheata(ji,jk)* zfac
          zpreca(ji,jk)= zpreca(ji,jk)* zfac
  612   continue

        if (ntraca.le.ntracb) then
          do 6120 n=ntraca,ntracb
            do 6122 ji=1,nlon
              zfac = cflmax / max (cflmax, zcfl(ji))
              aplu(ji,jk,n) = aplu(ji,jk,n) * zfac
 6122       continue
 6120     continue
        endif
  610 continue 

c     zz = 0.
c     do 8880 ji=1,nlon
c       zz = zz +  cflmax / max (cflmax, zcfl(ji))
c8880 continue
c     write(*,8881) nstep, jj, ji, cflmax, zz/nlon
c8881 format('plume: zonal mean cfl: nstep,jj,cflmax,zfac=',2i4,2f6.3)

c        Compute flux convergence into each layer and add to reservoirs
c        fplum* ("flux*time per layer", in commun) for use in vdif.
c        Also increment tracer fields themselves (no reservoir).
c        Use dtplume since reservoirs are not leapfrog variables
c        (see comments in reserv).

c        For simplicity, shift fluxes at layer midpoints down to next
c        lowest interface, (so give zconda,zheata,zpreca to layer below)
 
      do 700 jk=1,nlev
        do 702 ji=1,nlon
          fplumt(ji,jj,jk) = fplumt(ji,jj,jk)
     *                     + (  tplu(ji,jk)-tplu(ji,jk-1) 
     *                        + zheata(ji,jk-1)           ) * dtplume
          fplumq(ji,jj,jk) = fplumq(ji,jj,jk)
     *                     + (  qplu(ji,jk)-qplu(ji,jk-1)
     *                        - zconda(ji,jk-1)           ) * dtplume
          fplumu(ji,jj,jk) = fplumu(ji,jj,jk)
     *                     + (  uplu(ji,jk)-uplu(ji,jk-1) ) * dtplume
          fplumv(ji,jj,jk) = fplumv(ji,jj,jk)
     *                     + (  vplu(ji,jk)-vplu(ji,jk-1) ) * dtplume

          zqtot = zqtot
     *          + (qplu(ji,jk)-qplu(ji,jk-1))*dtplume*cosbud(jj)/nlon
          zptot = zptot + zconda(ji,jk-1)*dtplume*cosbud(jj)/nlon
          if (ifpbl .and. jk.eq.nlev)
     *      zetot = zetot + yqflx(ji,jj)*dtplume*cosbud(jj)/nlon
  702   continue

        if (ntraca.le.ntracb) then
          do 7020 n=ntraca,ntracb
            do 7022 ji=1,nlon
              zm = ypthic(ji,jj,jk)/gravit
              tracer(ji,jj,jk,n) = tracer(ji,jj,jk,n)
     *                      + (aplu(ji,jk,n)-aplu(ji,jk-1,n))*dtplume/zm
 7022       continue
 7020     continue
        endif
  700 continue

c        Increment large-scale liquid water contents (ylwc_*, kg/kg)
c        set within-plume precipitation (plumeprec, kg/m2/s), increment
c        global budget terms (d[m,h]clouf, kg/m2)

      do 720 ji=1,nlon
        ifanv(ji) = .true.
  720 continue

      do 730 jk=1,nlev                                   ! top to bottom
        do 732 ji=1,nlon

          zm = ypthic(ji,jj,jk)/gravit
          zc = cpluc (ji,jk)   - cpluc(ji,jk-1)
     *       + zconda(ji,jk-1) - zpreca(ji,jk-1)
          zlwcold(ji,jk) = ylwc_c(ji,jj,jk)

c         Give to anvil if plume top is high enough, else to conv
          if (zc.ne.0.) then
            if (ifanv(ji) .and. sig(jk).lt.siganvil) then 
              ylwc_a(ji,jj,jk)= ylwc_a(ji,jj,jk) + dtplume*zc/zm
            else
              ylwc_c(ji,jj,jk)= ylwc_c(ji,jj,jk) + dtplume*zc/zm
            endif
            ifanv(ji) = .false.
          endif

          ylwc_a(ji,jj,jk) = ylwc_a(ji,jj,jk)
     *                     + (cplua(ji,jk)-cplua(ji,jk-1)) * dtplume/zm
          ylwc_s(ji,jj,jk) = ylwc_s(ji,jj,jk)
     *                     + (cplus(ji,jk)-cplus(ji,jk-1)) * dtplume/zm

          plumeprec(ji,jj,jk) = zpreca(ji,jk-1)

          dmclouf = dmclouf + zconda(ji,jk-1)*(cosbud(jj)/nlon)*dtplume
          dhcloud = dhcloud + zheata(ji,jk-1)*(cosbud(jj)/nlon)*dtplume
  732   continue
  730 continue

c        Crudely fix very slightly negative cloud liquid water contents
c        (otherwise causes blowups in cldcmp/radctl)

      call zero (ztotc1, nlon)
      call zero (ztota1, nlon)
      call zero (ztots1, nlon)
      call zero (ztotc2, nlon)
      call zero (ztota2, nlon)
      call zero (ztots2, nlon)
      do 750 jk=1,nlev
        do 752 ji=1,nlon
          zm = ypthic(ji,jj,jk)/gravit
          ztotc1(ji) = ztotc1(ji) + ylwc_c(ji,jj,jk)*zm
          ztota1(ji) = ztota1(ji) + ylwc_a(ji,jj,jk)*zm
          ztots1(ji) = ztots1(ji) + ylwc_s(ji,jj,jk)*zm
          ylwc_c(ji,jj,jk) = max (ylwc_c(ji,jj,jk), 0.) 
          ylwc_a(ji,jj,jk) = max (ylwc_a(ji,jj,jk), 0.) 
          ylwc_s(ji,jj,jk) = max (ylwc_s(ji,jj,jk), 0.) 
          ztotc2(ji) = ztotc2(ji) + ylwc_c(ji,jj,jk)*zm
          ztota2(ji) = ztota2(ji) + ylwc_a(ji,jj,jk)*zm
          ztots2(ji) = ztots2(ji) + ylwc_s(ji,jj,jk)*zm
  752   continue
  750 continue

      do 754 ji=1,nlon
c       Since z*2 > z*1, if z*1>0, so is z*2:
        ztotc1(ji) = cvmgt ( ztotc1(ji)/max(ztotc2(ji),1.e-20), 0.,
     *                       ztotc1(ji).gt.0. )
        ztota1(ji) = cvmgt ( ztota1(ji)/max(ztota2(ji),1.e-20), 0.,
     *                       ztota1(ji).gt.0.)
        ztots1(ji) = cvmgt ( ztots1(ji)/max(ztots2(ji),1.e-20), 0.,
     *                       ztots1(ji).gt.0.)
  754 continue

      do 756 jk=1,nlev
        do 758 ji=1,nlon
          ylwc_c(ji,jj,jk) = ztotc1(ji) * ylwc_c(ji,jj,jk)
          ylwc_a(ji,jj,jk) = ztota1(ji) * ylwc_a(ji,jj,jk)
          ylwc_s(ji,jj,jk) = ztots1(ji) * ylwc_s(ji,jj,jk)
  758   continue
  756 continue

c*****
c     Diagnostic printout: cloud liquid water
c     do 8882 ji=1,nlon
c       do 8883 jk=1,nlev
c         if (ylwc_c(ji,jj,jk).lt.0.) then
c           write (*,*) 'plume: ji=',ji,' jj=',jj,' jk=',jk
c           write (*,8887)
c           do 8884 k=1,nlev 
c             zm = ypthic(ji,jj,k)/gravit
c             zca = (cpluc(ji,k)-cpluc(ji,k-1))*dtplume/zm
c             zcb = zconda(ji,k-1)*dtplume/zm
c             zcc = -zpreca(ji,k-1)*dtplume/zm
c             write(*,8888) k, 
c    *                      zw(ji,k),
c    *                      zmf(ji,k)*dtplume*gravit/ypthic(ji,jj,k),
c    *                      ylwc_t(ji,k)*1.e3,
c    *                      (zlw(ji,k)+ylwc_t(ji,k))*1.e3,
c    *                      zlwcold(ji,k)*1.e3,
c    *                      ylwc_c(ji,jj,k)*1.e3,
c    *                      zca*1.e3,
c    *                      zcb*1.e3,
c    *                      zcc*1.e3
c8887         format(' k', 2x, '         zw', '        zmf',
c    *                        '      ylwc_t', '        zlw',
c    *                        '    ylw_cold', '    ylw_cnew',
c    *                        '     d(cplu)',
c    *                        '   d(zconda)', '   d(zpreca)')
c8888         format(i2, 2x, 9f12.8)
c8884       continue
c           call endrun
c         endif
c8883   continue
c8882 continue
c*****

c     iuprinz(1) = nout
      iuprinz(1) = 6
      iuprinz(2) = iuout

c*****
c     Diagnostic printout: vertical profile at one point
c     nz = nint(86400./dtime)
c     if (nstep.gt.nstop-nz .and. jj.eq.norec/2) then
c       zzmax = -1.e20
c       do 8000 ji=1,nlon
c         zz = 0.
c         do 8002 jk=1,nlev
c           if (zconda(ji,jk-1).gt.0.) zz = zz + 1.
c8002     continue
c         if (zz.gt.zzmax) then
c           jim =ji
c           zzmax = zz
c         endif
c8000   continue

c       if (jj.eq.14) then
c         jim = 18
        if (jj.eq.1) then                                          !cloc
          jim = 1                                                  !cloc
          do 8010 ipz=1,nprinz
c           write(iuprinz(ipz),8012) nstep*dtime/86400.,jim,jj
c8012       format(/' day=',f8.3,'  ji=',i3,'  jj=',i3)
            write(iuprinz(ipz),*)
            write(iuprinz(ipz),8014) 
     *       (jk, sigkmh(jk),
     *            zarea(jim,jk), zw(jim,jk),
     *            zt(jim,jk),    zq(jim,jk), zrho(jim,jk),
     *            yh(jim,jk),    yt(jim,jj,jk),
     *            yq(jim,jj,jk), yrho(jim,jk),
     *            tplu(jim,jk),  qplu(jim,jk)*.864e5,
     *            zconda(jim,jk)*.864e5, 
     *            fplumt(jim,jj,jk)/(cpair*ypthic(jim,jj,jk)/gravit),
     *            fplumq(jim,jj,jk)/(ypthic(jim,jj,jk)/gravit),
c    *      jk=nlev,1,-1)
     *      jk=1,nlev)
 8014       format(' jk',4x,'sig',
     *             5x,'zarea',8x,'zw',
     *             8x,'zt',8x,'zq',6x,'zrho',
     *             8x,'yh',8x,'yt',
     *             8x,'yq',6x,'yrho',
     *             6x,'tplu',6x,'qplu',
     *             4x,'zconda',
     *             4x,'fplumt',4x,'fplumq'
     *           /(i3,f7.3,14f10.5))
 8010     continue
        endif
c     endif
c*****

c*****
c     Diagnostic printout: one line for one point
c     nz = nint(86400./dtime)
c     if (nstep.gt.nstop-2*nz .and. nstep.le.nstop-nz.and.jj.eq.norec/2)
c       then
c       if (mod(nstep,nz).eq.1) write(iuout,8020) ji,jj
c       ji = nlon/2
      if (jj.eq.1) then                                            !cloc
        ji = 1                                                     !cloc
c 8020   format(/' ji,jj=',2i3) !EWS - not used
        do 8022 jk=1,nlev
          chaplu(jk) = '   '
 8022   continue
        chaplu(laytop(ji)) = '***'
 
        do 8030 ipz=1,nprinz
c         write(iuprinz(ipz),8032)
c8032     format(5x,'day',2x,'heipbl',3x,'ytflx',
c    *           5x,5(6x,'yh',1x),5x,5(6x,'yq',1x) )
c         write(iuprinz(ipz),8034) nstep*dtime/86400., 
c    *      heipbl(ji), ytflx(ji,jj),
c    *      (yh(ji,jk),chaplu(jk)(1:1), 
c    *      jk=nlev,nlev-4,-1),
c    *      (yq(ji,jj,jk),chaplu(jk)(1:1), jk=nlev,nlev-4,-1)
c8034     format(f8.3,f8.1,f8.1,5x,5(f8.2,a1),5x,5(f8.5,a1))
          write(iuprinz(ipz),'(a,f8.1)') 'heipbl=',heipbl(ji)
 8030    continue
      endif
c*****

c        End of overall loop over latitude

c*************
 1000 continue
c*************

c*****
c     Diagnostic printout: number of negative q's globally vs level
      zz = dtime / (3600.)
      qnegmax = 1.e20
      nqnegtot = 0
      do 8040 jk=1,nlev
        qrat(jk) = 0.
        qneg(jk) = 0.
        nqrat(jk) = 0
        nqneg(jk) = 0
        do 8042 jj=1,norec
        do 8042 ji=1,nlon
          dqnew = zz * fplumq(ji,jj,jk) / (ypthic(ji,jj,jk)/gravit)
          if (abs(yq(ji,jj,jk)).gt.0.00005) then
            qrat(jk) = qrat(jk) + abs(dqnew) / abs(yq(ji,jj,jk))
            nqrat(jk) = nqrat(jk) + 1
          endif
          qnew = yq(ji,jj,jk) + dqnew
          if (qnew.lt.0.) then
            qneg(jk) = qneg(jk) + qnew
            nqneg(jk) = nqneg(jk) + 1
            nqnegtot = nqnegtot + 1
            if (qnew.lt.qnegmax) then
              qnegmax = qnew
              jkmax = jk
              jimax = ji
              jjmax = jj
            endif
          endif
 8042   continue
        if (nqrat(jk).gt.0) qrat(jk) = qrat(jk)/nqrat(jk)
        if( nqneg(jk).gt.0) qneg(jk) = qneg(jk)/nqneg(jk)
 8040 continue
 
      if (nqnegtot.gt.0) then
        do 8050 ipz=1,nprinz
          write(iuprinz(ipz),8052)
     *       nstep*dtime/86400., zptot, zetot, zqtot,
     *       qnegmax, jimax,jjmax,jkmax,
     *       nqrat(jkmax),qrat(jkmax),nqneg(jkmax),qneg(jkmax)
 8052     format(/' day=',f8.3,
     *            '   zptot=',f10.5,'   zetot=',f10.5,'   zqtot=',f10.5
     *           /' qnegmax=',f10.5,'  at=',3i4
     *           /' nqrat=',i4,'   qrat=',f9.6,'  nqneg=',i4,
     *            '  qneg=',f9.6)
 8050   continue
      endif
c*****

c*****
c     Diagnostic printout: max heimax, numfret
c     nz = nint(86400./dtime)
c     if (nstep.gt.nstop-nz) then
c       zz = -1.e20
c       do 8060 jj=1,norec
c       do 8060 ji=1,nlon
c         if (heimax(ji,jj).gt.zz) then
c           izz = ji
c           jzz = jj
c           zz = heimax(ji,jj)
c         endif
c8060   continue
c       write(*,8062) nstep*dtime/86400., zz, izz, jzz,
c    *                (numfret(jk),jk=nlev,1,-1) 
c8062   format(' day=',f8.3,'  heimax=',f8.1,'  at ji,jj=',2i3,
c    *         '  numfret=',(18i5))
c     endif
c*****
      
c*****
c     Diagnostic printout: number of pbl and plume tops in each layer
c     (written to fort.87) 
      if (mod(nstep,nint(1.*86400./dtime)) .eq. 0) then
        iu = 87
        write(iu,8070) nstep*dtime/86400.
        write(iu,8072) 'pbl  ', (numpbl (jk),jk=nlev,1,-1)
        write(iu,8072) 'fret ', (numfret(jk),jk=nlev,1,-1)
        write(iu,8072) 'freb ', (numfreb(jk),jk=nlev,1,-1)
        CALL flush(iu)
 8070   format('day: ',f8.3)
 8072   format(4x,'num',a,'=',(18i5))
      endif
c*****

c     DO J=1,ND
c       IF (ZMF(1,J) .NE. 0.0) THEN
c         YT(1,1,J) = YT(1,1,J) + TPLU(1,J) / (ZMF(1,J)*CPAIR)
c         YQ(1,1,J) = YQ(1,1,J) + QPLU(1,J) / (ZMF(1,J))
c       END IF
c     END DO


      return
      end
