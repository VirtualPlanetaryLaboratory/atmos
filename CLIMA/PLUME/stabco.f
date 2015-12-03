      subroutine stabco (yt, yq, ypm, ypthic, fplumt, fplumq, ylwc_s)

c     Does stable condensation, iterating niter times 

c     yt       = temperature field (modified)
c     yq       = specific humidity (modified)
c     ypm      = mid-layer pressure (supplied)
c     ypthic   = layer pressure thickness (supplied)
c     fplumt,q = reservoirs of "flux convergence*time" (3-D) (modified)
c     ylwc_s   = stratiform liquid water content (3-D) (modified)

c----------------------------------------------------------------------
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      INCLUDE 'complume'
c----------------------------------------------------------------------
      dimension 
     *  yt(nlon,norec,nlev),      yq(nlon,norec,nlev),
     *  ypm(nlon,norec,nlev),     ypthic(nlon,norec,nlev),
     *  fplumt(nlon,norec,nlev),  fplumq(nlon,norec,nlev),
     *  ylwc_s(nlon,norec,nlev)

      dimension
     *  zt(nlon),      zq(nlon),
     *  delq(nlon),    za(nlon),
     *  zsvp(nlon),    zsvq(nlon)

      parameter (niter=1)
c----------------------------------------------------------------------
      INCLUDE 'comsatplume'
c----------------------------------------------------------------------


c        Overall loops over latitude and level

c=======================
      do 1000 jj=1,norec
c=======================
      do 1100 jk=1,nlev
c=======================

c          Copy t,q to temporary arrays for this longitude circle
c          (so don't change coml30 arrays for time n)

        call scopy (nlon, yt(1,jj,jk), 1, zt, 1)
        call scopy (nlon, yq(1,jj,jk), 1, zq, 1)

c          Iterate for this longitude circle

c---------------------------
        do 2000 iter=1,niter
c---------------------------

c         call estabv (zsvp, zsvq, zt, ypm(1,jj,jk), nlon)
          do 190 ji=1,nlon
            zsvp(ji) = esat (zt(ji))
            zsvq(ji) = qsat (zsvp(ji), ypm(ji,jj,jk))
  190     continue

          do 200 ji=1,nlon

c              za is latent heat of condensation, determined by ambient
c              temperature (ie, rain or snow), plus correction
c              for specific heats to be consistent with lsx (imagining
c              all phase changes occur at tmelt)

c           za(ji) = cvmgt ( latvap + cpwv*(zt(ji)-tmelt)
c    *                              - ch2o*(zt(ji)-tmelt),
c    *                       latsub + cpwv*(zt(ji)-tmelt)
c    *                              - cice*(zt(ji)-tmelt),
c    *                       zt(ji).ge.tmelt )
c           For now, all cloud water is liquid with no specific heat
            za(ji) = latvap + cpwv*(zt(ji)-tmelt)

c              Test if exceed zqmax, which ranges from alw*saturation
c              for zstrat=0 and 100%saturation for zstrat=1, where
c              zstrat = existing stratus cloud fractional area (as
c              in cldcmp). This assumes uniformly distributed specific
c              humidities in non-stratus areas with range of +/- 1.-alw,
c              and saturation in stratus areas. Condense amount that 
c              yields zqmax, allowing fopr linearized change in svp with
c              temperature.Store condensed amount in delq.
c              (alw is a namelist param).

            zstrat = min (1., max (0., ylwc_s(ji,jj,jk)/clw_s ))
c 777        zqmax = zsvq(ji) * ( zstrat  + (1.-zstrat)*alw ) !EWS - label not used
c777        zqmax = zsvq(ji) * alw  ! (early v2)

            if (zq(ji).gt.zqmax) then
              dz = (za(ji)/rh2o) * zqmax / (zt(ji)**2)
              cpz = cpair*(1.+cpvir*zq(ji))
              delq(ji) = (zq(ji) - zqmax) / (1.+dz*za(ji)/cpz)
              zq(ji) = zq(ji) - delq(ji)
              zt(ji) = zt(ji) + delq(ji)*za(ji)
     *                                   / (cpair*(1.+cpvir*zq(ji)))
            else
              delq(ji) = 0.
            endif

c              Accumulate reservoirs fplum[t,q] (flux convergences*time
c              per layer). Use factor of 0.5 since reservoirs are not
c              leapfrog variables (see comments in reserv). Also
c              increment stratiform liquid water content (Kg/Kg) and
c              increment global budget terms (d[m,h]clouf)

            zz = delq(ji) * ypthic(ji,jj,jk)/gravit
            fplumt(ji,jj,jk) = fplumt(ji,jj,jk) + zz*za(ji)*0.5
            fplumq(ji,jj,jk) = fplumq(ji,jj,jk) - zz*0.5
            ylwc_s(ji,jj,jk) = ylwc_s(ji,jj,jk) + delq(ji)*0.5

            dmclouf = dmclouf + zz       *(cosbud(jj)/nlon)*0.5
            dhcloud = dhcloud + zz*za(ji)*(cosbud(jj)/nlon)*0.5
  200     continue

c---------------
 2000   continue
c---------------

c=============
 1100 continue
c=============

c=============
 1000 continue
c=============

      return
      end
