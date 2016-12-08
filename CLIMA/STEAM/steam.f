      PROGRAM RIDGES
      PARAMETER (NT=6, NP=6, NAD=11)
C
C      This program calculates the heat transport properties of water
C   in the vicinity of the midocean ridges.  The viscosity formulation
C   is from Haar and Gallagher, except that the coefficients a3 and
C   b14 are corrected to the values given in Sengers and Kamgar-Parsi
C   (J. Chem. Phys. Ref. Data, 13, 185-205, 1984).
C
C      The essence of this program comes from appendix B 
C   of the NBS/NRC steam tables,  by L. Haar, J.S.Gallagher, 
C   and G.S. Kell, Hemisphere Publishing Corporation, 
C   Washington, D.C., 1984.
C
C   Internal units are temperatures in degrees K and densities in
C   g/cm3.  Pressures are in mpa (1 mpa = 10 bars = 10**6 dynes/cm2);
C   free energies are in j/g; and the sound velocity is in m/sec.
C
      COMMON /qqqq/ q0, q5
      COMMON /fcts/ ad, gd, sd, ud, hd, cvd, cpd, dpdt, 
     1 dvdt, cjtt, cjth
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, 
     2 uref, sref
      COMMON /nconst/ g(40), ii(40), jj(40), nc
C
      DIMENSION TC(NT),PBAR(NP),RHO(NT,NP),CPS(NT,NP),ALPHAS(NT,NP),
     2  ETAKS(NT,NP),RAT1(NT,NP),RAT2(NT,NP),RAT3(NT,NP),SATP(NT)
      DIMENSION TAD(NAD),PAD(NAD),RHOAD(NAD),CPAD(NAD),ALPHAD(NAD),
     2  ETAKAD(NAD),RAT1AD(NAD),RAT2AD(NAD),RAT3AD(NAD),ZAD(NAD)
C
      DATA SATP/NT*0./
      DATA TTOP/374.0/, PTOP/220.55/, ZAD(1)/0./
C
C   Find the following thermodynamic quantities as functions of 
C   temperature (degrees C) and pressure (bar):
C       D = density (g/cm3)
C       CP = specific heat at constant pressure (J/g-K)
C       ALPHAV = cubical expansion coefficient, i.e. (1/V)dV/dt
C                at constant P (1/K)
C       ETAK = kinematic viscosity (cm2/s)
C
C   Also find the ratios:
C       RAT1 = ALPHAV/ETAK
C       RAT2 = ALPHAV*D*CP/ETAK
C       RAT3 = (D*CP)**(-1/3)
C
C   DELT = width of temperature interval
C   DELP = width of pressure interval
      DELT = 1000./(NT - 1)
      DELP = 1000./(NP - 1)
      GRAV = 980.
C
      OPEN(UNIT=1,FILE='STEAM.OUT',STATUS='NEW')
      write(1,4999)
 4999 format(2x,'I',6x,'CP',7X,'ALPHAV',5X,'D',8X,'ETAK',6X,'RAT2')

C
C ***** Skip for now *****
!      GO TO 3500
      DO 10 J=1,NP
  10  PBAR(J) = (J - 1)*DELP
      PBAR(1) = 0.01
      TC(1) = 0.
C
C   Loop over temperature
      DO 1 I=5,NT
      TC(I) = (I - 1)*DELT
      print 1000, i,tc(i)
 1000 format(1x,'I =',i2,2x,'TC =',f5.0)
      IF (TC(I) .LT. 1.) TC(I) = 1.
      T = TC(I) + 273.15
      CALL BB(T)
      RT = GASCON * T
C
C   Calculate approximate saturation vapor pressure (in bars)
      IF (TC(I) .GT. 374.) GO TO 11
      ISAT = I
      SATP(I) = 10. * VP(T)
  11  CONTINUE
C
C   Loop over pressure
      DO 2 J=1,NP
      print 1001, j,pbar(j)
 1001 format(6x,'J =',i2,2x,'PBAR =',f5.0)
      P = 0.1*PBAR(J)
C
      DGSS = P/T/.4
      PSAT = 2.E4
      DL = 0.
      DV = 0.
      IF(T.LT.TZ) CALL PCORR(T,PSAT,DL,DV)
      IF(P.GT.PSAT) DGSS = DL
      CALL DFIND(D,P,DGSS,T,DQ)
C
      CALL QQ(T,D)
      AB = BASE(D,T)
      CALL THERM(D,T)
      CALL VISCOS(D,T,ETA,ETAK)
C
C   Four main thermodynamic functions
C      A = AD*RT
C      U = UD*RT
C      H = HD*RT
C      GIBBS = GD*RT
C
C      S = SD*GASCON
C      CV = CVD*GASCON
C      VL = 1./D
C      DPDD = DQ
C      DPDT1 = DPDT
C
      CP = CPD*GASCON
      CPCGS = 1.E7*CP
      ALPHAV = D * DVDT
      ETKCGS = 1.E4 * ETAK
C
C   Store values in matrix format
      RHO(I,J) = D
      CPS(I,J) = CPCGS
      ALPHAS(I,J) = ALPHAV
      ETAKS(I,J) = ETKCGS
      RAT1(I,J) = ALPHAV/ETKCGS
      RAT2(I,J) = ALPHAV*D*CPCGS/ETKCGS
      RAT3(I,J) = (D*CPCGS)**(-0.33333333)
   2  CONTINUE
   1  CONTINUE
 3500 Continue
C   End loops over temperature and pressure
C
C   Calculate adiabat extending downward from the critical point
      NAD1 = NAD - 1
      PRINT 400
 400  FORMAT(1X,'Adiabat calculation:'/'  Enter 0 to start from'
     2  ' critical point',/'        1 to start from elsewhere')
      READ (5,401) ISTART
 401  FORMAT(I1)
      IF (ISTART .EQ. 0) GO TO 41
C
      PRINT 402
 402  FORMAT(1X,'Enter PTOP (bars), TTOP (deg C)',/1X,
     2  'XXX.X XXX.X')
      READ (5,403) PTOP,TTOP
 403  FORMAT(F5.1,1X,F5.1)
C
  41  CONTINUE
      PAD(1) = PTOP
      TAD(1) = TTOP
      DZAD = 0.1
C
      DO 20 I=1,NAD1
      T = TAD(I) + 273.15
      P = 0.1 * PAD(I)
      CALL BB(T)
      RT = GASCON * T
      DGSS = P/T/.4
      PSAT = 2.E4
      DL = 0.
      DV = 0.
      IF(T.LT.TZ) CALL PCORR(T,PSAT,DL,DV)
      IF(P.GT.PSAT) DGSS = DL
      CALL DFIND(D,P,DGSS,T,DQ)
C
      CALL QQ(T,D)
      AB = BASE(D,T)
      CALL THERM(D,T)
      CALL VISCOS(D,T,ETA,ETAK)
C
      CP = CPD*GASCON
      CPCGS = 1.E7*CP
      BETAS = - (CJTT - 1./D)/CP * 0.1
      ALPHAV = D * DVDT
      ETKCGS = 1.E4 * ETAK
      print 5000, i,cp,betas,zad(i),d
 5000 format(1x,'i =',i2,2x,'cp =',1pe10.3,2x,'betas =',e10.3,2x,
     2  'zad =',e10.3,2x,'d =',e10.3)
      RHOAD(I) = D
      CPAD(I) = CP
      ALPHAD(I) = ALPHAV
      ETAKAD(I) = ETKCGS
      RAT1AD(I) = ALPHAV/ETKCGS
      RAT2AD(I) = ALPHAV*D*CPCGS/ETKCGS
      RAT3AD(I) = (D*CPCGS)**(-0.33333333)
      write(1,5001) i,cp,alphav,d,etkcgs,rat2ad(i)
 5001 format(1x,i2,2x,1p5e10.3)
C
      ZAD(I+1) = ZAD(I) + DZAD
      DP = 0.1*D*GRAV*DZAD
      PAD(I+1) = PAD(I) + DP
  20  TAD(I+1) = TAD(I) + DP*BETAS
C   End adiabat calculation
C
C   Calculate average quantities over adiabat
      RHOAV = 0.
      CPAV = 0.
      ALPAV = 0.
      ETAKAV = 0.
      RAT1AV = 0.
      RAT2AV = 0.
      RAT3AV = 0.
C
      DO 30 I=1,NAD1
      RHOAV = RHOAV + RHOAD(I)/NAD1
      CPAV = CPAV + CPAD(I)/NAD1
      ALPAV = ALPAV + ALPHAD(I)/NAD1
      ETAKAV = ETAKAV + ETAKAD(I)/NAD1
      RAT1AV = RAT1AV + RAT1AD(I)/NAD1
      RAT2AV = RAT2AV + RAT2AD(I)/NAD1
      RAT3AV = RAT3AV + RAT3AD(I)/NAD1
  30  CONTINUE
C
      PRINT 95
  95  FORMAT(/1X,'ZAD, TAD, PAD')
      PRINT 200, ZAD
 200  FORMAT(1X,1P7E10.3)
      PRINT 200, TAD
      PRINT 200, PAD
      PRINT 299
 299  FORMAT(/1X,'AVERAGE PROPERTIES ALONG ADIABAT')
      PRINT 300, RHOAV,CPAV,ALPAV,ETAKAV,RAT1AV,
     2  RAT2AV,RAT3AV
 300  FORMAT(1X,'RHOAV =',1PE10.3,2X,'CPAV =',E10.3,2X,
     2  'ALPAV =',E10.3,2X,'ETAKAV =',E10.3,/1X,'RAT1AV =',
     3  E10.3,3X,'RAT2AV =',E10.3,2X,'RAT3AV =',E10.3)
      PRINT 301, PTOP,TTOP
 301  FORMAT(1X,'PTOP =',1PE10.3,2X,'TTOP =',E10.3)
C
      TC(1) = 0.
      PBAR(1) = 0.
C
C   Don't print any of this for now
      STOP
      PRINT 99
  99  FORMAT(/1X,'DENSITY'/)
      PRINT 100, TC
 100  FORMAT(11X,6F10.0)
      DO 3 J=1,NP
   3  PRINT 101, PBAR(J),(RHO(I,J),I=1,NT)
 101  FORMAT(1X,F10.0,1P6E10.3)
C
      PRINT 98
  98  FORMAT(/1X,'KINEMATIC VISCOSITY'/)
      PRINT 100, TC
      DO 4 J=1,NP
   4  PRINT 101, PBAR(J),(ETAKS(I,J),I=1,NT)
C
      PRINT 97
  97  FORMAT(/1X,'ALPHAV'/)
      PRINT 100, TC
      DO 5 J=1,NP
   5  PRINT 101, PBAR(J),(ALPHAS(I,J),I=1,NT)
C
      PRINT 96
  96  FORMAT(/1X,'SATP'/)
      PRINT 100, TC
      PRINT 101, SATP
C
      STOP
      END

      BLOCK data
C      This blockdata subroutine supplies most of the fixed parameters
C   used in the rest of the routines.
C
      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /nconst/ g(40), ii(40), jj(40), nc
      COMMON /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
      COMMON /bconst/ bp(10), bq(10)
      COMMON /addcon/ atz(4), adz(4), aat(4), aad(4)

      DATA atz/2*64.e1, 641.6e0, 27.e1/, adz/3*.319e0, 1.55e0/, 
     1 aat/2*2.e4, 4.e4, 25.e0/, aad/34.e0, 4.e1, 3.e1, 1.05e3/
      DATA wm/18.0152e0/, gascon/.461522e0/, tz/647.073e0/, 
     1 aa/1.e0/, nc/36/
      DATA uref, sref/-4328.455039e0, 7.6180802e0/
      DATA g1, g2, gf/11.e0, 44.33333333333e0, 3.5e0/
      DATA bp/.7478629e0, -.3540782e0, 2*0., .7159876e-2, 0., 
     1  -.3528426e-2, 3*0./, bq/1.1278334e0, 0., -.5944001e0, 
     2 -5.010996e0, 0., .63684256e0, 4*0./
      DATA g/-.53062968529023e3, .22744901424408e4, .78779333020687e3, 
     1  -.69830527374994e2, .17863832875422e5, -.39514731563338e5, 
     2  .33803884280753e5, -.13855050202703e5, -.25637436613260e6, 
     3  .48212575981415e6, -.34183016969660e6, .12223156417448e6, 
     4  .11797433655832e7, -.21734810110373e7, .10829952168620e7, 
     5  -.25441998064049e6, -.31377774947767e7, .52911910757704e7, 
     6  -.13802577177877e7, -.25109914369001e6, .46561826115608e7, 
     7  -.72752773275387e7, .41774246148294e6, .14016358244614e7, 
     8  -.31555231392127e7, .47929666384584e7, .40912664781209e6, 
     9  -.13626369388386e7, .69625220862664e6, -.10834900096447e7, 
     a  -.22722827401688e6, .38365486000660e6, .68833257944332e4, 
     b  .21757245522644e5, -.26627944829770e4, -.70730418082074e5, 
     c  -.225e0, -1.68e0,  .055e0, -93.0e0/
      DATA ii/4*0, 4*1, 4*2, 4*3, 4*4, 4*5, 4*6, 4*8, 2*2,0,4,3*2,4/
      DATA jj/2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 
     1 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 1, 3*4, 0, 2, 0, 0/
      END

      SUBROUTINE pcorr(t, p, dl, dv)
C  Calculates the vapor pressure p and the liq and
C  vapor densities corresponding to the input t, corrected such that
C  gl-gv=0.  The function vp is required which will give a reasonably
C  good approximation to the vapor pressure to be used as the starting
C  point for the iteration.
C
      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      p = VP(t)
   2  call CORR(t, p, dl, dv, delg)
      dp = delg*gascon*t/(1./dv - 1./dl)
      p = p + dp
      IF( ABS(delg) .LT. 1.e-4) RETURN
      GO TO 2
      END

      FUNCTION vp(t)
C   This function calculates an approximation to the vapor pressure, vp,
C   as a function of the input temperature.  The vapor pressure
C   calculated agrees with the vapor pressure predicted by the surface
C   to within .02% to within a degree or so of the critical temperature,
C   and can serve as an initial guess for further refinement by
C   imposing the condition that g1 = gv.
c
      DIMENSION a(8)
      DATA a/-7.8889166e0,  2.5514255e0,  -6.716169e0, 
     1  33.239495e0,  -105.38479e0,  174.35319e0,  -148.39348e0, 
     2  48.631602e0/
      IF( t .GT. 314.) GO TO 2
      pl = 6.3573118 - 8858.843/t + 607.56335 * t**(-.6)
      vp = .1*EXP(pl)
      RETURN
C
   2  v = t/647.25
      w = ABS(1.-v)
      b = 0.
      DO 4 i=1, 8
      z = i
   4  b = b + a(i) * w**((z+1.)/2.)
      q = b/v
      vp = 22.093 * EXP(q)
      RETURN
      END

      SUBROUTINE corr(t, p, dl, dv, delg)
C  Calculates, for an input t and p at or near the
C  vapor pressure, the corresponding liquid and vapor densities and also
C  delg = (gl-gv)/rt for use in calculating the correction to the vapor
C  pressure for delg = 0.
C
      COMMON /qqqq/ q00, q11
      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
      COMMON /fcts/ ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, 
     1 cjtt, cjth

      rt = gascon * t
      IF( t .GT. 646.3) GO TO 101
      dliq = dl
      IF( dl .LE. 0.) dliq = 1.11 - 0.0004*t
      call BB(t)
      call DFIND(dl, p, dliq, t, dq)
      call THERM(dl, t)
      gl = gd
      dvap = dv
      IF( dv .LE. 0.) dvap = p/rt
      call DFIND(dv, p, dvap, t, dq)
      IF( dv .LT. 5.e-7) dv = 5.e-7
      call THERM(dv, t)
      gv = gd
      delg = gl - gv
      RETURN
C
 101  p = 0.
      IF( t .GT. 647.126) RETURN
      delg = 0.
      call BB(t)
      tau = .657128 * (1. - t/647.126)**.325
      dl = .322 + tau
      dv = .322 - tau
      zdum = BASE(dv, t)
      call QQ(t, dv)
      p = rt*dv*zb + q00
      RETURN
      END

      SUBROUTINE bb(t)
C       This subroutine calculate the b's of eqs. 3,4 using coefficients
C   from blockdata, calculating also the first and second derivs w.r.
C   to temp.  The b's calculated here are in cm3/g.
C
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
      COMMON /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
      COMMON /bconst/ bp(10), bq(10)
      DIMENSION v(10)

      v(1) = 1.

      DO 2 i=2, 10
   2  v(i) = v(i-1)*tz/t

      b1 = bp(1) + bp(2) * ALOG(1./v(2))
      b2 = bq(1)
      b1t = bp(2) * v(2)/tz
      b2t = 0.
      b1tt = 0.
      b2tt = 0.

      DO 4 i=3, 10
      fi2 = i - 2
      b1 = b1 + bp(i)*v(i-1)
      b2 = b2 + bq(i)*v(i-1)
      b1t = b1t - fi2*bp(i)*v(i-1)/t
      b2t = b2t - fi2*bq(i)*v(i-1)/t
      b1tt = b1tt + bp(i)*fi2*fi2*v(i-1)/t/t
   4  b2tt = b2tt + bq(i)*fi2*fi2*v(i-1)/t/t

      b1tt = b1tt - b1t/t
      b2tt = b2tt - b2t/t
      RETURN
      END

      SUBROUTINE dfind(dout, p, d, t, dpd)
C   Routine to find density corresponding to input pressure p(mpa), and
C   temperature t(k), using initial guess density d(g/cm3).  The output
C   density is in g/cm3, also, as a byproduct, dpdrho is calculated
C   ("dpd", mpa cm3/g)
C
      COMMON /qqqq/ q0, q5
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref

      dd = d
      rt = gascon*t
      IF( dd .LE. 0.) dd = 1.e-8
      IF( dd .GT. 1.9) dd = 1.9
      l = 0
   9  l = l + 1
  11  IF( dd .LE. 0.) dd = 1.e-8
      IF( dd .GT. 1.9) dd = 1.9
      call QQ(t,  dd)
      pp = rt * dd * BASE(dd, t) + q0
      dpd = rt*(z + y*dz) + q5
C   The following 3 lines check for negative dp/drho, and if so assume
C   guess to be in 2-phase region, and adjust guess accordingly.
C
      IF( dpd .GT. 0.) GO TO 13
      IF( d .GE. .2967) dd = dd*1.02
      IF( d .LT. .2967) dd = dd*.98
      IF( l .LE. 10) GO TO 9
  13  dpdx = dpd*1.1
      IF( dpdx .LT. .1) dpdx = .1
      dp = ABS(1. - pp/p)
      IF( dp .LT. 1.e-8) GO TO 20
      IF( d .GT. .3 .AND. dp .LT. 1.e-7) GO TO 20
      IF( d .GT. .7 .AND. dp .LT. 1.e-6) GO TO 20
      x = (p-pp)/dpdx
      IF( ABS(x).GT. .1) x = x*.1/ABS(x)
      dd = dd + x
      IF( dd .LE. 0.) dd = 1.e-8
  19  IF( l .LE. 30) GO TO 9
  20  Continue
      dout = dd
      RETURN
      END

      SUBROUTINE qq(t, d)
C   This routine calculates, for a given t(k) and d(g/cm3), the residual
C   contributions to:  pressure (q), helmholtz fct (ar),dpdrho (q5),
C   and also to the gibbs function, entropy, internal energy, enthalpy,
C   isochoric heat capacity, and dpdt.  (eq.5)
C       Terms 37 thru 39 are the additional terms affecting only the
C   immediate vicinity of the critical point, and term 40 is the
C   additional term improving the low t, high p region.
C
      COMMON /resf/ ar, gr, sr, ur, hr, cvr, dpdtr
      COMMON /qqqq/ q, q5
      DIMENSION qr(11), qt(10), qzr(9), qzt(9)
      EQUIVALENCE (qr(3), qzr(1)), (qt(2), qzt(1))
      COMMON /nconst/ g(40), ii(40), jj(40), n
      COMMON /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
      COMMON /addcon/ atz(4), adz(4), aat(4), aad(4)

      rt = gascon*t
      qr(1) = 0.
      q5 = 0.
      q = 0.
      ar = 0.
      dadt = 0.
      cvr = 0.
      dpdtr = 0.
      e = EXP(-aa*d)
      q10 = d*d*e
      q20 = 1. - e
      qr(2) = q10
      v = tz/t
      qt(1) = t/tz

      DO 4 i=2, 10
      qr(i+1) = qr(i)*q20
   4  qt(i) = qt(i-1)*v

      DO 10 i=1, n
      k = ii(i) + 1
      l = jj(i)
      zz = k
      qp = g(i)*aa*qzr(k-1)*qzt(l)
      q = q + qp
      q5 = q5 + aa*(2./d - aa*(1.-e*(k-1)/q20))*qp
      ar = ar + g(i)*qzr(k)*qzt(l)/q10/zz/rt
      dfdt = q20**k * (1-l)*qzt(l+1)/tz/k
      d2f = l*dfdt
      dpt = dfdt*q10*aa*k/q20
      dadt = dadt + g(i)*dfdt
      dpdtr = dpdtr + g(i)*dpt
  10  cvr = cvr + g(i)*d2f/gascon

      qp = 0.
      q2a = 0.
c
      DO 20 j=37, 40
      IF( g(j) .EQ. 0.) GO TO 20
      k = ii(j)
      km = jj(j)
      ddz = adz(j-36)
      del = d/ddz - 1.
      IF( ABS(del) .LT. 1.e-10) del = 1.e-10
      dd = del*del
      ex1 = -aad(j-36)*del**k
      dex = EXP(ex1)*del**km
      att = aat(j-36)
      tx = atz(j-36)
      tau = t/tx - 1.
      ex2 = -att*tau*tau
      tex = EXP(ex2)
      q10 = dex*tex
      qm = km/del - k*aad(j-36)*del**(k-1)
      fct = qm*d*d*q10/ddz
      q5t = fct*(2./d + qm/ddz) - (d/ddz)*(d/ddz)*q10*(km/del/del +
     1  k*(k-1)*aad(j-36)*del**(k-2))
      q5 = q5 + q5t*g(j)
      qp = qp + g(j)*fct
      dadt = dadt - 2.*g(j)*att*tau*q10/tx
      dpdtr = dpdtr - 2.*g(j)*att*tau*fct/tx
      q2a = q2a + t*g(j)*(4.*att*ex2 + 2.*att)*q10/tx/tx
      ar = ar + q10*g(j)/rt
  20  Continue

      sr = - dadt/gascon
      ur = ar + sr
      cvr = cvr + q2a/gascon
      q = q + qp
      RETURN
      END

      FUNCTION base(d, t)
C   Calculates base portion of Helmholtz equation
C
      COMMON /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
      COMMON /basef/ ab, gb, sb, ub, hb, cvb, dpdtb
      COMMON /aconst/ wm, gascon, tz, a, z, dz, y, uref, sref

      y = 0.25*b1*d
      y2 = y*y
      x = 1. - y
      x2 = x*x
      x3 = x*x2
      x4 = x*x3
      t2 = t*t
      z0 = (1. + g1*y + g2*y2)/x3
      z = z0 + 4.*y*(b2/b1 - gf)
      dz0 = (g1 + 2.*g2*y)/x3 + 3.*(1. + g1*y + g2*y2)/x4
      dz = dz0 + 4.*(b2/b1 - gf)
      ab = -ALOG(x) - (g2 - 1.)/x + 28.16666667/x2 + 4.*y*(b2/b1 - gf)
     1  + 15.16666667 + ALOG(d*t*gascon/.101325)
      gb = ab + z
      base = z
      bb2tt = t2*b2tt
      ub = -t*b1t*(z-1.-d*b2)/b1 - d*t*b2t
      hb = z + ub
      tbt2 = (t*b1t/b1)*(t*b1t/b1)
      cvb = 2.*ub + (z0 - 1.)*(tbt2 - t2*b1tt/b1)
     1  - d*(bb2tt - gf*b1tt*t2) - tbt2 *y*dz0
      dpdtb = base/t + base*d/z*(dz*b1t/4. + b2t - b2/b1*b1t)
      sb = ub - ab
      RETURN
      END

      SUBROUTINE therm(d, t)
C   This subroutine calculates the thermodynamic functions in
C   dimensionless units (ad=a/rt, gd=g/rt, sd=s/r, ud=u/rt,
C   hd=h/rt, cvd=cv/r, and cpd=cp/r)

      COMMON /aconst/ wm, gascon, tz, aa, zb, dzb, y, uref, sref
      COMMON /qqqq/ qp, qdp
      COMMON /basef/ ab, gb, sb, ub, hb, cvb, dpdtb
      COMMON /resf/ ar, gr, sr, ur, hr, cvr, dpdtr
      COMMON /idf/ ai, gi, si, ui, hi, cvi, cpi
      COMMON /fcts/ ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, 
     1 cjtt, cjth

      call IDEAL(t)
      rt = gascon*t
      z = zb + qp/rt/d
      dpdd = rt*(zb + y*dzb) + qdp
      ad = ab + ar + ai - uref/t + sref
      gd = ad + z
      ud = ub + ur + ui - uref/t
      dpdt = rt*d*dpdtb + dpdtr
      cvd = cvb + cvr + cvi
      cpd = cvd + t*dpdt*dpdt/(d*d*dpdd*gascon)
      hd = ud + z
      sd = sb + sr + si - sref
      dvdt = dpdt/dpdd/d/d
      cjtt = 1./d - t*dvdt
      cjth = - cjtt/cpd/gascon
      RETURN
      END

      SUBROUTINE ideal(t)
C   This subroutine calculates the thermodynamic properties for 
C   water in the ideal gas state from function of H.W. Woolley
C
      COMMON /idf/ ai, gi, si, ui, hi, cvi, cpi
      DIMENSION c(18)
      DATA c /.19730271018e2, .209662681977e2, -.483429455355e0, 
     1  .605743189245e1, 22.56023885e0, -9.87532442e0, 
     2  -.43135538513e1, .458155781e0, -.47754901883e-1, 
     3  .41238460633e-2, -.27929052852e-3, .14481695261e-4, 
     4  -.56473658748e-6, .16200446e-7, -.3303822796e-9, 
     5  .451916067368e-11, -.370734122708e-13, .137546068238e-15/

      tt = t/1.e2
      tl = ALOG(tt)
      gi = - (c(1)/tt + c(2)) * tl
      hi = (c(2) + c(1)*(1.-tl)/tt)
      cpi = c(2) - c(1)/tt

      DO 8 i=3, 18
      tti6 = tt**(i-6)
      gi = gi - c(i)*tti6
      hi = hi + c(i)*(i-6)*tti6
   8  cpi = cpi + c(i)*(i-6)*(i-5)*tti6

      ai = gi - 1.
      ui = hi - 1.
      cvi = cpi - 1.
      si = ui - ai
      RETURN
      END

      SUBROUTINE VISCOS(D,T,ETA,ETAK)
      DIMENSION A(4),B(6,5),sav(6,5)
C
      DATA A/0.0181583, 0.0177624, 0.0105287, -0.0036744/
      DATA B/
     9  0.501938,  0.162888,  -0.130356, 0.907919,-0.551119,  0.146543,
     1  0.235622,  0.789393,   0.673665, 1.207552, 0.0670665,-0.0843370,
     2 -0.274637, -0.743539,  -0.959456,-0.687343,-0.497089,  0.195286,
     3  0.145831,  0.263129,   0.347247, 0.213486, 0.100754, -0.032932,
     4 -0.0270448,-0.0253093,-0.0267758,-0.0822904,0.0602253,-0.0202595/
C
      DMKS = 1000. * D
      DSTAR = 317.763
      TSTAR = 647.27
      SUM1 = 0.
      SUM2 = 0.
C
      DO 1 K=1,4
      K1 = K - 1
   1  SUM1 = SUM1 + A(K)*(TSTAR/T)**K1
      ETA0 = SQRT(T/TSTAR) * 1.E-6/SUM1
C
      DO 2 I=1,6
      I1 = I - 1
      DO 2 J=1,5
      J1 = J - 1
   2  SUM2 = SUM2 + B(I,J)*(TSTAR/T - 1.)**I1 * (DMKS/DSTAR - 1.)**J1
      ETA = ETA0 * EXP(DMKS/DSTAR * SUM2)
C
C   Convert to kinematic viscosity
      ETAK = ETA/DMKS
      RETURN
      END
C  /EOF 
