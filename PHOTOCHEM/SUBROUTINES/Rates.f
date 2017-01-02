      SUBROUTINE RATES
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      CHARACTER*8 REACTYPE,PLANET,CHEMJ
      real*8 mass
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'
      real*8 k0,kinf

c-mc rate constant units are cm^3/mol/s
       ! chemical reaction file
       open(9, file='PHOTOCHEM/INPUTFILES/reactions.rx',status='OLD')


667   FORMAT(58X,E9.2,3X,F8.2)            !for two body reaction rates
668   FORMAT(58X,E9.2,3X,E9.2,2X,2F5.2)   !for three body reaction rates

669   FORMAT(70X,E9.3,2X,E10.3,2X,E10.3,2X,E9.3,2X,E9.3,2X,E9.3,2X,
     $      E10.3,2X,E10.3,2X,E9.3,2X,E9.3) !two, three, & capture rate reaction rates

      KIDA = 0

       do J=1,NR
             if (kida.eq.1) then
                read(9,*)
           ! Read in two body reaction rates
                  if (REACTYPE(J) .EQ. '2BODY') then
           !read in Arhenius and Temperature factor (note TFAC contains the negative sign. not standard practice)
                   read (9,669) alpha_0,beta_0,gamma_0,f_0,g_0

                     do i=1,nz
                      A(J,I) = alpha_0 * EXP (-1.*beta_0/T(I))   !two body reaction rates
                     enddo
           ! Read in 2BCRT reaction rates
                  else if (REACTYPE(J) .EQ. '2BCRT') then
           !read in Arhenius and Temperature factor (note TFAC contains the negative sign. not standard practice)
                     read (9,669) alpha_0,beta_0,gamma_0,f_0,g_0
                      do i=1,nz
                       A(J,I) = alpha_0 * EXP (beta_0/T(I))
                      enddo
           ! Read in three body reaction rates
                  else if (REACTYPE(J) .EQ. '3BODY') then        
                    read(9,669) alpha_0,beta_0,gamma_0,f_0,g_0,alpha_inf
     $        ,beta_inf,gamma_inf,f_inf,g_inf


                      if (PLANET .EQ. 'MARS') then
                        B=B*2.5      !multiply low density rate by 2.5 to account for CO2 rather than N2 as background gas (Nair, 94)
                      endif

                      do i=1,nz
                       A(J,I)= TBDY(B,C,D,E,T(I),DEN(I))  !computed three body reaction rate
                      enddo
                end if

             else   ! if kida <> 1 then proceed with the old school reactions.rx read ins.
C READ IN TWO BODY REACTION RATES
          if (REACTYPE(J) .EQ. '2BODY') then
           !read in Arhenius and Temperature factor (note TFAC contains the negative sign. not standard practice)
           read (9,667) ARH, TFAC


          do i=1,nz
            A(J,I) = ARH * EXP (TFAC/T(I))   !two body reaction rates
           enddo

C READ IN THREE BODY REACTION RATES
          else if (REACTYPE(J) .EQ. '3BODY') then
             read(9,668) B,C,D,E          !read in K0, Kinf, T0exp, Tinfexp
                                          !A(J,I) = TBDY(K0,KINF,T0exp,Tinfexp,T,DEN)
c             print *, J, B,C,D,E
             if (PLANET .EQ. 'MARS') then
                B=B*2.5      !multiply low density rate by 2.5 to account for CO2 rather than N2 as background gas (Nair, 94)
             endif
             do i=1,nz
              A(J,I)= TBDY(B,C,D,E,T(I),DEN(I))  !computed three body reaction rate
             enddo

C SPECIFY 'WEIRD' REACTION RATES
C  FILL UP WEIRD RATE CONSTANTS
c-mc the below are rate constants which don't fit in the 2BODY or 3BODY category


         else if (REACTYPE(J) .EQ. 'WEIRD') then
          read(9,*) !move to the next line in the .rx file
          DO I=1,NZ

!   H2 + O -> OH + H
      if (CHEMJ(1,J) .EQ. 'H2' .AND. CHEMJ(2,J) .EQ. 'O') THEN
       A(J,I) = 1.34E-15*(T(I)/298.)**6.52 *EXP(-1460./T(I))   ! Low T fit NIST 05
      endif

!   HO2 + HO2 -> H2O2 + O2
      if (CHEMJ(1,J).EQ.'HO2'.AND.CHEMJ(2,J).EQ.'HO2') THEN
       !JPL-02
       A(J,I) = 2.3E-13*EXP(590./T(I)) + 1.7E-33*EXP(1000./T(I))*DEN(I)
      endif

!   O + O + M -> O2 + M
      if (CHEMJ(1,J).EQ.'O'.AND.CHEMJ(2,J).EQ.'O') THEN
         JOO_O2=J  !this is used below in S+S->S2
         !note - re-hard coding into S+S->S2 evil I know...
       A(J,I) = 9.46E-34*EXP(480./T(I))*DEN(I)  ! NIST 05 low Temp in N2.  Its 2x bigger in H2 or O2
      endif

!   CO + OH -> CO2 + H (also CS + HS -> CS2 + H)  !SORG
      if ((CHEMJ(1,J).EQ.'CO'.AND.CHEMJ(2,J).EQ.'OH') .OR.
     $    (CHEMJ(1,J).EQ.'CS'.AND.CHEMJ(2,J).EQ.'HS')) THEN
       PATM = DEN(I)*1.38E-16*T(I)/1.013E6   !a good pressure
       A(J,I) = 1.5E-13 * (1. + 0.6*PATM)     !JPL-02
      endif

!   CO + O + M -> CO2 + M
      if (CHEMJ(1,J).EQ.'CO'.AND.CHEMJ(2,J).EQ.'O') THEN
       A(J,I) = 2.2E-33*EXP(-1780./T(I))*DEN(I)  ! I use NIST 05 for 257-277 K
      endif

!gna - updated according to d-g 2011
!   H + CO + M -> HCO + M (also H + CS + M -> HCS + M) !sorg
      if ((CHEMJ(1,J).EQ.'H'.AND.CHEMJ(2,J).EQ.'CO') .OR.
     $    (CHEMJ(1,J).EQ.'H'.AND.CHEMJ(2,J).EQ.'CS')) THEN
       A(J,I) = 2.0E-33*EXP(-850./T(I))*DEN(I)  ! I use NIST 05 for 333-1000 K, theory
      endif

!   H2CO + H -> H2 + HCO
      if (CHEMJ(1,J).EQ.'H2CO'.AND.CHEMJ(2,J).EQ.'H') THEN
       A(J,I) = 2.14E-12*(T(i)/298.)**1.62*EXP(-1090./T(I))    ! NIST 2005 Baulch et al 2002
      endif

!   H + H + M -> H2 + M
      if (CHEMJ(1,J).EQ.'H'.AND.CHEMJ(2,J).EQ.'H') THEN
       A(J,I) = 8.85E-33*(T(I)/287)**(-0.6) * DEN(I)  !gna Baluch 1994
      endif

!   H + OH + M -> H2O + M
      if (CHEMJ(1,J).EQ.'H'.AND.CHEMJ(2,J).EQ.'OH') THEN
       A(J,I) = 6.9E-31*(298./T(I))**2 *DEN(I)     ! Baulch et al 2002 in N2
      endif

!   CH3 + CH3 + M  ->  C2H6 + M
      if (CHEMJ(1,J).EQ.'CH3'.AND.CHEMJ(2,J).EQ.'CH3') THEN
       A71_3 = 1.17e-25*exp(-500./T(I))
       if (PLANET .EQ. 'MARS') A71_3=A71_3*2.5    !CO2 rather than N2
       A(J,I) = TBDY(A71_3,3.0E-11,3.75E0,1.0E0,T(I),DEN(I))   ! what a mess NIST 2005 -
c      A(J,I) = 1.7E-17/T(I)**2.3 * DEN(I)
      endif


!   CH3 + H2CO  ->  CH4 + HCO
      if (CHEMJ(1,J).EQ.'CH3'.AND.CHEMJ(2,J).EQ.'H2CO') THEN
       A(J,I) = 4.9E-15*(T(i)/298.)**4.4 *EXP(-2450./T(I))    !NIST 2005
      endif

!   H + NO + M  ->  HNO + M
      if (CHEMJ(1,J).EQ.'H'.AND.CHEMJ(2,J).EQ.'NO') THEN
       A77_3 = 1.2E-31*EXP(-210./T(I))                       ! there was a bug here - one that mattered
       A(J,I) = TBDY(A77_3,2.4E-10, 1.17E0,0.41E0,T(I),DEN(I))   ! taken from NIST 2005 Tsang/Hampson 91
c      A(J,I) = 1.2E-31*(298./T(I))**1.17 *EXP(-210./T(I)) * DEN(I)
      endif

!   N + N + M  ->  N2 + M
      if (CHEMJ(1,J).EQ.'N'.AND.CHEMJ(2,J).EQ.'N') THEN
C       A(J,I) = 8.3E-34*EXP(500./T(I)) * DEN(I)
c       A(J,I) = 1.25E-32                   !NIST 2005
       A(J,I) = 1.25E-32* DEN(I)       !NIST 2005 BUG HERE - SHOULD BE X DEN (KEEP AS WEIRD)
      endif


!   HNO3 + OH  ->  H2O + NO3
      if (CHEMJ(1,J).EQ.'HNO3'.AND.CHEMJ(2,J).EQ.'OH') THEN
        AK0 = 2.4E-14*EXP(460./T(I))
        AK2 = 2.7E-17*EXP(2199./T(I))
        AK3M = 6.5E-34*EXP(1335./T(I))*DEN(I)
       A(J,I) = AK0 + AK3M/(1. + AK3M/AK2)  !JPL-06
      endif

!   H + HNO  ->  H2 + NO
      if (CHEMJ(1,J).EQ.'H'.AND.CHEMJ(2,J).EQ.'HNO') THEN
       A(J,I) = 3.0E-11 * EXP(-500./T(I))   ! NIST 2005  ??? shouldn't this be 2body? (check on rx list)
      endif

!   C2H6 + O  ->  C2H5 + OH
      if (CHEMJ(1,J).EQ.'C2H6'.AND.CHEMJ(2,J).EQ.'O') THEN
       A(J,I) = 8.54E-12*(T(I)/300.)**1.5 *EXP(-2920./T(I))   ! NIST 05
      endif

!   SO + O -> SO2
      if (CHEMJ(1,J).EQ.'SO'.AND.CHEMJ(2,J).EQ.'O') THEN
       A(J,I) = 6.0E-31 * DEN(I)         ! NIST 2005 updated from D-G 2011 gna
      endif

!    S + S -> S2
      if (CHEMJ(1,J).EQ.'S'.AND.CHEMJ(2,J).EQ.'S') THEN
c      A(J,I) = 1.2e-29 * DEN(I)   ! in H2S, but its much 1e4 slower in Ar
c       A(J,I) = min(5.E-11, 3* A(JOO_O2,I))      ! reported rate is 3X larger for S+S in Ar than for O+O in Ar

       A(J,I) = 1.87E-33 * EXP(-206/T(I))*DEN(I) !updated from D-G 2011 gna
      endif

!   S + S2 + M -> S3 + M
      if (CHEMJ(1,J).EQ.'S'.AND.CHEMJ(2,J).EQ.'S2') THEN
       A(J,I) = MIN(5.0E-11,2.5E-30*DEN(I)/5.E1)    ! NIST reported In CO2 with factor 5 error
      endif

!    S2 + S2 + M -> S4 + M
      if (CHEMJ(1,J).EQ.'S2'.AND.CHEMJ(2,J).EQ.'S2') THEN
       A(J,I) = MIN(5.0E-11,2.5E-30*DEN(I)/1.E1)    ! NIST reported In CO2 with factor 5 error I'm taking lower bound
      endif

!     S + S3 + M -> S4 + M
      if (CHEMJ(1,J).EQ.'S'.AND.CHEMJ(2,J).EQ.'S3') THEN
       A(J,I) = MIN(5.0E-11,2.5E-30*DEN(I)/5.E1)    ! NIST reported In CO2 with factor 5 error
       !assumed rate equal to the S+S2 rate
      endif

!     S4 + S4 + M -> S8 + M
!       OR
!     S4 + S4 -> S8AER
      if (CHEMJ(1,J).EQ.'S4'.AND.CHEMJ(2,J).EQ.'S4') THEN
       A(J,I) = MIN(5.0E-11,2.5E-30*DEN(I)/1.E1)    ! NIST reported In CO2 with factor 5 error I'm taking lower bound
       !assumed rate equal to S2+S2
      endif

!     S + CO + M -> OCS + M
      if (CHEMJ(1,J).EQ.'S'.AND.CHEMJ(2,J).EQ.'CO') THEN
       A(J,I) = 1.*2.2E-33*EXP(-1780./T(I))*DEN(I)  ! no information
c  I'm guessing that S + CO goes at about the same rate as O + CO
c  with the same activation energy and a 3X bigger A factor
c-mc but yet factor of 1 out from, so this is the same as o+co

       !gna - DG 2011 has:
       !A(J,I) = 6.5E-33*EXP(-2180./T(I))*DEN(I)  ! assumed same as (CO+O)

      endif

!    OCS + S + M -> OCS2 + M   ! NIST 8.3e-33*Den in Ar
      if (CHEMJ(1,J).EQ.'OCS'.AND.CHEMJ(2,J).EQ.'S') THEN
       A(J,I) = 8.3E-33*Den(i)     ! reduce by 1000 to turn it off
      endif


!    CO + O1D + M -> CO2 + M
      if (CHEMJ(1,J).EQ.'CO'.AND.CHEMJ(2,J).EQ.'O1D') THEN
       A(J,I) = 0.0 * 8.00E-11           ! reported but no sign of density dependence
      endif

!    O1D + N2 + M -> N2O + M   !this is really damn slow...
      if (CHEMJ(1,J).EQ.'O1D'.AND.CHEMJ(2,J).EQ.'N2') THEN
       A(J,I) = 2.80E-36*DEN(I)*(T(I)/300.)**(-0.9) !JPL-06  no high pressure limit reported
      endif

!   CL + O2 + M -> CLOO + M    !IUPAC 2007 rec - only a low density rate (could modify TBDY to deal with this...)
      if (CHEMJ(1,J).EQ.'CL'.AND.CHEMJ(2,J).EQ.'O2') THEN
       A(J,I) = 1.4E-33*DEN(I)*(T(I)/300.)**(-3.9)    !measured in N2
      endif


!   CL + NO + M -> CLNO + M
      if (CHEMJ(1,J).EQ.'CL'.AND.CHEMJ(2,J).EQ.'NO') THEN
       A(J,I) = 7.6E-32 * (300./T(I))**1.8 * DEN(I) !JPL-06 (no high pressure limit given)
      endif

!   CCL3NO4 + M -> CCL3O2 + NO2               WEIRD     !IUPAC-07
      if (CHEMJ(1,J).EQ.'CCL3NO4'.AND.CHEMJ(2,J).EQ.'M') THEN
       k0 =4.3E-3*exp(-10235/T(I))*DEN(I)
       kinf=4.8E16*exp(-11820/T(I))
       Fc=0.32
       A(J,I) = ((k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2)))/DEN(I)
c       A(J,I) = 4.8e16*exp(-11820/T(I))     !old
      endif

!   NO + CH3O -> HNO + H2CO              WEIRD     !IUPAC worksheet - 2 body with the 300/T form... (yuk 123)
      if (CHEMJ(1,J).EQ.'NO'.AND.CHEMJ(2,J).EQ.'CH3O') THEN
       A(J,I) = 2.3E-12 * (300./T(I))**0.7
      endif

!   NO2 + CH3O + M -> CH3ONO2 + M  IUPAC
      if (CHEMJ(1,J).EQ.'NO2'.AND.CHEMJ(2,J).EQ.'CH3O') THEN
       k0 =(8.1E-29*(300./T(I))**4.5)*DEN(I)
       kinf=2.1E-11
       Fc=0.44
       A(J,I) = ((k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2)))/DEN(I)
      endif


!   NO2 + CH3O2 + M -> CH3O2NO2 +M                   WEIRD     ! IUAPC-06 (yuk 131)
      if (CHEMJ(1,J).EQ.'NO2'.AND.CHEMJ(2,J).EQ.'CH3O2') THEN
       k0 =(2.5E-30*(300./T(I))**5.5)*DEN(I)
       kinf=1.8E-11
       Fc=0.36
       A(J,I) = ((k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2)))/DEN(I)
      endif

!   CH3O2NO2  M -> CH3O2 + NO2               WEIRD     !IUPAC
      if (CHEMJ(1,J).EQ.'CH3O2NO2'.AND.CHEMJ(2,J).EQ.'M') THEN
       k0 =9.0E-5*exp(-9690./T(I))*DEN(I)
       kinf=1.1E16*exp(-10560./T(I))
       Fc=0.6
       A(J,I) = ((k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2)))/DEN(I)
      endif


!   CH2O2 + M -> CO + H2O               WEIRD     ! IUPAC-06 thermal decomp
      if (CHEMJ(1,J).EQ.'CH2O2'.AND.CHEMJ(2,J).EQ.'M') THEN
c first order rate constant of 6e4
       A(J,I)=6E4/Den(i) !because it will be multiplied by DEN(i) so net is 6e4*[CH2O2]
c      A(J,I)=9.96E-20 !old version
      endif


!   CH2OOH + M -> OH + H2CO              WEIRD     !Vaghjinai et al. 1989 (via NIST) as 1st order -(yuk 278) has 1e-10 as 2body
      if (CHEMJ(1,J).EQ.'CH2OOH'.AND.CHEMJ(2,J).EQ.'M') THEN
       A(J,I)=5E4/Den(i) !1st order rate, will be multiplied by DEN(i) so net is 5e4*[CH2OOH]
c       A(J,I) = 1.e-10   !yuks rate
      endif


!    N2O5 + M -> NO3 + NO2               WEIRD   !IUPAC
      if (CHEMJ(1,J).EQ.'N2O5'.AND.CHEMJ(2,J).EQ.'M') THEN
       k0 =1.3E-3*exp(-11000./T(I))*(300./T(I))**3.5*DEN(I)
       kinf=9.7E14*exp(-11080./T(I))*(T(I)/300.)**0.1
       Fc=0.35
       A(J,I) = ((k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2)))/DEN(I)
      endif


!   HO2NO2 + M -> HO2 + NO2
      if (CHEMJ(1,J).EQ.'HO2NO2'.AND.CHEMJ(2,J).EQ.'M') THEN
       k0 =4.1E-5*exp(-10650./T(I))*DEN(I)
       kinf=4.8E15*exp(-11170./T(I))     !there was an error here. should be -11170
       Fc=0.6                             !this should be Fc=0.6
       A(J,I) = ((k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2)))/DEN(I)
      endif

!   CL2O2 + M -> CLO + CLO
      if (CHEMJ(1,J).EQ.'CL2O2'.AND.CHEMJ(2,J).EQ.'M') THEN
       k0 =3.7E-7*exp(-7690./T(I))*DEN(I)
       kinf=7.9E15*exp(-8820./T(I))
       Fc=0.45
c       A(J,I) = (k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
c     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2))

       A(J,I) = ((k0*kinf)/(k0 + kinf)* 10**(log10(Fc)/
     $  (1 +(log10(k0/kinf)/(0.75-1.27*log10(Fc)))**2)))/DEN(I)
      endif


!   CH3ONO + M -> CH3O + NO                WEIRD     1.06e15     -18282.                ! Fernandez-Ramos et al. 1998 (300K-1500K) 1st order reaction
      if (CHEMJ(1,J).EQ.'CH3ONO'.AND.CHEMJ(2,J).EQ.'M') THEN
       A(J,I) = 1.06e15*exp(-18041./T(I))!/DEN(I)  !or is it first order?
c       A(J,I) = 1.06e15*exp(-18282./T(I))  !old
      endif


!   CLO + O2 + M ->  CLO3 + M !  WEIRD     9.00E-28  (214 yuk) - THIS IS NOT IN JPL-06 ???
      if (CHEMJ(1,J).EQ.'CLO'.AND.CHEMJ(2,J).EQ.'O2'.AND.
     $    CHEMJ(3,J).EQ.'CLO3') THEN
       A362_0=9.00E-28*T(i)**(-2.0)
       A362_inf=4.50E-7*T(i)**(-2.0)

!tuning
       A362_0=9.00E-38*T(i)**(-2.0)
       A362_inf=4.50E-19*T(i)**(-2.0)

!rego - from DeMore 1990
c      A362_0=1.00D-32*T(i)**(-2.0)
c      A362_inf=5.00D-12*T(i)**(-2.0)



c       A(J,I) = TBDY(A362_0,A362_inf, 0.0E0,0.0E0,T(I),DEN(I))   ! from yuk (214)

c temptest
       A(J,I) = 0.0   !zeroing for now. eventually rebuild this as CLO.O2

c       A(J,I) = TBDY(1.0D-32,5.0E-12, 2.0E0,2.0E0,T(I),DEN(I))   ! from yuk (214)
      endif

!   CLO + CLO3 + M -> CL2O4 +M              WEIRD     !Xu and Lin 2003
!   CLO + CLO3  -> CLOO + OCLO              WEIRD     !Xu and Lin 2003 (slow)
!   CLO + CLO3  -> OCLO +OCLO               WEIRD     !Xu and Lin 2003 (slow)
      if (CHEMJ(1,J).EQ.'CLO'.AND.CHEMJ(2,J).EQ.'CLO3') THEN
        if (CHEMJ(3,J).EQ.'CL2O4') then 
         A0=8.62E15*exp(-1826./T(I))*T(i)**(-9.75)
         Ainf=1.43E-10*exp(-82./T(I))*T(i)**(0.094)
         A(J,I) = TBDY(A0,Ainf, 0.E0,0.E0,T(I),DEN(I))
        endif
        if (CHEMJ(3,J).EQ.'CLOO') then 
         A(J,I) = 1.85E-18*T(I)**2.28*exp(-2417./T(I))
        endif
        if (CHEMJ(3,J).EQ.'OCLO') then 
         A(J,I) = 1.42E-18*T(I)**2.11*exp(-2870./T(I))
        endif
      endif


!    O3 + CL + M -> CLO3                        WEIRD     !Simonaitis 1975
      if (CHEMJ(1,J).EQ.'O3'.AND.CHEMJ(2,J).EQ.'CL') THEN
       A(J,I) =  3E-30*Den(i)     ! rate constant estimated for 300K (some words about t-dependence - could look to analogs...)
c       A(J,I) =  1E-31*Den(i)     !lowering this uncertain rate by a factor of 30 1e-13 is nominial case
      endif


c

!    CLO + O2 -> CLOOO                       WEIRD     ! DeMore 1990 rate (Via Shindell)
      if (CHEMJ(1,J).EQ.'CLO'.AND.CHEMJ(2,J).EQ.'O2'.AND.
     $    CHEMJ(3,J).EQ.'CLOOO') THEN
!from DeMore 1990      !could take this to classic 3body
       A366_0=1.00E-36    !going for -3
       A366_inf=5.00E-12!
       A(J,I) = TBDY(A366_0,A366_inf, 2.0E0,2.0E0,T(I),DEN(I))
      endif

!    CLOOO + M -> CLO + O2                WEIRD - thermal decomp (divide above by eq rate in Prassad) - Vogel has 3.1e-18 at 200K

      if (CHEMJ(1,J).EQ.'CLOOO'.AND.CHEMJ(2,J).EQ.'M') THEN
       A367_EQ = 2.9E-26*exp(3500./T(I))

       A(J,I) = TBDY(A366_0,A366_inf, 2.0E0,2.0E0,T(I),DEN(I))
     $           /A367_EQ/1E5/DEN(I)
!    note use of coefficients from formation reaction so this is a dependency
! this species was used a test in Catling et al. 2010 for the CLO.O2 adduct - see sensitivity analysis there...
      endif

!   O + CLO + M  ->  OCLO + M
      if (CHEMJ(1,J).EQ.'O'.AND.CHEMJ(2,J).EQ.'CLO') THEN
         Alow=8.60E-21*T(I)**(-4.1)*EXP(-420./T(I))
         Ainf=4.33E-11*T(I)**(-0.03)*EXP(43./T(I))
        
       A(J,I) = TBDY(Alow,Ainf,0.0E0,0.0E0,T(I),DEN(I))   ! taken from Zhu and Lin 2003
c       print *, I, A(J,I),Alow,Ainf,T(I),DEN(I)
      endif


!   OH + CLO3 + M  ->  HCLO4 or HO2 + OCLO  (Zhu and Lin 2001) - note this will interfere if we go back to Simonitatis!
      if (CHEMJ(1,J).EQ.'OH'.AND.CHEMJ(2,J).EQ.'CLO3') THEN
         if (CHEMJ(3,J).EQ.'HCLO4') then
          A0=1.94E36*T(I)**(-15.3)*EXP(-5542./T(I))
          Ainf=3.2E-10*T(I)**(0.07)*EXP(-25./T(I))

          A(J,I) = TBDY(A0,Ainf,0.0E0,0.0E0,T(I),DEN(I))   ! taken from Zhu and Lin 2001

!should perhaps put the Simonaitis rate in here for completeness so we don't have to switch back and forth on reactions.rx

         endif

         if (CHEMJ(3,J).EQ.'HO2') then
          A(J,I)=2.1E-10*T(I)**(0.09)*EXP(-18./T(I))   !ditto
         endif
      endif



!   CS2 + S  ->  CS + S2  !SORG
      if (CHEMJ(1,J).EQ.'CS2'.AND.CHEMJ(2,J).EQ.'S') THEN
       ! Woiki et al 1995
       A(J,I) = 1.9E-14 * EXP(-580./T(I)) * (T(I)/300.)**3.97
      endif

!   C2H6S + H ->  CH3SH + CH3  !SORG

      if (CHEMJ(1,J).EQ.'C2H6S'.AND.CHEMJ(2,J).EQ.'H'.AND. !3 needed because there is another C2H6S + H  2 body reaction in the table
     $    CHEMJ(3,J).EQ.'CH3SH') THEN


!gna - possible mistake here: in D.-G. et al 2011  C2H6S + H ->  H2 + C2H4 + HS is listed with different reaction rate coefficients
!  8.34E-12    -2212.     1.6
! can't parse the Zhang paper so not sure what is correct.
       A(J,I) = 4.81E-12 * EXP(-1100./T(I)) * (T(I)/300.)**1.70    ! theory - Zhang et al. 2005
      endif

!this used to be same as  C2H6S + H ->  CH3SH + CH3  !SORG but that is not what shawn has in 2011 paper table
!   C2H6S + H ->  H2 + C2H4 + HS !HC   Theory. Zhang et al. [2005]. Produces C2H5S, which then can split into C2H4 + HS
          if(CHEMJ(1,J).EQ.'C2H6S'.AND.CHEMJ(2,J).EQ.'H'.AND.
     $    CHEMJ(3,J).EQ.'H2') THEN
             A(J,I) = 8.34E-12 * EXP(-2212./T(I)) * (T(I)/300.)**1.60 ! theory - Zhang et al. 2005
          endif


!gna (this was missing from sorg reactions.rx too so had to add it)
!C2H6S + O -> CH3 + CH3 + SO
!CH3SH + O -> CH3 + HSO
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'C2H6S'.AND.CHEMJ(2,J).EQ.'O'.AND.
     $    CHEMJ(3,J).EQ.'CH3') .OR.
     $   (CHEMJ(1,J).EQ.'C2H6SH'.AND.CHEMJ(2,J).EQ.'O'.AND.
     $    CHEMJ(3,J).EQ.'CH3') )THEN
      A(J,I) = 1.30E-11 * EXP(-410./T(I)) * (T(I)/298.)**1.1    ! Sander 2006
      endif

!gna
!C2H6S2 + O -> CH3 + CH3S + SO
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'C2H6S2'.AND.CHEMJ(2,J).EQ.'O'.AND.
     $    CHEMJ(3,J).EQ.'CH3') )THEN
      A(J,I) = 3.90E-11 * EXP(290./T(I)) * (T(I)/298.)**1.1    ! Sander 2006
      endif

!gna
!C2H6S + OH -> CH21 + CH3S + CS2O
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'C2H6S'.AND.CHEMJ(2,J).EQ.'OH'.AND.
     $    CHEMJ(3,J).EQ.'CH21') )THEN
      A(J,I) = 1.10E-11 * EXP(400./T(I)) * (T(I)/298.)**1.1    ! Sander 2006
      endif

!gna
!C2H6S2 + OH -> CH3 + CH3SH + S
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'C2H6S2'.AND.CHEMJ(2,J).EQ.'OH'.AND.
     $    CHEMJ(3,J).EQ.'CH3') )THEN
      A(J,I) = 6.00E-11 * EXP(400./T(I)) * (T(I)/298.)**1.2    ! Sander 2006
      endif

!gna
!C2H6S + O --> CH3 + Ch3 + SO
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'C2H6S2'.AND.CHEMJ(2,J).EQ.'OH'.AND.
     $    CHEMJ(3,J).EQ.'CH3') )THEN
      A(J,I) = 6.00E-11 * EXP(400./T(I)) * (T(I)/298.)**1.2    ! Sander 2006
      endif

!gna
!CH3 + OH -> CH3O + H
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'CH3'.AND.CHEMJ(2,J).EQ.'OH'.AND.
     $    CHEMJ(3,J).EQ.'CH3O') )THEN
      A(J,I) = 9.3E-11 * EXP(-1606/T(I)) * (T(I)/298.)**1.    ! Jasper 2007
      endif

!gna
!CH3 + HNO -> CH4 + NO
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'CH3'.AND.CHEMJ(2,J).EQ.'HNO'.AND.
     $    CHEMJ(3,J).EQ.'CH4') )THEN
      A(J,I) = 1.85E-11 * EXP(-176/T(I)) * (T(I)/298.)**0.6    ! Choi and Lin 2005
      endif


!gna
!H2S + H -> H2 + HS
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'H2S'.AND.CHEMJ(2,J).EQ.'H'.AND.
     $    CHEMJ(3,J).EQ.'H2') )THEN
      A(J,I) = 3.66E-12 * EXP(-455/T(I)) * (T(I)/298.)**1.94    ! Choi and Lin 2005
      endif

!gna
!SO + HCO -> HSO + CO
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'SO'.AND.CHEMJ(2,J).EQ.'HCO'.AND.
     $    CHEMJ(3,J).EQ.'HSO') )THEN
      A(J,I) = 5.6E-12 *  * (T(I)/298.)**0.4    ! Kasting 1990
      endif

!gna
!NH2 + H -> NH3
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'NH2'.AND.CHEMJ(2,J).EQ.'H'.AND.
     $    CHEMJ(3,J).EQ.'NH3') )THEN
      A(J,I) = (6.E-30*DEN(I))/(1.+3.E-20*DEN(I))    ! Gordon 1971
      endif

!gna
!NH + H -> NH2
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'NH'.AND.CHEMJ(2,J).EQ.'H'.AND.
     $    CHEMJ(3,J).EQ.'NH2') )THEN
      A(J,I) =  (6.E-30*DEN(I))/(1.+3.E-20*DEN(I))   ! Kasting 1982
      endif

!gna
!CS + HS -> CS2 + H
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'CS'.AND.CHEMJ(2,J).EQ.'CS'.AND.
     $    CHEMJ(3,J).EQ.'CS2') )THEN
      A(J,I) =  1.5E-13*(1.+0.6*DEN(I))   ! assumed samed as k(CO+OH)
      endif

!gna
!CH3SH + OH --> CH3S + H2O
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'CH3SH'.AND.CHEMJ(2,J).EQ.'OH'.AND.
     $    CHEMJ(3,J).EQ.'CH3S') )THEN
      A(J,I) = 9.90E-12 * EXP(360/T(I)) * (T(I)/298.)**1.07    ! Sander 2006
      endif


!gna
!HCO + M --> H + CO + M
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'HCO'.AND.CHEMJ(2,J).EQ.'M'.AND.
     $    CHEMJ(3,J).EQ.'H') )THEN
      A(J,I) = 6.0E-11 * EXP(-7721/T(I)) * DEN(I)   !Krasnoperov et al 2004
      endif

!gna
!HNO + M --> NO + H + M
!in reactions list as "WEIRD" but wasn't here...
      if ((CHEMJ(1,J).EQ.'HNO'.AND.CHEMJ(2,J).EQ.'M'.AND.
     $    CHEMJ(3,J).EQ.'NO') )THEN
       A(J,I) = 1.04E-6 * EXP(28618/T(I))*(T(I)/298.)**(-1.61)*DEN(I) !Tsang 1986
      endif


!gna -- ordering was wrong of CHEMJ indices compared to reactions.rx for sorg template

!   CH3S + HCS ->CS + CH3SH !SORG
      if (CHEMJ(2,J).EQ.'HCS'.AND.CHEMJ(1,J).EQ.'CH3S') THEN

         A(J,I) = 1.18E-12*EXP(-910./T(I))*(T(I)/300.)**0.65  !Liu et al. 2006 (via NIST)

      endif

!   C + H2 -> CH23
      if (CHEMJ(1,J).EQ.'C'.AND.CHEMJ(2,J).EQ.'H2') THEN
      B0 = 8.75E-31 * EXP(524./T(I))
      BI = 8.3E-11
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!   CH + H2 -> CH3            !apparantly the same as C+H2->CH23
      if (CHEMJ(1,J).EQ.'CH'.AND.CHEMJ(2,J).EQ.'H2'.AND.CHEMJ(3,J).EQ.
     $    'CH3') THEN
       B0 = 8.75E-31 * EXP(524./T(I))
       BI = 8.3E-11
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif


!   CH23 + H -> CH3
      if (CHEMJ(1,J).EQ.'CH23'.AND.CHEMJ(2,J).EQ.'H'.AND.CHEMJ(3,J).EQ.
     $    'CH3') THEN
      B0 = 3.1E-30 * EXP(457./T(I))
      BI = 1.5E-10
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  C2H + H -> C2H2
      if (CHEMJ(1,J).EQ.'C2H'.AND.CHEMJ(2,J).EQ.'H') THEN
      B0 = 1.26E-18 * EXP(-721./T(I)) / T(I)**3.1
      BI = 3.E-10
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)

      !gna - shawn has diff:
      !B0 = 2.64E-26 * EXP(-721./T(I)) / (T(I)/300.)**3.1
      !BI = 3.E-10
      !A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  CH23 + CO -> CH2CO

      if (CHEMJ(1,J).EQ.'CH23'.AND.CHEMJ(2,J).EQ.'CO') THEN
      B0 = 1.E-28
      BI = 1.E-15
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  CH3 + CO -> CH3CO
      if (CHEMJ(1,J).EQ.'CH3'.AND.CHEMJ(2,J).EQ.'CO') THEN
      A(J,I) = 1.4E-32 * EXP(-3000./T(I))*DEN(I)
      endif

!  C2H2 + H -> C2H3
      if (CHEMJ(1,J).EQ.'C2H2'.AND.CHEMJ(2,J).EQ.'H') THEN
      B0 = 2.6E-31
      BI = 3.8E-11 * EXP(-1374./T(I)) !gna shawn has 8.3E-11 in paper
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  C2H3 + CH4 ->  C2H4 + CH3
      if (CHEMJ(1,J).EQ.'C2H3'.AND.CHEMJ(2,J).EQ.'CH4') THEN
      A(J,I) = 2.4E-24 * T(I)**4.02 * EXP(-2754./T(I))
      endif

!  C2H4 + H -> C2H5
!  C3H6 + H -> C3H7
      if ((CHEMJ(1,J).EQ.'C2H4'.AND.CHEMJ(2,J).EQ.'H') .OR.
     $    (CHEMJ(1,J).EQ.'C3H6'.AND.CHEMJ(2,J).EQ.'H')) THEN
       B0 = 2.15E-29 * EXP(-349./T(I))
       BI = 4.95E-11 * EXP(-1051./T(I))
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  CH + CH4 -> C2H4 + H
      if (CHEMJ(1,J).EQ.'CH'.AND.CHEMJ(2,J).EQ.'CH4') THEN
      A270A = 2.5E-11 * EXP(200./T(I))
      A270B = 1.7E-10
      A(J,I) = AMIN1(A270A,A270B)
      endif

!  C2H5 + CH3 -> C2H4 + CH4
      if (CHEMJ(1,J).EQ.'C2H5'.AND.CHEMJ(2,J).EQ.'CH3') THEN
      A(J,I) = 3.25E-11/T(I)**0.5
      endif

!  C2H2 + OH -> CH2CO + H
      if (CHEMJ(1,J).EQ.'C2H2'.AND.CHEMJ(2,J).EQ.'OH'.AND.CHEMJ(3,J).EQ.
     $    'CH2CO') THEN
      B0 = 5.8E-31 * EXP(1258./T(I))
      BI = 1.4E-12 * EXP(388./T(I))
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  C2H5 + CH3 -> C3H8 !gna shawn has diff
      if (CHEMJ(1,J).EQ.'C2H5'.AND.CHEMJ(2,J).EQ.'CH3') THEN
      B0 = 2.519E-16 / T(I)**2.458
      BI = 8.12E-10 / T(I)**0.5
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  C3H8 + O -> C3H7 + OH
      if (CHEMJ(1,J).EQ.'C3H8'.AND.CHEMJ(2,J).EQ.'O') THEN
      A(J,I)=1.6E-11 * EXP(-2900./T(I)) + 2.2E-11 * EXP(-2200./T(I))
      endif

! C2H3 + CH3 -> C3H6
      if (CHEMJ(1,J).EQ.'C2H3'.AND.CHEMJ(2,J).EQ.'CH3') THEN
      B0 = 1.3E-22
      BI = 1.2E-10
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  CH + C2H4 -> CH2CCH2 + H                 WEIRD           ! Romani et al. [1993]
!  CH + C2H4 -> CH3C2H + H                 WEIRD                                                            ! Romani et al. [1993]
!this should be OK to capture both channels and proceed them at the same rate (which is what shawn did)
      if (CHEMJ(1,J).EQ.'CH'.AND.CHEMJ(2,J).EQ.'C2H4') THEN
       A272A = 5.5E-11 * EXP(173./T(I))
       A272B = 3.55E-10
      A(J,I) = AMIN1(A272A,A272B)
      endif

!  CH2CCH2 + H -> CH3 + C2H2              WEIRD                                                            ! Yung et al. [1984]
!  CH2CCH2 + H -> C3H5                        WEIRD                                                            ! Yung et al. [1984]
      if (CHEMJ(1,J).EQ.'CH2CCH2'.AND.CHEMJ(2,J).EQ.'H') THEN
      B0 = 8.E-24/T(I)**2 * EXP(-1225./T(I))
       if (CHEMJ(3,J).EQ.'CH3')   BI=9.7E-13 * EXP(-1550./T(I))
       if (CHEMJ(3,J).EQ.'CH3H5') BI=1.4E-11 * EXP(-1000./T(I))

      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  C2H3 + C2H5 -> CH3 + C3H5              WEIRD                                                            ! Romani et al. [1993]
      if (CHEMJ(1,J).EQ.'C2H3'.AND.CHEMJ(2,J).EQ.'C2H5'.AND.CHEMJ(3,J)
     $ .EQ.'CH3') THEN
        B0 = 1.9E-27
        BI = 2.5E-11
      A(J,I) = BI - B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif


!  C3H5 + H -> C3H6                        WEIRD                                                            ! Yung et al. [1984] !has another branch!
      if (CHEMJ(1,J).EQ.'C3H5'.AND.CHEMJ(2,J).EQ.'H'.AND.CHEMJ(3,J)
     $ .EQ.'C3H6') THEN
        B0 = 1.E-28
        BI = 1.E-11
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  CH + C2H2 -> C3H2 + H                 WEIRD                                                            ! Romani et al. [1993]
      if (CHEMJ(1,J).EQ.'CH'.AND.CHEMJ(2,J).EQ.'C2H2') THEN
      A271A = 1.75E-10 * EXP(61./T(I))
      A271B = 5.3E-10
      A(J,I) = AMIN1(A271A,A271B)
      endif

!  CH23 + C2H2 -> CH3C2H                      WEIRD                                                            ! Laufer et al. [1983] and Laufer [1981]
      if (CHEMJ(1,J).EQ.'CH23'.AND.CHEMJ(2,J).EQ.'C2H2') THEN
       B0 = 3.8E-25
       BI = 2.2E-12
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  CH3C2H + H -> CH3 + C2H2              WEIRD                                                            ! Whytock et al [1976] and Von Wagner and Zellner [1972]
!  CH3C2H + H -> C3H5                        WEIRD                                                            ! Yung et al. [1984], same as RXN 303
      if (CHEMJ(1,J).EQ.'CH3C2H'.AND.CHEMJ(2,J).EQ.'H') THEN
       B0 = 8.E-24/T(I)**2 * EXP(-1225./T(I))
       BI = 9.7E-12 * EXP(-1550./T(I))
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  C3H2 + H -> C3H3                        WEIRD                                                            ! Yung et al. [1984]
!  C3H3 + H -> CH3C2H                      WEIRD                                                            ! Yung et al. [1984], same as RXN 37
!  C3H3 + H -> CH2CCH2                     WEIRD                                                            ! Yung et al. [1984], same as RXN 337
      if ((CHEMJ(1,J).EQ.'C3H2'.AND.CHEMJ(2,J).EQ.'H').OR.
     $    (CHEMJ(1,J).EQ.'C3H3'.AND.CHEMJ(2,J).EQ.'H')) THEN
       B0 = 1.7E-26
       BI = 1.5E-10
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  CH23 + C2H2 -> CH2CCH2                     WEIRD                                                            ! Laufer et al. [1983] and Laufer [1981]
      if (CHEMJ(1,J).EQ.'CH23'.AND.CHEMJ(2,J).EQ.'C2H2') THEN
       B0 = 3.8E-25
       BI = 3.7E-12
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif

!  C2H5 + H -> C2H6                        WEIRD                                                            ! Gladstone [1983]
      if (CHEMJ(1,J).EQ.'C2H5'.AND.CHEMJ(2,J).EQ.'H'.AND.CHEMJ(3,J)
     $ .EQ.'C2H6') THEN
       B0 = 5.5E-23/T(I)**2 * EXP(-1040./T(I))
       BI = 1.5E-13 * EXP(-440./T(I))
      A(J,I) = B0*BI*DEN(I)/(B0*DEN(I) + BI)
      endif



c mc some testing detritus...
c      do i=280,180,-5.
c         print *, i,A(295,1)/(9.3E-28*exp(8835./i)),1.3E-27*exp(8744./i)  !testing CL2O2
c         print *, i,A(252,1)/(5.7E-25*exp(2500./i))*1e-9!,5.7E-25*exp(2500./i)
c      enddo
c      stop
c      print *, (A(284,I),I=1,NZ)



          ENDDO
         else
             read(9,*) !do nothing with PHOTO...
          endif
       end if
       end do
       close(9)



      RETURN
      END

      real*8 FUNCTION TBDY(A0,AI,CN,CM,T,D)
      real*8 B0,BI,Y,X,A0,AI,CN,CM,T,D
      !this form of the three body equations allows us to copy the n and m factor directly from the JPL recomendations.
      !note that it is presented there as: (T/300)^(-n) which is done here as (300/T)^(n) either way. it should be n and m directly.
      B0 = A0*(300./T)**CN
      BI = AI*(300./T)**CM
      Y = LOG10(B0*D/BI)
      X = 1./(1. + Y**2)
      TBDY = B0*D/(1. + B0*D/BI) * 0.6**X
      RETURN
      END
