      SUBROUTINE RATESWASP12B
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      CHARACTER*8 REACTYPE,PLANET,CHEMJ
      real*8 mass
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The data for Enthalpy coefficients and Shomate equation is obtained from 
! NIST chemistry webbook:
!
! http://webbook.nist.gov/chemistry/form-ser.html
!
! Whenever the coefficients are not available, corresponding compound's 
! Enthalpy formation function data was generated from NASA's thermo-build 
! web-page to generate corresponding analytical fitting formulae.
!
!  http://www.grc.nasa.gov/WWW/CEAWeb/ceaThermoBuild.htm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      real *8 Df_Ho(NSP,NZ),Df_Ho298(NSP),HH298(NSP,NZ)
      real *8 Df_Go(NSP,NZ),So(NSP,NZ),df,S_ent
      real *8 Df_HoC,HH298C(NZ),SoC(NZ)
      real *8 coeff(8),RR,kB,PFAC,akeq,akr0,akri

!CONTENT BELOW IS RELEVANT TO THE WASP12B/HOT JUPITER BACKWARD REACTIONS ONLY
!***********************************************************************
! Df_Ho - Standard Enthalpy of formation (increment in the enthalpy
!         for forming the given compound from its elements in 
!         reference states
!
! Df_Ho298 - Enthalpy of formation at the reference temperature (298.15 K)
!
! HH298 - Enthalpy of the compound (or) the elements (H-H298)
! 
! Df_Go - Gibbs free energy of formation
! So    - Entropy

! Shomate equation coefficients:
! coeff(*,1) = A, coeff(*,2) = B, coeff(*,3) = C, coeff(*,4) = D
! coeff(*,5) = E, coeff(*,6) = F, coeff(*,7) = G, coeff(*,8) = H

!      open(30,file='Gibbs-species.dat')
!      open(31,file='GibbsMajorSpecies.dat')
      RR           = 8.3144 ! Gas constant

! Enthalpies of formation at the reference temperature = 298.15 K
! Note: HE is part of NSP species, but not needed here. This is for NSP-IN species.

      Df_Ho298(LO)      = 249.18
       print*,'LO',LO
      Df_Ho298(LO2)      = 0.0
       print*,'LO2',LO2
      Df_Ho298(LH2O)    = -241.826
       print*,'LH2O',LH2O
      Df_Ho298(LH)      = 217.998
       print*,'LH',LH
      Df_Ho298(LOH)     = 38.99
       print*,'LOH',LOH
      Df_Ho298(LCO2)    = -393.52
       print*,'LCO2',LCO2
      Df_Ho298(LCO)     = -110.53 
       print*,'LCO',LCO
      Df_Ho298(LHCO)    =  43.51
       print*,'LHCO',LHCO
!      Df_Ho298(LH2CO)   = -115.90  ! From NIST database
      Df_Ho298(LH2CO)   = -108.580  ! From NASA thermobuild
       print*,'LH2CO',LH2CO
      Df_Ho298(LCH4)    = -74.87
       print*,'LCH4',LCH4
      Df_Ho298(LCH3)    = 145.69
       print*,'LCH3',LCH3
      Df_Ho298(LCH3OH)  = -200.940
       print*,'LCH3OH',LCH3OH
      Df_Ho298(LCH)     = 594.13
       print*,'LCH',LCH
      Df_Ho298(LCH23)    = 386.39
       print*,'LCH23',LCH23
      Df_Ho298(LH2)     = 0.0
       print*,'LH2',LH2
      Df_Ho298(LH2COH)  = -17.8 ! From NASA thermobuild
       print*,'LH2COH',LH2COH
      Df_Ho298(LC)      = 716.68
       print*,'LC,Df_Ho298(LC) ',LC,Df_Ho298(LC) 
      Df_Ho298(LCH3O)   = 13.0    ! From NASA thermobuild
       print*,'LCH3O',LCH3O
      Df_Ho298(LO1D)    = 0.0
       print*,'LO1D',LO1D
      Df_Ho298(LCH21)   = 417.5  ! Yung & Demore, P 145, Table 5.9
       print*,'LCH21',LCH21

      

       

      do I = 1,NZ

!*******************************************************************
! For elemental carbon (graphite)


         HH298C(I) = 9.165e-25 *T(I)**8-2.59e-20*T(I)**7 +
     .                3.101e-16*T(I)**6  -2.052e-12*T(I)**5+
     .                8.22e-09*T(I)**4-2.055e-05*T(I)**3+
     .                0.03244*T(I)**2 -5.931*T(I) -642.4

         SoC(I)    = 1.3250e-27*T(I)**8- 3.534e-23*T(I)**7 +
     .                3.906e-19*T(I)**6 - 2.297e-15*T(I)**5+
     .                7.590e-12*T(I)**4-1.292e-08*T(I)**3+
     .                4.162e-06*T(I)**2 + 0.03122*T(I) - 3.704

         HH298C(I) = HH298C(I)/1000.0e0

!*******************************************************************
! For 'O'



         HH298(LO,I) = -6214.0 + 21.1*T(I) -6.334E-05*T(I)**2
          So(LO,I) = -8.168e-19*T(I)**6+ 9.303e-15*T(I)**5
     .               -4.314e-11*T(I)**4 + 1.054e-07*T(I)**3
     .               -0.0001482*T(I)**2 + 0.1318*T(I) + 132.5

         HH298(LO,I) = HH298(LO,I)/1000.0e0
!*******************************************************************
******
! FOR O2 
! T(I)= 700-2000 K
         if( (T(I) >=700.0E0) .and. (T(I) <=2000.0E0))then
             coeff(1) = 30.03235
             coeff(2) = 8.772972
             coeff(3) = -3.988133
             coeff(4) = 0.788313
             coeff(5) = -0.741599
             coeff(6) = -11.32468
             coeff(7) = 236.1663
             coeff(8) = 0.0E0
         else if(T(I)>2000.0E0)then
!T(I) = 2000-6000 K
             coeff(1) = 20.91111
             coeff(2) = 10.72071
             coeff(3) = -2.020498
             coeff(4) = 0.146449
             coeff(5) = 9.245722
             coeff(6) = 5.337651
             coeff(7) = 237.6185
             coeff(8) = 0.0
         endif
 
         call shomate(coeff,df,S_ent,T,i,NZ)
         HH298(LO2,I) = df
         So(LO2,I)    = S_ent
!*******************************************************************
**********
! FOR H2O
!  TEMP = 500-1700 K
       if( (T(I) >=500.0E0) .and. (T(I) <=1700.0E0))then
             coeff(1) = 30.09200
             coeff(2) = 6.832514
             coeff(3) = 6.793435
             coeff(4) = -2.534480
             coeff(5) =  0.082139
             coeff(6) =  -250.8810
             coeff(7) = 223.3967
             coeff(8) = -241.8264

!   TEMP = 1700-6000 K
        else if(T(I)>1700.0E0)then
             coeff(1) = 41.96426
             coeff(2) = 8.622053
             coeff(3) = -1.499780
             coeff(4) = 0.098119
             coeff(5) = -11.15764
             coeff(6) = -272.1797
             coeff(7) = 219.7809
             coeff(8) = -241.8264
         endif

        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LH2O,I) = df
        So(LH2O,I)    = S_ent
!*******************************************************************
! FOR H

        coeff(1) = 20.78603 
        coeff(2) = 4.850638E-10
        coeff(3) = -1.582916E-10
        coeff(4) = 1.525102E-11
        coeff(5) = 3.196347E-11
        coeff(6) = 211.8020
        coeff(7) = 139.8711
        coeff(8) = 217.9994

        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LH,I) = df
        So(LH,I)    = S_ent

!*******************************************************************
! FOR OH

! Temp = 298-1300 K
       if( (T(I) >=298.0E0) .and. (T(I) <=1300.0E0))then
          coeff(1) = 32.27768
          coeff(2) = -11.36291
          coeff(3) = 13.60545
          coeff(4) = -3.846486
          coeff(5) = -0.001335
          coeff(6) = 29.75113
          coeff(7) = 225.5783
          coeff(8) = 38.98706

! Temp = 1300-6000 K
       else if(T(I)>1300.0E0)then
          coeff(1) = 28.74701
          coeff(2) = 4.714489
          coeff(3) = -0.814725
          coeff(4) = 0.054748
          coeff(5) = -2.747829
          coeff(6) = 26.41439
          coeff(7) = 214.1166
          coeff(8) = 38.98706
      endif

      call shomate(coeff,df,S_ent,T,i,NZ)
      HH298(LOH,I) = df
      So(LOH,I)    = S_ent
!*******************************************************************
! FOR CO2

! Temp = 298-1200 K
        if( (T(I) >=298.0E0) .and. (T(I) <=1200.0E0))then
            coeff(1) = 24.99735
            coeff(2) = 55.18696
            coeff(3) = -33.69137
            coeff(4) =  7.948387
            coeff(5) =  -0.136638
            coeff(6) = -403.6075
            coeff(7) = 228.2431
            coeff(8) = -393.5224

! Temp = 1200-6000 K
        else if(T(I)>1200.0E0)then
            coeff(1) = 58.16639
            coeff(2) = 2.720074
            coeff(3) = -0.492289
            coeff(4) =  0.038844
            coeff(5) =  -6.447293
            coeff(6) =  -425.9186
            coeff(7) = 263.6125
            coeff(8) = -393.5224
        endif

        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LCO2,I) = df
        So(LCO2,I)    = S_ent
!*******************************************************************
! FOR CO

! Temp = 298-1300 K
       if( (T(I) >=298.0E0) .and. (T(I) <=1300.0E0))then
            coeff(1) =  25.56759
            coeff(2) =  6.096130
            coeff(3) =  4.054656
            coeff(4) =  -2.671301
            coeff(5) =  0.131021
            coeff(6) =  -118.0089
            coeff(7) =  227.3665
            coeff(8) =  -110.5271

!Temp = 1300-6000 K
        else if(T(I)>1300.0E0)then
            coeff(1) = 35.15070  
            coeff(2) = 1.300095
            coeff(3) =  -0.205921
            coeff(4) = 0.013550
            coeff(5) = -3.282780
            coeff(6) = -127.8375 
            coeff(7) = 231.7120
            coeff(8) = -110.5271
        endif

        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LCO,I) = df
        So(LCO,I)    = S_ent
!*******************************************************************
! FOR HCO
! Temp = 298-1200
        if( (T(I) >=298.0E0) .and. (T(I) <=1200.0E0))then
             coeff(1) = 21.13803
             coeff(2) = 40.43610
             coeff(3) = -14.71337
             coeff(4) = 0.969010
             coeff(5) = 0.239639
             coeff(6) = 36.34712
             coeff(7) = 240.1695
             coeff(8) = 43.51402

!Temp = 1200-6000
        else if(T(I)>1200.0E0)then
             coeff(1) = 52.79371
             coeff(2) = 2.666155
             coeff(3) = -0.392339
             coeff(4) = 0.023808
             coeff(5) = -7.457018
             coeff(6) = 11.37797
             coeff(7) = 267.2798
             coeff(8) = 43.51402
        endif
        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LHCO,I) = df
        So(LHCO,I)    = S_ent
!*******************************************************************
! FOR H2CO

!Temp = 298-1200 K
        if( (T(I) >=298.0E0) .and. (T(I) <=1200.0E0))then
            coeff(1) = 5.193767
            coeff(2) = 93.23249
            coeff(3) = -44.85457
            coeff(4) = 7.882279
            coeff(5) = 0.551175
            coeff(6) = -119.3591
            coeff(7) = 202.4663
            coeff(8) = -115.8972

!Temp = 1200 - 6000 K
         else if(T(I)>1200.0E0)then
            coeff(1) = 71.35268
            coeff(2) = 6.174497
            coeff(3) = -1.191090
            coeff(4) = 0.079564
            coeff(5) = -15.58917
            coeff(6) = -170.6327
            coeff(7) = 262.3180
            coeff(8) = -115.8972
        endif

        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LH2CO,I) = df
        So(LH2CO,I)    = S_ent
!*******************************************************************
! FOR CH4

! Temp = 298 - 1300 K
        if( (T(I) >=298.0E0) .and. (T(I) <=1300.0E0))then
             coeff(1) = -0.703029
             coeff(2) = 108.4773
             coeff(3) = -42.52157
             coeff(4) = 5.862788
             coeff(5) = 0.678565
             coeff(6) = -76.84376
             coeff(7) = 158.7163
             coeff(8) = -74.87310

! Temp = 1300 - 6000 K
        else if(T(I)>1300.0E0)then
             coeff(1) = 85.81217
             coeff(2) = 11.26467
             coeff(3) = -2.114146
             coeff(4) = 0.138190
             coeff(5) = -26.42221
             coeff(6) = -153.5327
             coeff(7) = 224.4143
             coeff(8) = -74.87310
        endif

        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LCH4,I) = df
        So(LCH4,I)    = S_ent
!*******************************************************************
! FOR CH3
        coeff(1) =  -8.107e-25
        coeff(2) =   2.155e-20
        coeff(3) =  -2.344e-16
        coeff(4) =    1.312e-12
        coeff(5) =   -3.688e-09 
        coeff(6) =    2.366e-06
        coeff(7) =  0.01655
        coeff(8) = 28.17 

        HH298(LCH3,I) = coeff(1)*T(I)**8+coeff(2)*T(I)**7 +
     .                 coeff(3)*T(I)**6 + coeff(4)*T(I)**5+
     .                 coeff(5)*T(I)**4+coeff(6)*T(I)**3+
     .                 coeff(7)*T(I)**2 + coeff(8)*T(I) -9913.0
        
        So(LCH3,I)    = -8.870e-27*T(I)**8+ 2.396e-22*T(I)**7 
     .                 -2.713e-18*T(I)**6 + 1.674e-14*T(I)**5
     .                 -6.143e-11*T(I)**4+1.382e-07*T(I)**3
     .                 -1.958e-04*T(I)**2 + 0.2114*T(I) + 144.7
        HH298(LCH3,I) = HH298(LCH3,I)/1000.0E0
!*******************************************************************
! FOR CH3OH
        coeff(1) =  -5.828e-26
        coeff(2) =  -8.813e-22
        coeff(3) =   5.288e-17
        coeff(4) =   -7.611e-13
        coeff(5) =    5.612e-09 
        coeff(6) =   -2.474e-05
        coeff(7) =  0.0688
        coeff(8) =  7.159

        HH298(LCH3OH,I) = coeff(1)*T(I)**8+coeff(2)*T(I)**7 +
     .                 coeff(3)*T(I)**6 + coeff(4)*T(I)**5+
     .                 coeff(5)*T(I)**4+coeff(6)*T(I)**3+
     .                 coeff(7)*T(I)**2 + coeff(8)*T(I) -7615.0

        So(LCH3OH,I)    = -1.678e-27*T(I)**8+ 4.593e-23*T(I)**7 
     .                 -5.317e-19*T(I)**6 + 3.421e-15*T(I)**5
     .                 -1.367e-11*T(I)**4+3.725e-08*T(I)**3
     .                 -8.041e-05*T(I)**2 + 0.1788*T(I) + 192.8

        HH298(LCH3OH,I) = HH298(LCH3OH,I)/1000.0E0

!*******************************************************************
! FOR CH

! Temp = 298 - 1100 K
       if( (T(I) >=298.0E0) .and. (T(I) <=1100.0E0))then
            coeff(1) = 32.94210
            coeff(2) = -16.71056
            coeff(3) = 24.18595
            coeff(4) = -7.784709
            coeff(5) = -0.065198
            coeff(6) = 584.6303
            coeff(7) = 226.5138
            coeff(8) = 594.1280

! Temp = 1100 - 6000 K
       else if(T(I)>1100.0E0)then
            coeff(1) = 30.15367
            coeff(2) = 8.455112
            coeff(3) = -1.969644
            coeff(4) = 0.154270
            coeff(5) = -4.980090
            coeff(6) = 576.6891
            coeff(7) = 209.3531
            coeff(8) = 594.1280
        endif

      call shomate(coeff,df,S_ent,T,i,NZ)
      HH298(LCH,I) = df
      So(LCH,I)    = S_ent
!*******************************************************************
! FOR CH23

! Temp = 298 - 1400 K
      if( (T(I) >=298.0E0) .and. (T(I) <=1400.0E0))then
             coeff(1) = 31.96823
             coeff(2) = 6.783603
             coeff(3) = 12.51890
             coeff(4) = -5.696265
             coeff(5) = -0.031115
             coeff(6) = 376.3558
             coeff(7) = 229.9150
             coeff(8) = 386.3924

! Temp = 1400 - 6000 K
       else if(T(I)>1400.0E0)then
             coeff(1) = 51.55901
             coeff(2) = 3.876975
             coeff(3) = -0.649608
             coeff(4) = 0.037901
             coeff(5) = -10.72589
             coeff(6) = 350.6715
             coeff(7) = 232.3212
             coeff(8) = 386.3924
 
        endif

        call shomate(coeff,df,S_ent,T,i,NZ)
         HH298(LCH23,I) = df
         So(LCH23,I)    = S_ent

!*******************************************************************
! FOR H2

!Temp = 298-1000 K
        if( (T(I) >=298.0E0) .and. (T(I) <=1000.0E0))then
             coeff(1) = 33.066178
             coeff(2) = -11.363417
             coeff(3) = 11.432816
             coeff(4) = -2.772874
             coeff(5) = -0.158558
             coeff(6) = -9.980797
             coeff(7) = 172.707974
             coeff(8) = 0.0


!Temp = 1000-2500 K
        else if( (T(I)>1000.0E0) .and. (T(I) <=2500.0E0))then
             coeff(1) = 18.563083
             coeff(2) = 12.257357
             coeff(3) = -2.859786
             coeff(4) = 0.268238
             coeff(5) = 1.977990
             coeff(6) = -1.147438
             coeff(7) = 156.288133
             coeff(8) = 0.0


!Temp = 2500-6000 K
         else
              coeff(1) = 43.413560
              coeff(2) = -4.293079
              coeff(3) = 1.272428
              coeff(4) = -0.096876
              coeff(5) = -20.533862
              coeff(6) = -38.515158
              coeff(7) = 162.081354
              coeff(8) = 0.0
         endif
          
        call shomate(coeff,df,S_ent,T,i,NZ)
        HH298(LH2,I) = df
        So(LH2,I)    = S_ent

!*******************************************************************
! Nothing for O1D
         coeff(1) = 0.0E0
         coeff(2) = 0.0E0
         coeff(3) = 0.0E0
         coeff(4) = 0.0E0
         coeff(5) = 0.0E0
         coeff(6) = 0.0E0
         coeff(7) = 0.0E0
         coeff(8) = 0.0E0
         call shomate(coeff,df,S_ent,T,i,NZ)
         HH298(LO1D,I) = df
         So(LO1D,I)    = S_ent

!*******************************************************************
! Nothing for CH21
         coeff(1) = 0.0E0
         coeff(2) = 0.0E0
         coeff(3) = 0.0E0
         coeff(4) = 0.0E0
         coeff(5) = 0.0E0
         coeff(6) = 0.0E0
         coeff(7) = 0.0E0
         coeff(8) = 0.0E0
!         call shomate(coeff,df,S_ent,T,i,NZ)
!         HH298(16,I) = df
!         So(16,I)    = S_ent
!*********************************************
! Difference between CH23 & CH21 is 31.1     

         HH298(LCH21,I) = HH298(LCH21,I) + 31.1
         So(LCH21,I)    = So(LCH21,I)

!*******************************************************************
! FOR H2COH

          coeff(1) =  6.03e-25 
          coeff(2) =  -1.735e-20
          coeff(3) =   2.154e-16 
          coeff(4) =  -1.529e-12
          coeff(5) =   6.969e-09
          coeff(6) =  -2.183e-05
          coeff(7) =  0.04864
          coeff(8) =  23.95

          HH298(LH2COH,I) = coeff(1)*T(I)**8+coeff(2)*T(I)**7 +
     .                   coeff(3)*T(I)**6 + coeff(4)*T(I)**5+
     .                   coeff(5)*T(I)**4+coeff(6)*T(I)**3+
     .                   coeff(7)*T(I)**2 + coeff(8)*T(I) -10950.0

          So(LH2COH,I) =    -4.169e-27*T(I)**8+ 1.166e-22*T(I)**7 
     .                   -1.382e-18*T(I)**6 + 9.067e-15*T(I)**5
     .                   -3.622e-11*T(I)**4+9.241e-08*T(I)**3
     .                   -0.0001592*T(I)**2 + 0.2251*T(I) + 189.1

           HH298(LH2COH,I) = HH298(LH2COH,I)/1000.0E0

!*******************************************************************
! FOR C
          coeff(1) = 21.17510
          coeff(2) = -0.812428
          coeff(3) = 0.448537
          coeff(4) = -0.043256
          coeff(5) = -0.013103
          coeff(6) = 710.3470
          coeff(7) = 183.8734
          coeff(8) = 716.6690
          call shomate(coeff,df,S_ent,T,i,NZ)
          HH298(LC,I) = df
          So(LC,I)    = S_ent
!*******************************************************************
! FOR CH3O

           coeff(1) =   6.948e-25
           coeff(2) =   -2.069e-20
           coeff(3) =   2.669e-16
           coeff(4) =  -1.966e-12
           coeff(5) =   9.205e-09
           coeff(6) =    -2.882e-05
           coeff(7) =   0.06176
           coeff(8) =   18.09

           HH298(LCH3O,I) = coeff(1)*T(I)**8+coeff(2)*T(I)**7 +
     .                    coeff(3)*T(I)**6 + coeff(4)*T(I)**5+
     .                    coeff(5)*T(I)**4+coeff(6)*T(I)**3+
     .                    coeff(7)*T(I)**2 + coeff(8)*T(I) -10210.0
     
           HH298(LCH3O,I) = HH298(LCH3O,I)/1000.0E0

           So(LCH3O,I) = -3.068e-27*T(I)**8+ 8.627e-23*T(I)**7 +
     .                 (-1.030e-18*T(I)**6) + 6.838e-15*T(I)**5+
     .                  (-2.791e-11*T(I)**4)+7.413e-08*T(I)**3 + 
     .                  (-0.0001377*T(I)**2) + 0.2195*T(I) + 181.2

      enddo    ! End of the loop

!*******************************************************************

      do I = 1,NZ

! Since the element is 'O2', Divide HH298 by 2
          Df_Ho(LO,I) = Df_Ho298(LO) + HH298(LO,I) - (HH298(LO2,I)
     .                  /2.0E0)
          Df_Go(LO,I) = Df_Ho(LO,I)*1000.0E0 - T(I)*(So(LO,I) - 
     .                 (So(LO2,I)/2.0E0))


! Enthalpy of formation for element 'O2' is zero
          Df_Ho(LO2,I) = 0.0E0  
          Df_Go(LO2,I) = 0.0E0


! For 'H2O', divide 'O2' by 2 and not element 'H2' by 2 because
! there are two hydrogen .
          Df_Ho(LH2O,I) = Df_Ho298(LH2O) + HH298(LH2O,I)
     .                 - ( HH298(LH2,I) + (HH298(LO2,I)/2.0E0) )
          Df_Go(LH2O,I) = Df_Ho(LH2O,I)*1000.0E0 - T(I)*(So(LH2O,I) - 
     .                 (So(LH2,I) + (So(LO2,I)/2.0E0)))


! For 'H', element is 'H2', so divide by 2
          Df_Ho(LH,I) = Df_Ho298(LH) + HH298(LH,I)-(HH298(LH2,I)/2.0E0)
          Df_Go(LH,I) = Df_Ho(LH,I)*1000.0E0 - T(I)*(So(LH,I) -
     .                 (So(LH2,I)/2.0E0))
               


! For 'OH', divide HH298 by 2 for both O2 & H2
          Df_Ho(LOH,I) = Df_Ho298(LOH) + HH298(LOH,I) -( 
     .                   (HH298(LO2,I)/2.0E0) + (HH298(LH2,I)/2.0E0))
          Df_Go(LOH,I) = Df_Ho(LOH,I)*1000.0E0 - T(I)*(So(LOH,I) - 
     .                 ((So(LO2,I)/2.0E0) + (So(LH2,I)/2.0E0)))
 

! For 'CO2', do not divide element 'O2' by 2 because 'CO2' has two
! atoms.
          Df_Ho(LCO2,I) = Df_Ho298(LCO2) + HH298(LCO2,I) - 
     .                 (HH298C(I) + HH298(LO2,I))
          Df_Go(LCO2,I) = Df_Ho(LCO2,I)*1000.0E0 - T(I)*(So(LCO2,I) - 
     .                 (SoC(I) + So(LO2,I)))
     
! For 'CO', divide element 'O2' by 2
          Df_Ho(LCO,I) = Df_Ho298(LCO) + HH298(LCO,I) - (HH298C(I) + 
     .                 (HH298(LO2,I)/2.0E0))
           Df_Go(LCO,I) = Df_Ho(LCO,I)*1000.0E0 - T(I)*(So(LCO,I) -
     .                  ((So(LO2,I)/2.0E0) + SoC(I)) )


! For 'HCO'  
           Df_Ho(LHCO,I) = Df_Ho298(LHCO) + HH298(LHCO,I) - 
     .                     ((HH298(LH2,I)/2.0E0)+
     .                  HH298C(I)   +  (HH298(LO2,I)/2.0E0))
           Df_Go(LHCO,I) = Df_Ho(LHCO,I)*1000.0E0 - T(I)*(So(LHCO,I) -
     .                  ( (So(LH2,I)/2.0E0)+SoC(I) + 
     .                  (So(LO2,I)/2.0E0) ) )



! For 'H2CO'
            Df_Ho(LH2CO,I) = Df_Ho298(LH2CO) + HH298(LH2CO,I) 
     .                       - (HH298(LH2,I) +
     .                   HH298C(I) + (HH298(LO2,I)/2.0E0))
            Df_Go(LH2CO,I) = Df_Ho(LH2CO,I)*1000.0E0 - T(I)
     .                       *(So(LH2CO,I) - 
     .                     (So(LH2,I) + SoC(I) + (So(LO2,I)/2.0E0)) )

 

! For 'CH4'
            Df_Ho(LCH4,I) = Df_Ho298(LCH4) + HH298(LCH4,I) 
     .                      -(HH298C(I) + 
     .                      2.00E0*HH298(LH2,I))
            Df_Go(LCH4,I) = Df_Ho(LCH4,I)*1000.0E0 - T(I)
     .                      *(So(LCH4,I) - 
     .                      (SoC(I) + 2.0E0*So(LH2,I)) )


! For 'CH3'
            Df_Ho(LCH3,I) = Df_Ho298(LCH3) + HH298(LCH3,I) 
     .                      -(HH298C(I) + 
     .                      1.50E0*HH298(LH2,I))
            Df_Go(LCH3,I) = Df_Ho(LCH3,I)*1000.0E0 - 
     .                      T(I)*(So(LCH3,I) - 
     .                    (SoC(I) + 1.50E0*So(LH2,I) ))


! For 'CH3OH'
             Df_Ho(LCH3OH,I) = Df_Ho298(LCH3OH) + HH298(LCH3OH,I)
     .                         -(HH298C(I) + 
     .                         2.00E0*HH298(LH2,I)+ 
     .                         (HH298(LO2,I)/2.0E0))
             Df_Go(LCH3OH,I) =     Df_Ho(LCH3OH,I)*1000.0E0 - 
     .                         T(I)*(So(LCH3OH,I) - 
     .                         (SoC(I)+ 2.0E0*So(LH2,I) 
     .                         + (So(LO2,I)/2.0E0)))


! For 'CH'
             Df_Ho(LCH,I) = Df_Ho298(LCH) + HH298(LCH,I) -(HH298C(I) + 
     .                     (HH298(LH2,I)/2.0E0))
             Df_Go(LCH,I) = Df_Ho(LCH,I)*1000.0E0 - T(I)*(So(LCH,I) - 
     .                     (SoC(I) + (So(LH2,I)/2.0E0)) )


! For 'CH23'
             Df_Ho(LCH23,I) = Df_Ho298(LCH23) + HH298(LCH23,I) 
     .                       -(HH298C(I) + HH298(LH2,I))
             Df_Go(LCH23,I) = Df_Ho(LCH23,I)*1000.0E0 - 
     .                       T(I)*(So(LCH23,I) - 
     .                     (SoC(I) + So(LH2,I)) )

! For 'H2'
             Df_Ho(LH2,I) = 0.0E0
             Df_Go(LH2,I) = 0.0E0

! For 'O1D'
             Df_Ho(LO1D,I) = 0.0E0
             Df_Go(LO1D,I) = 0.0E0
         


! For 'CH21'

             Df_Ho(LCH21,I) = Df_Ho298(LCH21) + HH298(LCH21,I)
     .                        -(HH298C(I) +  HH298(LH2,I))
             Df_Go(LCH21,I) = Df_Ho(LCH21,I)*1000.0E0 - 
     .                        T(I)*(So(LCH21,I) -
     .                     (SoC(I) + So(LH2,I)) )



! For 'H2COH'
             Df_Ho(LH2COH,I) = Df_Ho298(LH2COH) + HH298(LH2COH,I)
     .                         - (HH298C(I) +
     .                        1.50E0*HH298(LH2,I)+ (HH298(LO2,I)/2.0E0))
             Df_Go(LH2COH,I) = Df_Ho(LH2COH,I)*1000.0E0 - 
     .                         T(I)*(So(LH2COH,I) -
     .                         (SoC(I) + 1.50E0*So(LH2,I) + 
     .                         (So(LO2,I)/2.0E0)))
     
! For 'C'
             Df_Ho(LC,I) = Df_Ho298(LC) + HH298(LC,I) - HH298C(I)
             Df_Go(LC,I) = Df_Ho(LC,I)*1000.0E0 - T(I)*(So(LC,I) - 
     .                     (SoC(I)) ) 

! For 'CH3O'
             Df_Ho(LCH3O,I) = Df_Ho298(LCH3O) + HH298(LCH3O,I)
     .                        - (HH298C(I) + 
     .                     1.50E0*HH298(LH2,I)+ (HH298(LO2,I)/2.0E0))
             Df_Go(LCH3O,I) = Df_Ho(LCH3O,I)*1000.0E0 - 
     .                        T(I)*(So(LCH3O,I) - 
     .                        (SoC(I) + 1.50E0*So(LH2,I) 
     .                        + (So(LO2,I)/2.0E0)))
      enddo

*****************************************************************
c-mc rate constant units are cm^3/mol/s

       ! chemical reaction file
       open(9, file='PHOTOCHEM/INPUTFILES/reactions.rx',status='OLD')
 
  	print*,'Running rateswasp12b.f...'
  	
 !old (Atmos) rates.f has a different format for 2body and 3body cases:
 !old 2body: TWO values read: A0 in E format, e0 in F format
 !old 3body: FOUR values read: A0 and A1 in E format, n0 and n1 in F (for the 2 pressure limits)
 !667   FORMAT(58X,E9.2,3X,F8.2)            !for two body reaction rates
 !668   FORMAT(58X,E9.2,3X,E9.2,2X,2F5.2)   !for three body reaction rates
 
 !the hot jupiter reaction.rx will have the same format for both...
 !new 2body/3body: SIX values read: A0, A1 in E, n0, n1 in F, then e0 and e1 in E (activation energy):
 !for 2body cases, A1,n1 and E1 are padded in as 0's (not used in computation).
  !In other words, in the the existing reaction.rx, this needs to be updated (by 0 padding).
 667   FORMAT(58X,E9.3,3X,E9.3,2X,2F5.2,1X,E10.3,1X,E10.3)            !for two body reaction rates

        kB = 1.38054E-16 !Boltzmann constant in cgs (erg K-1)
        PFAC = 1.013E+06 !Standard pressure in cgs (dyne/cm^2)
        print*,'kB,PFAC =',kB,PFAC
                     
       do J=1,NR
C READ IN REACTION INFORMATION
           
           !WHAT THE TERMS MEAN, in rate eq. A = a*T^tn*exp(-e/T)
           ! a0, a1 = rate coefficient high P and low P limits, respectively
           ! tn0, tn1 = temperature coefficient limits
           ! e0, e1 = activation energy limits
           ! Only one set of a, tn, and e's for 2body: 0's padded for rest
C COMPUTE TWO BODY REACTION RATES
          if (REACTYPE(J) .EQ. '2BODY') then
           read (9,667) a0,a1,tn0,tn1,e0,e1         
           !PRINT*, "J,a0,a1,tn0,tn1,e0,e1 = ",J,a0,a1,tn0,tn1,e0,e1
           do I=1,NZ
            A(J,I)=a0*(T(I)/298.)**tn0*EXP(-e0/T(I))   !two body reaction rates
!!            A(5,I)=0.350E-12*(T(I)/298.0E0)**( 0.267E+01)
!!     .       *exp(-0.316E+04/T(I))
            !IF(I.EQ.1)print*,'J, A(NZ=1) (2body rxn)',J,A(J,I)
            IF(J.EQ.151)A(J,I)=0.0
C            IF(I.EQ.50)print*,
C     .      'J, A(highest P), A(Z=50) (2body rxn)',J,A(J,1),A(J,50)
           enddo

C COMPUTE THREE BODY REACTION RATES
          else if (REACTYPE(J) .EQ. '3BODY') then
           read (9,667) a0,a1,tn0,tn1,e0,e1
           !PRINT*, "J,a0,a1,tn0,tn1,e0,e1 = ",J,a0,a1,tn0,tn1,e0,e1
             if (PLANET .EQ. 'MARS') then
                B=B*2.5      !multiply low density rate by 2.5 to account for CO2 rather than N2 as background gas (Nair, 94)
             endif
           do I=1,NZ
              A(J,I)= TBDYR(a0,a1,tn0,tn1,e0,e1,T(I),DEN(I))  !computed three body reaction rate
              !IF(I.EQ.1)print*,'J, A(NZ=1) (3body rxn)',J,A(J,I)
C              IF(I.EQ.50)print*,
C     .        'J, A(highest P), A(Z=50) (3body rxn)',J,A(J,1),A(J,50)
           enddo

C COMPUTE BACKWARD REACTION RATES (J) BASED ON PREVIOUS FORWARD REACTION (J-1) (see below).
C This code ASSUMES that the preceding reaction to J, i.e. J-1, corresponds to reaction J's forward version, thus 
C a0, a1, etc. obtained from (J-1) reaction are being substituted directly. 
C If user cannot uncertain this for their specific reactions.rx file, uncomment the print statement.
        
          else if (REACTYPE(J) .EQ. '2BACK') then !BACKWARD REACTIONS TO 2 BODY
           read (9,*)        
C           PRINT*, "J, Reactants (Back)= ",J,CHEMJ(1,J),CHEMJ(2,J)
C           PRINT*,"Products (Back)= ",CHEMJ(3,J),CHEMJ(4,J),CHEMJ(5,J)   
           DO I=1,NZ
             DG    = Df_Go(JCHEM(3,J),I)+
     .        Df_Go(JCHEM(4,J),I)-(Df_Go(JCHEM(1,J),I)+
     .        Df_Go(JCHEM(2,J),I))
              IF(J.EQ.86)DG=-DG !reversed for some reason in Ravi's
             akeq  = exp(DG/(RR*T(I))) 
             A(J,I)=A(J-1,I)/akeq
            IF(J.EQ.152)A(J,I)=0.0
C             IF(I.EQ.40) print*,
C     .        'J, I, A(Z=40) (2 back),DG,RR,T(I)',
C     .        J,I,A(J,I),DG,RR,T(I)
           ENDDO
C           PRINT*," "   
           
           !Here is an hardcoded example from Ravi's rates.f file for 2body backward:            
             ! DG    = Df_Go(LH2,I)+Df_Go(LO,I)-(Df_Go(LOH,I)+Df_Go(LH,I))
            !A(6,I)=A(5,I)/exp(DG/(RR*T(I)))
                 
          else if (REACTYPE(J) .EQ. '3BACK') then !BACKWARD REACTIONS TO 3 BODY
           read (9,*)        
C           PRINT*, "J, Reactants (Back)= ",J,CHEMJ(1,J),CHEMJ(2,J)
C           PRINT*,"Products (Back)= ",CHEMJ(3,J),CHEMJ(4,J),CHEMJ(5,J)          
           DO I=1,NZ 
           !Using general equation - setting Df_Go = 0.0 for blank species
            IF(JCHEM(2,J).EQ.0)Df_Go(JCHEM(2,J),I)=0.0
            IF(JCHEM(4,J).EQ.0)Df_Go(JCHEM(4,J),I)=0.0              
             DG    = Df_Go(JCHEM(3,J),I)+
     .        Df_Go(JCHEM(4,J),I)-(Df_Go(JCHEM(1,J),I)+
     .        Df_Go(JCHEM(2,J),I))
             akeq  = exp(DG/(RR*T(I)))
             akr0  = (a0/akeq)*(PFAC/(kB*T(I)))
             akri  = (a1/akeq)*(PFAC/(kB*T(I)))
             ! akri  = (a1/akeq)!*(PFAC/(kB*T(I)))
 !            IF(I.EQ.50) print*,
 !    .    'J, A(highestP), A(Z=50) (3body back rxn to above W/ PFAC)',
 !    .     J,A(J,1),A(J,50)
 
      !!!!!!!!! DEBUG !!!!!
C            IF (J.EQ.118.AND.I.EQ.40) THEN
C             print*,'DG,RR,T(I),akeq,akr0,akri =',
C     .        DG,RR,T(I),akeq,akr0,akri
C             print*,'J,I,CHEMJ(1),JCHEM(1),Df_Go(JCHEM(1,J),I):',
C     .        J,I,CHEMJ(1,J),JCHEM(1,J),Df_Go(JCHEM(1,J),I)
C             print*,'J,I,CHEMJ(2),JCHEM(2),Df_Go(JCHEM(2,J),I):',
C     .        J,I,CHEMJ(2,J),JCHEM(2,J),Df_Go(JCHEM(2,J),I)
C             print*,'J,I,CHEMJ(3),JCHEM(3),Df_Go(JCHEM(3,J),I):',
C     .        J,I,CHEMJ(3,J),JCHEM(3,J),Df_Go(JCHEM(3,J),I)
C             print*,'J,I,CHEMJ(4),JCHEM(4),Df_Go(JCHEM(4,J),I):',
C     .        J,I,CHEMJ(4,J),JCHEM(4,J),Df_Go(JCHEM(4,J),I)
C            ENDIF
     !!!!!!!!
C           IF (I.EQ.1) print*,'J,T,D,akeq,a0,akr0,a1,akri,DG,A(s)',
C     .        J,,akeq,a0,akr0,a1,akri,DG,A(J,I)
           !IF (I.EQ.1) print*,'J,akeq,DG',J,akeq,DG               
             !print*,'J,A (hardcoded)',J,A(J,I)
             IF(JCHEM(4,J).EQ.0) THEN
c-mab: I have a reservation about this. I believe original scheme is missing a DEN factor consideration....
c-mab: Perhaps we should switch order of the forward and backward?
              A(J,I)=(A(J-1,I)/akeq)*(DEN(I))**(-1) !IMPORTANT! Read the comments below.
             ELSE
              A(J,I)=A(J-1,I)*(DEN(I))/akeq !IMPORTANT! Read the comments below.

             !A(J,I)=A(J-1,I)/akeq !IMPORTANT! Read the comments below.
             ENDIF
                         
C             IF(I.EQ.50) print*,
C     .    'J, A(highestP) (3body back rxn to above w/ PFAC)',J,A(J,50)
C             IF(I.EQ.50) print*,'ARAVI = ',ARAVI


           !Attempts to simplify Ravi's expressions for the back reactions
          !B0 = (A0*(300./T)**CN*exp(-E0/T)/akeq)*(PFAC/(kB*T(I)))
          !BI = (A1*(300./T)**CM*exp(-EI/T)/akeq)*(PFAC/(kB*T(I)))
          !A(J-1,I) = B0*D/(1. + B0*D/BI)
          
           !Here is an hardcoded example from Ravi's rates.f file for 2body backward:            
           !DG    = Df_Go(LH,I)+Df_Go(LOH,I)-(Df_Go(LH2O,I))
           !akeq  = exp(DG/(RR*T(I)))
           !akr0  = (0.660E-31/akeq)*(1.013E+06/(BK*T(I)))
           !akri  = (0.270E-09/akeq)*(1.013E+06/(BK*T(I)))
           !A(8,I)=TBDYR(akr0,akri,-2.10, 0.00, 0.000E+00,0.750E+02,TT,DN)
           ENDDO
C           PRINT*," " 
     
          else if (REACTYPE(J) .EQ. 'PHOTO') then
           read (9,*)
C           PRINT*, "J, Reactants = ",J,CHEMJ(1,J),CHEMJ(2,J)
C           PRINT*,"J, Products = ",J,CHEMJ(3,J),CHEMJ(4,J)                    
          	 print*, "Nothing. Photolysis rxs not handled in rates.f." 
          else
             print*,"Unfamiliar reaction type"
             print*,"Aborting..."
             stop
          endif
       enddo

!          ENDDO
!         else
!             read(9,*) !do nothing with PHOTO...
!          endif
!       enddo

       close(9)

      RETURN
      END

       !The stuff below has been ported from Ravi's WASP12b work rates.f file 
        subroutine shomate(coeff,df,S_ent,T,i,NZ)
       real *8 df,S_ent,coeff(8),T(NZ)

       T(i) = T(i)/1000.0E0

       df = coeff(1)*T(i) + coeff(2)*(T(i)**2/2.0) + 
     .      coeff(3)*(T(i)**3/3.0E0) +
     .      coeff(4)*(T(i)**4/4.0E0) - (coeff(5)/T(i)) 
     .      + coeff(6) - coeff(8)

       S_ent = coeff(1)*log(T(i)) + coeff(2)*T(i) +
     .         coeff(3)*(T(I)**2/2.0E0) +
     .         coeff(4)*(T(I)**3/3.0E0) - 
     .         (coeff(5)/(2.0E0*T(I)**2)) + coeff(7)

        T(i) = T(i)*1000.0E0
        end
                     
        FUNCTION TBDYR(A0,AI,CN,CM,E0,EI,T,D)
        real*8 tbdyr
          B0 = A0*(300./T)**CN*exp(-E0/T)
          BI = AI*(300./T)**CM*exp(-EI/T)
!          print *,'a0 = ',a0,'ai = ',ai
!          print *,'cn = ',cn,'cm = ',cm
!          print *,'e0 = ',e0,'ei = ',ei
!          print *,TT,D
!          print *,B0,BI
!          B0 = A0*(300./T)**CN
!          BI = AI*(300./T)**CM
!           Y = ALOG10(B0*D/BI)
!           X = 1./(1. + Y**2)
!           TBDYR = B0*D/(1. + B0*D/BI) * 0.6**X
           TBDYR = B0*D/(1. + B0*D/BI) 
        RETURN
      
      END
