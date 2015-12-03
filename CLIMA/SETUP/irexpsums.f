
c     SUBROUTINE IREXPSUMS(weight,xkappa)
      SUBROUTINE IREXPSUMS

c-mm  This subroutine will READ in the exponential sum datafiles in the
c-mm  ir for both h2o and co2 and output the appropriate value for each
c-mm  wavelength interval.  This will crash if our values exceed the
c-mm  boundaries set forth in the data, but this shouldn't happen under
c-mm  normal Martian conditions.
c-mm  Indices for xkappa and kappa are as follows
c-mm  Index #1:  Temp   1=less than T   2=more than T
c-mm  Index #2:  Pres   1=less than p   2=more than p
c-mm  Index #3:  lam    wavelength bin
c-mm  Index #4:  i      gauss point
c-mm  Index #5:  Spec   1=co2  2=h2o  3=ch4

      include 'CLIMA/INCLUDE/header.inc'
c      PARAMETER (NF=55,NGS=7)
      PARAMETER(NF = 55)
      INTEGER i,lam,junk,pind,j,k, Tind2  ! Added Tind2 and removed Tinddum 3/23/2012
                                                        ! EWS - removed unused l, LINE, m, and Tind 8/26/2015
      REAL  CIA1(7,NF), kappa_irh2o, kappa_irco2 ! Declaring CO2 and H2O coefficient arrays as real. 8/27/2012
                                                   ! EWS - rmeoved unused variable zeroed 8/26/2015 
      PARAMETER (IK=8)
c     REAL xkappa(8,12,55,8,3)
c     DIMENSION WEIGHT(8,3)
           CHARACTER :: DIRINOUT*8,DIRDATA*10
      COMMON/DIR/DIRINOUT,DIRDATA
C new common block, von Paris, 21/04/2006
      COMMON/IRDATA/WEIGHTCH4(6),xkappa(3,12,NF,IK),
     & CIA(7,NF), CPRW(ND,NF)!c-rr !3/23/11 put CIA matrix in IRDATA. Added weight co2_h2O 3/20/2012

      COMMON/VARIR/kappa_irh2o(NF, 8, 8, IK), kappa_irco2(NF, 8, 8, IK)! Added kappa matrix in IR for kpsectrum Co2 and H2O coefficients 8/26/2012  

      COMMON/HYDN2CIA/ H2N2CIA(5,NF), H2N2FIN(ND,NF) ! c-rr 5/23/12 H2-N2 CIA coefficients from Borysow et al., (2002).

      COMMON/HYDCIA/ H2H2CIA(6,NF), H2H2FIN(ND,NF) ! c-rr 7/02/2012 H2-H2 CIA coefficients from Borysow et al. 

      COMMON/OXYCIA/O2O2CIA(15,NF), O2O2FIN(ND,NF) ! c-rr 6/17/2012 O2-O2 CIA coefficients from HITRAN CIA database

      COMMON/BPS_IR/s_abir(NF), f_abir(NF), TDir(NF), Bsir(NF), 
     &    Bfir(NF)  ! Added COMMON BLOCK FOR BPS CONTINUUM FOR IR  8/30/2012
      COMMON/weightsIR/weightco2_h2oIR(IK)
      
c-rr  The 8 weights below are used for both the CO2 and H2O kspectrum coefficients 8/26/2012
      DATA weightco2_h2oIR/1.65231051440290516E-01, 
     &  3.09768948559709378E-01,  
     &  3.09768948559709378E-01, 
     &  1.65231051440290516E-01, 
     &  8.69637112843634277E-03, 
     &  1.63036288715636482E-02, 
     &  1.63036288715636482E-02, 
     &  8.69637112843634277E-03 /


      OPEN(unit=20,file= DIRDATA//'/ir_expsums.pdat')


c-rr  CH4 below this line--------------------------------
      OPEN(unit=23, file=DIRDATA//'/ch4_kcoeffs.out')

      WEIGHTCH4(1)=0.08566225
      WEIGHTCH4(2)=0.18038079
      WEIGHTCH4(3)=0.23395697
      WEIGHTCH4(4)=0.23395697
      WEIGHTCH4(5)=0.18038079
      WEIGHTCH4(6)=0.08566225

      do ii=1,3
       READ(23,*) ! skips first 3 lines 3/23/2012
      enddo
      do 9600 Tind2=1,3 !Removed Tinddum and using Tind2 to represent indices for T=150, 300, and 600K (1,2, and 3 respectively instead of 3,6, and 8).3/23/2012
         do 9601 pind=7,12
            do 9602 lam=1,NF  ! Changing intervals from 1 to 55. CH4 only has values though for lam = 5 to 42. 3/23/2012
               READ(23,*) junk,(xkappa(Tind2,pind,lam,ii), ii=1,6)
!               print *, (xkappa(Tind2,pind,lam,ii,3), ii=1,6)
!               pause
 9602          continue
              READ(23,*)
              READ(23,*)
              READ(23,*) ! skips three lines
 9601          continue
 9600          continue

        close(23)

! Extraneous stuff in code? 3/23/2012
!        do 9570 i=7,12
c        do 9571 j=1,38
!        do 9571 j=5,42
!        do 9572 k=1,6
!          xkappa(1,i,j,k,3)=xkappa(1,i,j,k,3)
c          xkappa(1,i,j,k,3)=xkappa(3,i,j,k,3)
! 9572 continue
! 9571 continue
! 9570 continue





C-rr READ in the new separate CO2 and H2O absorption coefficients for IR intervals 1-55 8/26/2012

C-----------------------------------------------------------------------

! Initializing k-coefficient array
        do i = 1,NF
                 do it = 1, 8 ! 8 temperatures
                            do ip = 1,8  ! 8 pressures
                                       do k = 1,IK
                                    kappa_irh2o(i,it,ip,k) = 1.e-60
                                    kappa_irco2(i,it,ip,k) = 1.e-60
                                       enddo ! k loop
                           enddo ! ends pressure loop
                enddo ! ends temperature loop
        enddo  ! ends interval loop




      do iii = 1,16
              read(17,*) ! Initially skip 16 lines
              read(18,*) !
      enddo


        do i = 1,NF

                 do it = 1, 8 ! 8 temperatures
                        read(17,*) ! skip 3 more lines to read data
                        read(17,*)
                        read(17,*)
                        read(18,*) ! skip 3 more lines to read data
                        read(18,*)
                        read(18,*)
                            do ip = 1,8  ! 8 pressures
                            read(17,*)a,(kappa_irh2o(i,it,ip,k),k=1,IK)
                            read(18,*)a,(kappa_irco2(i,it,ip,k),k=1,IK)


!                           IF ((i.ge.1).and.(i.le.8))then
!                                 do k = 1,IK
!                                     kappa_irco2(i,it,ip,k) = 1.d-60 ! No CO2 coefficients in these intervals according to irexpsums.pdat 8/29/2012
!                                 enddo
!                              ENDIF
          
!                              IF((i.ge.49).and.(i.le.NF))then
!                                do k = 1,IK  
!                                    kappa_irco2(i,it,ip,k) = 1.d-60 ! No CO2 coefficients in these intervals according to irexpsums.pdat 8/29/2012
!                                      enddo
!                           ENDIF

                          enddo ! ends pressure loop
                enddo ! ends temperature loop
        enddo  ! ends interval loop


C-rr-----------------------------------------

C-rr reads CO2 CIA matrix(3/22/2011)----------------------
        read(30,9800) !skips two lines
        do k=1,9
                read(30,9801) (CIA1(it,k),it=1,7)  ! Put CIA values in the right bins 
        enddo
        

        do k=20,27
                read(30,9801) (CIA1(it,k),it=1,7)  ! Put CIA values in the right bins
        enddo


        do it=1,7
                do k=1,NF
                     CIA(it,k)=amax1(CIA1(it,k),1.e-60) ! set all zero CIA values to 1.e-60
                enddo
        enddo


 9800   format(/)
 9801   format(19x, 1pe24.17,6(3x,1pe24.17))
c 9802   format(i3, 3x,1pe24.17,6(3x,1pe24.17)) !EWS - not used


c 9950   format(/)         !EWS - not used
c 9952   format(27X,E16.10)!EWS - not used
c 9956   format(83X,E9.3)  !EWS - not used
c 9958   format(56X,E12.6) !EWS - not used
c 9960   format(21X,E12.6) !EWS - not used


c-rr -------------------------------------------------------

c-rr   reads H2-N2 CIA matrix 4/29/2012

       OPEN(UNIT = 25, file=DIRDATA//'/H2N2CIA_TABLE.pdat')



        do i = 1, 3
         read(25,*) ! skips three lines before reading H2 CIA
        enddo

       do it =1,5
          do k = 1,NF
          H2N2CIA(it,k) = 1.e-60 ! initializing H2 CIA matrix first.
          enddo
       enddo

       
              do k = 1,29
               read(25, *)a, (H2N2CIA(it,k), it=1,5)
!               print *, (H2N2CIA(it,k),it=1,6)
!               pause
              enddo  ! ends k loop
       
c-rr-----------------------------------------------------

c-rr   read O2-O2 CIA matrix 6/17/2012

       OPEN(UNIT = 26, file=DIRDATA//'/O2O2CIA_TABLE.pdat')

        do i = 1, 3
         read(26,*) ! skips three lines before reading O2 CIA
        enddo

       do it =1,15
          do k = 1,NF
          O2O2CIA(it,k) = 1.e-60 ! initializing O2 CIA matrix first.
          enddo
       enddo


              do k = 21,28
               read(26, *)a, (O2O2CIA(it,k), it=1,15)
!               print *, (O2O2CIA(it,k),it=1,6)
!               pause
              enddo  ! ends k loop         


c-rr ----------------------------------------------------

c-rr    reads in H2-H2 CIA matrix 7/2/2012

           OPEN(unit = 27, file = DIRDATA//'/H2CIA_TABLE.pdat')

        do i = 1,3
         read(27,*) ! skips three lines before reading H2 CIA
        enddo

        do it = 1,6
            do k = 1,NF
            H2H2CIA(it,k) = 1.e-60 ! initializing H2 CIA matrix first

            enddo
        enddo


        do k = 1,NF
         read(27,*)a, (H2H2CIA(it,k), it =1,6)
!         print *, (H2H2CIA(it,k),it=1,6)
!         pause
        enddo
!-----------------------------------------------------------------


c-rr   reads in BPS CONTINUUM coefficients at T = 296K for IR  8/30/2012

            OPEN(unit = 36, file = DIRDATA//'/IR_BPS.dat')

  ! initialize arrays

         do i = 1,NF
                s_abir(i) = 0.0d0
                f_abir(i) = 0.0d0
                TDir(i) =0.0d0
                Bsir(i) = 0.0d0
                Bfir(i) =0.0d0
         enddo
 
      
        do i = 1,5
                read(36,*)
        enddo


        do j = 1,NF
        read(36,*)a,b,s_abir(j), f_abir(j), Bsir(j), Bfir(j), TDir(j)
!                print *, s_abir(j), f_abir(j), Bsir(j), Bfir(j), TDir(j)
        enddo

!--------------------------------------------------------------------         
        RETURN
        END
