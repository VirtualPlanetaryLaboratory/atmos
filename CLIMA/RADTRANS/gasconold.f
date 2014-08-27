      SUBROUTINE GASCON(T,PF,FO2,FI,CGAS,NST)
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NS=3, NS1=NS+1,NGS=5)
C
C  NGS is the number of gases in the solar code. The order of gases in
C  this subroutine is: O2, H2O,CO2, O3, CH4. Note that CH4 is treated in
C  a funny way because Lisa simply replaced O3 with it in the main part
C  of the code. Consequently, we have to recopy it out from gas 3 to gas
C  5.
C

      DIMENSION T(ND),PF(ND),FI(NS1,ND),F(NGS,ND),FF(NGS,ND),
     2  CGAS(ND,NGS),CON(NGS), FNC(ND)  ! Added noncondensible mixing fraction array c-rr 5/14/2011 
      COMMON/CONSS/C,BK,G,PI,SM,DM,DM2
      DATA CON/1., 8.030E-4, 1.963E-3, 1., 1.e-5/
C
C        THIS SUBROUTINE CONVERTS VOLUME MIXING RATIOS OF THE VARIOUS
C   GASES INTO COLUMN CONCENTRATIONS WITHIN EACH LAYER.  O2 AND O3
C   COLUMN DEPTHS ARE EXPRESSED IN ATM CM; H2O AND CO2 ARE IN GM/CM2.
C
      ND1 = ND - 1
      BKMG = BK*273.15/(SM*G)
C



!      DO 6 J=1, ND
!   6  FNC(J) = 1 - FI(1,J) - FI(2,J)! Define noncondensible mixing ratio at each layer as FNC(J) c-rr 5/14/2011. 8/31/2011 removed this loop
     
c These are the two condensibles, water and CO2(gases 1 and 2 in FI). c-rr 5/14/2011
      DO 1 I=1,2
      DO 1 J=1,ND
   1  F(I+1,J) = FI(I,J)
      
C
C *** ELIMINATE O3 FROM THE SOLAR CODE AND PUT CH4 INTO GAS 5 SLOT ***
      DO 5 J = 1,ND
!      F(4,J) = FI(4,J)*FNC(J)   ! Ethane mixing ratio obtained by multiplying FI(4,J) by noncondensible mixing ratio c-rr 5/14/2011.. 8/31/2011 removed this
!   5  F(5,J) = FI(3,J)*FNC(J)   ! Methane mixing ratio obtained by multiplyling FI(3,J) by noncondensible mixing ratio c-rr 5/14/2011 8/31/2011 removed this
       F(4,J) = FI(4,J)  ! removing FNC
    5  F(5,J) = FI(3,J)  ! removing FNC

      
      
      DO 2 J=1,ND
   2  F(1,J) = FO2*(1. - FI(1,J)) !removing FNC
!    2  F(1,J)= FO2*FNC(J)  ! oxygen mixing ratio obtained by multiplying FO2 by noncondensible mixing ratio c-rr 5/14/2011  
  
      DO 3 I=1,NGS
      FF(I,1) = F(I,1)
      FF(I,ND) = F(I,ND)
      DO 3 J=2,ND1
   3  FF(I,J) = SQRT(F(I,J)*F(I,J-1))
C

      DO 4 I=1,NGS
      FACT = CON(I)*BKMG
      DO 4 J=1,ND1
   4  CGAS(J,I) = FACT * (PF(J+1) - PF(J)) * (FF(I,J) + 4.*F(I,J)
     2  + FF(I,J+1))/(6.*DM) 


!      print 200, ((CGAS(J,I), I=1,NGS), J=1, ND1)

! 200  format(1x, 1p5e12.5)
!      stop

 100  FORMAT(1X,1P10E12.5)
 
      RETURN
      END
