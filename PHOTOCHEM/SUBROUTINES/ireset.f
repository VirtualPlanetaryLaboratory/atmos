        SUBROUTINE RESET
        INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
        CHARACTER*20 fmtstr
        DIMENSION USOL(NQ,NZ)
        DIMENSION T(NZ),EDD(NZ),DEN(NZ),O3(NZ),SL(NSP,NZ)
        DIMENSION AERSOL(NZ,NP),WFALL(NZ,NP),RPAR(NZ,NP)
        DIMENSION PARTICLES(NZ,NP)
        COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FH2
C  Using IRESET as a guide to overwrite mixing_ratios.dat with modern earth quantities
        open(69, file=
     &  'PHOTOCHEM/INPUTFILES/templates/modernEarthS/mixing_ratios.dat')

         READ(69,102) FAR, FCH4, FC2H6, FCO2,
     &                  FN2, FO2, FH2, FNO2, JCOLD
  102    FORMAT(1PE10.3,10x,/,1PE10.3,10x,/,
     &          1PE10.3,10x,/,1PE10.3,10x,/,
     &          1PE10.3,10x,/,1PE10.3,10x,/,
     &          1PE10.3,10x,/,1PE10.3,10x,
     &          /,I3,17x)


        DO J=1,NZ
        DO I=1,NQ
        if I=LCO2 then
          USOL(I,J)=FCO2
         else if I=LAR then
          USOL(I,J)=FAR
         else if I=LC2H6 then
          USOL(I,J)=FC2H6
         else if I=LCH4 then
          USOL(I,J)=FCH4
         else if I=LN2 then
          USOL(I,J)=FN2
         else if I=LNO2 then
          USOL(I,J)=FNO2
        ENDDO
        ENDDO

        READ(66,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz)
 881  format(5E17.8)

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP*3
        do i=1,nz
         READ(67,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
     $                  (RPAR(i,j),j=1,NP)           ! print aerosols into another file .aersol
        enddo

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP

      if(USETD.EQ.1) then
      do i=1,nz
       READ(68,fmtstr)  (PARTICLES(i,j),j=1,np)  !print tridag into another file .tridag
      enddo
      endif



      END
