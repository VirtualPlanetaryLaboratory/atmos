      program convdist

      implicit real*8(A-H,O-Z)
      character*15 fmtstr

      PARAMETER(NZ=200,NQ=62,NSP=56,NP=4)   !NSP doesn't matter here

      dimension USOL2(NQ,NZ),T2(NZ),USOL3(NQ-1,NZ),USOL4(NQ+1,NZ)
      dimension USOL5(NQ+NP,NZ)
      dimension USOL(NQ,NZ),T(NZ),EDD(NZ),DEN(NZ),O3(NZ),SL(NSP,NZ)
      dimension AERSOL(NZ,NP),WFALL(NZ,NP),RPAR(NZ,NP) 
      dimension SO4AER(NZ),S8AER(NZ), PARTICLES(NZ,NP)

C ***** READ THE INPUT DATAFILE *****

      !option 2 is to read in new style, delete one column
      !option 3 is to read in new style, swap two columns
      !option 4 is read in new style, multiply a column by a factor and add in
      !option 5 double the grid 
      !option 6 halve the grid
      !option 7 converts particles from TD to main loop

!make the above a user choice

      PRINT *, ''
      PRINT *, 'Which option would you like to run?'
      PRINT *, ''
      PRINT *, '2 = read in new style, delete one column'
      PRINT *, '3 = read in new style, swap two columns'
      PRINT *, '4 = read in new style, multiply a column'
      PRINT *, '    by a factor and add in'
      PRINT *, '5 = double the number of lines in the grid'
      PRINT *, '6 = halve the number of lines in the grid'
      PRINT *, '7 = convert particles from TD to main loop' 
      PRINT *, ''
      PRINT *, ''
      READ *, option   !assumes integer
      Moption = option !integer value of option for neater output
      
      PRINT *, ''
      PRINT *, 'You have selected option: ', Moption
      PRINT *, ''
      
      IF (option.EQ.2) THEN
         PRINT *, 'Which Column would you like to delete?'
         PRINT *, ''
         READ *, kskip
         Mkskip = kskip
         PRINT *, ''
         PRINT *, 'Column to delete: ', Mkskip
      END IF

      IF (option.EQ.3) THEN
         PRINT *, 'Enter the first column to be swapped: '
         PRINT *, ''
         READ *, Kswap1 
         Mkswap1 = Kswap1
         PRINT *, ''
         PRINT *, 'Frist column to be swapped: ', Mkswap1
         PRINT *, ''
         PRINT *, 'Enter the second column to be swapped: '
         PRINT *, ''
         READ *, Kswap2
         Mkswap2 = Kswap2
         PRINT *, ''
         PRINT *, 'Second column to be swapped: ', Mkswap2
      END IF
         
      IF (option.EQ.4) THEN
         PRINT *, 'Which column would you like to copy?'
         PRINT *, ''
         READ *, coltocopy      !assumes acceptable column value
         Mcoltocopy = coltocopy
         PRINT *, ''
         PRINT *, 'Column to be copied = ', Mcoltocopy
         PRINT *, ''
         PRINT *, 'Please enter multiplication factor: '
         PRINT *, ''
         READ *, factor
         PRINT *, ''
         PRINT *, 'Multiplication factor = ', factor
         PRINT *, ''
         PRINT *, 'Which column would you like to insert?'
         READ *, Mcoltoinsert
         PRINT *, 'new scaled column to be added in  = ', Mcoltocopy
         PRINT *, ''
         

      END IF

      if (moption.ne.7) then
      PRINT *, ''                                          
      PRINT *, 'Does the out.dist file have particles in the'
      PRINT *, '      tri-diagonal format rather than in the main loop?'
      PRINT *, 'Please enter 0 for no, 1 for yes'
      PRINT *, ''
      READ *, USETD !assumes either 0 or 1 will be entered
      MUSETD = USETD
      PRINT *, ''
      else  !option 7
         print *, 'converting tri-diag to main loop'
        USETD=1 
      MUSETD = USETD
      endif   

      PRINT *, 'You have chosen ', MUSETD


   !read in the new format file

      if (option.ge.2) then
     
      open(18, file='out.dist',status='OLD') 

      IROW = 10  !num columns in .dist file
      LR = NQ/IROW + 1      !NQ=80  -> 9
      RL = FLOAT(NQ)/IROW + 1  !9
      DIF = RL - LR  !0
      IF (DIF.LT.0.001) LR = LR - 1  !so LR=8
C

      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.LT.LR) then
        read(18, 880) ((USOL(k,i),K=K1,K2),i=1,nz) 
       ELSE
          K2 = NQ
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10
            read(18, 880) ((USOL(k,i),K=K1,K2),i=1,nz) 
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           read(18, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
          endif    
       ENDIF
      enddo 


        read(18,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz) !this final one is CO2 number density

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP*3

        do i=1,nz
         read (18,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
     $                  (RPAR(i,j),j=1,NP) 
        enddo  

      if(USETD.EQ.1) then !particles in tri-diag
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP
      do i=1,nz
       read(18,fmtstr)  (PARTICLES(i,j),j=1,np)  !ordering is that in species.dat
      enddo
      endif

      endif   !end read in of new style 10 column .dist file


! done with reading in, now move to writing out...
!this is option 1 stuff - needs a bit of re-testing

 

      NQwrite=NQ
      NPwrite=NP
      USOL2=USOL

      if (option.eq.2) then
         NQwrite=NQ-1

C        Kskip=31
         m=1

         do i=1,nq
          if (i.ne.kskip) then 
c             print *, i,m
             do j=1,nz
             USOL2(m,j)=USOL(i,j)
           enddo   
            m=m+1  
          endif
         enddo    
      endif   !end option 2


      if (option.eq.3) then
         NQwrite=NQ

c        Kswap1=1
c        Kswap2=1
         m=1

         do i=1,nq
         
            if (i.eq.Kswap1) m=Kswap2
             
            if (i.eq.Kswap2) m=Kswap1
           do j=1,nz
             USOL2(m,j)=USOL(i,j)
           enddo   
            m=i+1  
         enddo    
      endif   



      if (option.eq.4) then
         NQwrite=NQ+1
         print *, ''
         print *, 'new column added: ',Mcoltoinsert
         inc=0
         do i=1,nq+1
           do j=1,nz
              if (i.eq.Mcoltoinsert) then
               USOL4(i,j)=USOL(INT(coltocopy),j)*factor
               inc=1
              else   
               USOL4(i,j)=USOL(i-inc,j)
              endif
           enddo   
         enddo    
      endif   !end option 4


      if (option.eq.7) then
         NQwrite=NQ+NP
         NPwrite=NP
         print *, ''
         USETD=0
         do i=1,NQ+NP
           do j=1,nz
              if (i.gt.nq) then
               USOL5(i,j)=PARTICLES(j,i-NQ)
              else   
               USOL5(i,j)=USOL(i,j)
              endif
           enddo   
         enddo    
      endif   !end option 7

 880  format(10E17.8)
 881  format(5E17.8)

      IF (option.le.4 .or. option.eq.7) THEN

      open(88, file='out.distNEW',status='UNKNOWN')    ! formatted output

      IROW = 10  !num columns in .dist file
      LR = NQwrite/IROW + 1      
      RL = FLOAT(NQwrite)/IROW + 1  
      DIF = RL - LR 
      IF (DIF.LT.0.001) LR = LR - 1  


      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.lt.LR) then
        if (option.le.3) then  
         write(88, 880) ((USOL2(k,i),K=K1,K2),i=1,nz) 
         else if (option.eq.7) then
        write(88, 880) ((USOL5(k,i),K=K1,K2),i=1,nz) 
         else   
        write(88, 880) ((USOL4(k,i),K=K1,K2),i=1,nz) 
        endif
       ELSE
          K2 = NQWrite
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10

           if (option.le.3) then  
            write(88, 880) ((USOL2(k,i),K=K1,K2),i=1,nz) 
           else if (option.eq.7) then
            write(88, 880) ((USOL5(k,i),K=K1,K2),i=1,nz) 
           else
            write(88, 880) ((USOL4(k,i),K=K1,K2),i=1,nz) 
           endif

          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10

           if (option.le.3) then  
            write(88, fmtstr) ((USOL2(k,i),K=K1,K2),i=1,nz) 
           else if (option.eq.7) then
            write(88, fmtstr) ((USOL5(k,i),K=K1,K2),i=1,nz) 
           else
            write(88, fmtstr) ((USOL4(k,i),K=K1,K2),i=1,nz) 
           endif
          endif
       ENDIF
      enddo 



        write (88,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz) !this final one is CO2 number density


        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NPwrite*3   
        do i=1,nz
         write(88,fmtstr) (AERSOL(i,j),j=1,NPwrite),(WFALL(i,j),
     $     j=1,NPwrite),(RPAR(i,j),j=1,NPwrite) 
        enddo  

        if(USETD.EQ.1) then
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NPwrite   
        
      do i=1,nz
       write(88,fmtstr)  (PARTICLES(i,j),j=1,npwrite)  !ordering is that in species.dat
      enddo

      endif

      ENDIF

C     Start to add new code for doubling and halving

      IF (option.eq.5) THEN

      open(88, file='out.distNEW',status='UNKNOWN') ! formatted output

      IROW = 10  !num columns in .dist file
      LR = NQwrite/IROW + 1      
      RL = FLOAT(NQwrite)/IROW + 1  
      DIF = RL - LR 
      IF (DIF.LT.0.001) LR = LR - 1  


      
      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
 
      IF (L.lt.LR) then
         DO i=1,nz
           write(88, 880) (USOL(k,i),K=K1,K2)
           write(88, 880) (USOL(k,i),K=K1,K2)
         ENDDO  
       ELSE
          K2 = NQWrite
    
          if (K2-K1 .EQ. 9) then !this only occurs if NQ is a multiple of 10
            DO i=1,nz
               write(88, 880) (USOL(k,i),K=K1,K2)
               write(88, 880) (USOL(k,i),K=K1,K2)
            ENDDO
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1 !OK for one character as this should always be <10
           DO i=1,nz
              write(88, fmtstr) (USOL(k,i),K=K1,K2)  
              write(88, fmtstr) (USOL(k,i),K=K1,K2)  
           ENDDO
          endif
       ENDIF
      enddo 

      
        DO i=1,nz
           write (88,881) T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i) !this final one is CO2 number density
           write (88,881) T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i)
        ENDDO
   
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NPwrite*3   
        do i=1,nz
           write(88,fmtstr) (AERSOL(i,j),j=1,NPwrite),(WFALL(i,j),
     $          j=1,NPwrite),(RPAR(i,j),j=1,NPwrite) 
           write(88,fmtstr) (AERSOL(i,j),j=1,NPwrite),(WFALL(i,j),
     $          j=1,NPwrite),(RPAR(i,j),j=1,NPwrite) 

        enddo  

        if(USETD.EQ.1) then
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NPwrite   
        
        do i=1,nz
           write(88,fmtstr)  (PARTICLES(i,j),j=1,npwrite) !ordering is that in species.dat
           write(88,fmtstr)  (PARTICLES(i,j),j=1,npwrite) !ordering is that in species.dat
        enddo

      endif


      END IF




      IF (option.eq.6) THEN

      open(88, file='out.distNEW',status='UNKNOWN') ! formatted output

      IROW = 10  !num columns in .dist file
      LR = NQwrite/IROW + 1      
      RL = FLOAT(NQwrite)/IROW + 1  
      DIF = RL - LR 
      IF (DIF.LT.0.001) LR = LR - 1  


      
      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
 
      IF (L.lt.LR) then
          write(88, 880) ((USOL(k,i),K=K1,K2),i=1,nz,2)
       ELSE
          K2 = NQWrite
    
          if (K2-K1 .EQ. 9) then !this only occurs if NQ is a multiple of 10
             write(88, 880) ((USOL(k,i),K=K1,K2),i=1,nz,2)
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1 !OK for one character as this should always be <10
           write(88, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz,2)         
          endif
       ENDIF
      enddo 



        write (88,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz,2) !this final one is CO2 number density


        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NPwrite*3   
        do i=1,nz,2
         write(88,fmtstr) (AERSOL(i,j),j=1,NPwrite),(WFALL(i,j),
     $     j=1,NPwrite),(RPAR(i,j),j=1,NPwrite) 
        enddo  

        if(USETD.EQ.1) then
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NPwrite   
        
      do i=1,nz,2
       write(88,fmtstr)  (PARTICLES(i,j),j=1,npwrite)  !ordering is that in species.dat
      enddo

      endif


      END IF


      print *, ''
      print *, 'requested changes are in out.distNEW'


        END
