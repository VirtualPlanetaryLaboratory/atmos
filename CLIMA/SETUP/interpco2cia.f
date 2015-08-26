        SUBROUTINE interpco2cia(T,IL,IMS,IMS1,AFX)
c       In Fortran both input and output arguments are put in function and CALLS too!                
c        Written by Ramses Ramirez
c        Dont forget the include statements throughout the program files!
        include 'CLIMA/INCLUDE/header.inc'
        PARAMETER(NF=55)

        
        real tempgrid(7), FX(ND), tempr, Tnum, Tden
        dimension MS(ND), MS1(ND)
        
        data tempgrid/100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0/


        tempr=T  ! Temperature around which CIAs are interpolated.This will be done within first altitude loop and following
! will be done for all 100 layers (100 temperatures). MS, MS1, and FX will be stored for each layer (101 element array vectors for each).
! MS, MS1, and FX are not functions of frequency (I), only of height.
                        
        DO M=1,7
                
                MS(IL) = M
                if (tempr.lt.tempgrid(M))then
                       
                        exit
                else
                MS(IL) = 8
                endif
        ENDDO
       
       

        if (MS(IL) .le. 1) then   ! MS1 is left-most temperature grid point. MS is right-most temperature grid pointc
           MS1(IL)=MS(IL)
           FX(IL)=1   
          
        elseif (MS(IL) .ge. 8) then
           MS(IL) = MS(IL) - 1
           MS1(IL) = MS(IL)
           FX(IL) =1 
          
        else
          
           MS1(IL) = MS(IL)-1
           
                    
           Tnum = tempr - tempgrid(MS1(IL))  ! I am assuming temperature is not on a log scale? Otherwise Tnum and Tden are logged.
           Tden = tempgrid(MS(IL)) - tempgrid(MS1(IL))
          
           
           FX(IL) = Tnum/Tden
            
        endif
 
        IMS = MS(IL)
        IMS1 = MS1(IL)
        AFX  = FX(IL)
!        print *, T, MS(IL), MS1(IL), IL
!        pause
         
          ! Put these statements under second altitude loop (which is under a wavelength loop). Will use stored MS,MS1, and FX values at right altitude to get CIA at the right altitude.
! Might need another height variable for the CIAs like CIA(IL, MS1,I)
! I is the wavelength loop (max NF), IL is the height loop, and MS1(IL) and MS(IL) are the left-most and right-most temperatures at the temperature grid point, respectively. MS and MS1 have 100 values each but they only range from 1-7.


c        CIAMS1L = alog(CIA(IL,MS1(IL),I))
c        CIAMSL = alog(CIA(IL,MS(IL),I))
c        CPRL = CIAMS1L + FX(IL)*(CIAMSL-CIAMS1L)
c        CPRW(IL,I) = exp(CPRL)        

        END
