        SUBROUTINE interpo2cia(T,IL,IMS,IMS1,AFX)
c       In Fortran both input and output arguments are put in function and CALLS too!                
c        Written by Ramses Ramirez
c        Don't forget the include statements throughout the program files!
        include 'CLIMA/INCLUDE/header.inc'
        PARAMETER(NF=55)


        real tempgrid(15), FX(ND), tempr, Tnum, Tden
        dimension MS(ND), MS1(ND)
        
        data tempgrid/193.4, 206.8, 218.6, 220.8, 229.4, 229.6,240.0,
     &               249.4,270.7,297.5,297.8,320.5,340.8,344.8,353.4/


        tempr=T  ! Temperature around which N2-H2 CIA self broadening coefficients are interpolated.This will be done within first altitude loop and following
! will be done for all 100 layers (100 temperatures). MS, MS1, and FX will be stored for each layer (101 element array vectors for each).
! MS, MS1, and FX are not functions of frequency (I), only of height.
                        
        DO M=1,15
                
                MS(IL) = M
                if (tempr.lt.tempgrid(M))then
                       
                        exit
                else
                MS(IL) = 16
                endif
        ENDDO
       
       

        if (MS(IL) .le. 1) then   ! MS1 is left-most temperature grid point. MS is right-most temperature grid pointc
           MS1(IL)=MS(IL)
           FX(IL)=1   
          
        elseif (MS(IL) .ge. 16) then
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
         

         END
