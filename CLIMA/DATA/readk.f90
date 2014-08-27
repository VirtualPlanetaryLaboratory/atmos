PROGRAM READK
! reads separate CO2 and H2O coefficients. 
! Compile as : ifort -r8 readk.f90
IMPLICIT NONE

integer ::  i,j, k,t,p,solar_or_ir
real :: a, kappa_sol
integer, parameter :: intervals = 55 ! 38 for solar, 55 for IR
integer, parameter :: IK = 8 ! 8 total terms
dimension :: kappa_sol(intervals,8,8,IK)
 
!--------------------------------------------------------------

! Initializing k-coefficient array
do i = 1,intervals
         do t = 1, 8 ! 8 temperatures
            do p = 1,8  ! 8 pressures
               do k = 1,IK
            kappa_sol(i,t,p,k) = 1.e-60
               enddo ! k loop
           enddo ! ends pressure loop
	enddo ! ends temperature loop
enddo  ! ends interval loop


solar_or_ir = 1


!open(solar_or_ir, file='solar_38_H2O.dat')
open(solar_or_ir, file='ir_55_H2O.dat')

do j = 1, 16
read(solar_or_ir, *) ! skips 16 lines to get to where it is repeatable
enddo


do i = 1,intervals

         do t = 1, 8 ! 8 temperatures
            read(solar_or_ir,*) ! skip 3 more lines to read data
            read(solar_or_ir,*)
            read(solar_or_ir,*)
            do p = 1,8  ! 8 pressures
            read(solar_or_ir,*)a,(kappa_sol(i,t,p,k),k=1,IK)

           enddo ! ends pressure loop
	enddo ! ends temperature loop
enddo  ! ends interval loop



! ------ Write out k-coefficients test---------

do i = 1,intervals
         do t = 1, 8 ! 8 temperatures
            do p = 1,8  ! 8 pressures
               
            write(60,*)(kappa_sol(i,t,p,k), k =1,IK)
            ! pause
           enddo ! ends pressure loop
            write(60,*)
	enddo ! ends temperature loop
enddo  ! ends interval loop


END ! ENDS PROGRAM
