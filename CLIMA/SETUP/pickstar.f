      SUBROUTINE pickstar(STARR,SOLINT)
! Written by Ramses Ramirez 

! Stars after first 4 written in following manner: 
!The first letter are the Metallicities A,B,and C. These correspond to -0.5,+0.0, and +0.5, respectively. 
!The middle two numbers (30, 35, 40, 45, and 50)are logg. These correspond to 3.0, 3.5, 4.0, 4.5, and 5.0, respectively. 
!The final two numbers are the temperatures, 26 to 72, which correspond to 2600 to 7200K.'

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NSOL=38)
      DIMENSION SOLINT(NSOL), s(297) ! For 292 stars
      CHARACTER(5) :: STARR*5
      CHARACTER :: DIRDATA*10, DIRINOUT*8
      LOGICAL :: OPENINT
      real :: s, TOTAL 

      COMMON/DIR/DIRINOUT, DIRDATA

      ! Opening up Teal's muscles_stars.dat file

       OPEN(UNIT=82,FILE= DIRDATA//'/muscles_stars.dat')


!      read(10,*)  Removed labels so no need to skip line 
        do I =1,NSOL
          read(10,*)a,WL1,WL2,(s(k), k=1,297)
c          print *, read n,WL1,WL2,(s(k), k=1,296)
          IF (STARR== "Sun") THEN
             SOLINT(I)=s(1)
c             print *, SOLINT(I)
          ELSE IF (STARR=="GJ581")THEN
             SOLINT(I)=s(2)
!            print *, SOLINT(I)
          ELSE IF (STARR=="ADLEO")THEN
             SOLINT(I)=s(3)
!             print *, SOLINT(I)
          ELSE IF (STARR =="GJ644") THEN
             SOLINT(I)=s(4)
!             print *, SOLINT(I)

!------------------------------------ Starting -0.5 Metallicities

      ELSE IF (STARR =="A4026") THEN
             SOLINT(I)=s(5)
!            print *, SOLINT(I)
      ELSE IF (STARR =="A4526") THEN
             SOLINT(I)=s(6)
!            print *, SOLINT(I)
 



      ELSE IF (STARR =="A4028") THEN
             SOLINT(I)=s(7)
!             print *, SOLINT(I)
          ELSE IF (STARR =="A4528") THEN
              SOLINT(I)=s(8)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4030") THEN
             SOLINT(I)=s(9)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4530") THEN
             SOLINT(I)=s(10)
!             print *, SOLINT(I)                
 


 
      ELSE IF (STARR =="A4032") THEN
             SOLINT(I)=s(11)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4532") THEN
             SOLINT(I)=s(12)
!             print *, SOLINT(I)
 



      ELSE IF (STARR =="A4034") THEN
             SOLINT(I)=s(13)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4534") THEN
             SOLINT(I)=s(14)
!             print *, SOLINT(I)



      ELSE IF (STARR =="A4036") THEN
             SOLINT(I)=s(15)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4536") THEN
             SOLINT(I)=s(16)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4038") THEN
             SOLINT(I)=s(17)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4538") THEN
             SOLINT(I)=s(18)
!             print *, SOLINT(I)
 


      ELSE IF (STARR =="A4040") THEN
             SOLINT(I)=s(19)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4540") THEN
             SOLINT(I)=s(20)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4042") THEN
             SOLINT(I)=s(21)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4542") THEN
             SOLINT(I)=s(22)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4044") THEN
             SOLINT(I)=s(23)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4544") THEN
             SOLINT(I)=s(24)
!             print *, SOLINT(I)



      ELSE IF (STARR =="A4046") THEN
             SOLINT(I)=s(25)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4546") THEN
             SOLINT(I)=s(26)
!             print *, SOLINT(I)


      ELSE IF (STARR =="A4048") THEN
             SOLINT(I)=s(27)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4548") THEN
             SOLINT(I)=s(28)
!             print *, SOLINT(I)


      ELSE IF (STARR =="A4050") THEN
             SOLINT(I)=s(29)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4550") THEN
             SOLINT(I)=s(30)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4052") THEN
             SOLINT(I)=s(31)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4552") THEN
             SOLINT(I)=s(32)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4054") THEN
             SOLINT(I)=s(33)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4554") THEN
             SOLINT(I)=s(34)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4056") THEN
             SOLINT(I)=s(35)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4556") THEN
             SOLINT(I)=s(36)
!             print *, SOLINT(I)



      ELSE IF (STARR =="A4058") THEN
             SOLINT(I)=s(37)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4558") THEN
             SOLINT(I)=s(38)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4060") THEN
             SOLINT(I)=s(39)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4560") THEN
             SOLINT(I)=s(40)
!             print *, SOLINT(I)





      ELSE IF (STARR =="A4062") THEN
             SOLINT(I)=s(41)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4562") THEN
             SOLINT(I)=s(42)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4064") THEN
             SOLINT(I)=s(43)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4564") THEN
             SOLINT(I)=s(44)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4066") THEN
             SOLINT(I)=s(45)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4566") THEN
             SOLINT(I)=s(46)
!             print *, SOLINT(I)



      ELSE IF (STARR =="A4068") THEN
             SOLINT(I)=s(47)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4568") THEN
             SOLINT(I)=s(48)
!             print *, SOLINT(I)




      ELSE IF (STARR =="A4070") THEN
             SOLINT(I)=s(49)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4570") THEN
             SOLINT(I)=s(50)
!             print *, SOLINT(I)
 




      ELSE IF (STARR =="A4072") THEN
             SOLINT(I)=s(51)
!             print *, SOLINT(I)
      ELSE IF (STARR =="A4572") THEN
             SOLINT(I)=s(52)
!             print *, SOLINT(I)



!---------------------------------- For 0.0 Metallicities



      ELSE IF (STARR =="B3026") THEN
             SOLINT(I)=s(53)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3526") THEN
             SOLINT(I)=s(54)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4026") THEN
             SOLINT(I)=s(55)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4526") THEN
             SOLINT(I)=s(56)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5026") THEN
             SOLINT(I)=s(57)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3028") THEN
             SOLINT(I)=s(58)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3528") THEN
             SOLINT(I)=s(59)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4028") THEN
             SOLINT(I)=s(60)
!             print *, SOLINT(I)
          ELSE IF (STARR =="B4528") THEN
              SOLINT(I)=s(61)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5028") THEN
             SOLINT(I)=s(62)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3030") THEN
             SOLINT(I)=s(63)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3530") THEN
             SOLINT(I)=s(64)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4030") THEN
             SOLINT(I)=s(65)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4530") THEN
             SOLINT(I)=s(66)
!             print *, SOLINT(I)                
      ELSE IF (STARR =="B5030") THEN
             SOLINT(I)=s(67)
!             print *, SOLINT(I)                 


      ELSE IF (STARR =="B3032") THEN
             SOLINT(I)=s(68)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3532") THEN
             SOLINT(I)=s(69)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4032") THEN
             SOLINT(I)=s(70)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4532") THEN
             SOLINT(I)=s(71)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5032") THEN
             SOLINT(I)=s(72)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3034") THEN
             SOLINT(I)=s(73)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3534") THEN
             SOLINT(I)=s(74)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4034") THEN
             SOLINT(I)=s(75)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4534") THEN
             SOLINT(I)=s(76)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5034") THEN
             SOLINT(I)=s(77)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3036") THEN
             SOLINT(I)=s(78)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3536") THEN
             SOLINT(I)=s(79)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4036") THEN
             SOLINT(I)=s(80)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4536") THEN
             SOLINT(I)=s(81)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5036") THEN
             SOLINT(I)=s(82)
!             print *, SOLINT(I)      


      ELSE IF (STARR =="B3038") THEN
             SOLINT(I)=s(83)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3538") THEN
             SOLINT(I)=s(84)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4038") THEN
             SOLINT(I)=s(85)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4538") THEN
             SOLINT(I)=s(86)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5038") THEN
             SOLINT(I)=s(87)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3040") THEN
             SOLINT(I)=s(88)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3540") THEN
             SOLINT(I)=s(89)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4040") THEN
             SOLINT(I)=s(90)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4540") THEN
             SOLINT(I)=s(91)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5040") THEN
             SOLINT(I)=s(92)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3042") THEN
             SOLINT(I)=s(93)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3542") THEN
             SOLINT(I)=s(94)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4042") THEN
             SOLINT(I)=s(95)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4542") THEN
             SOLINT(I)=s(96)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5042") THEN
             SOLINT(I)=s(97)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3044") THEN
             SOLINT(I)=s(98)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3544") THEN
             SOLINT(I)=s(99)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4044") THEN
             SOLINT(I)=s(100)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4544") THEN
             SOLINT(I)=s(101)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5044") THEN
             SOLINT(I)=s(102)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3046") THEN
             SOLINT(I)=s(103)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3546") THEN
             SOLINT(I)=s(104)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4046") THEN
             SOLINT(I)=s(105)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4546") THEN
             SOLINT(I)=s(106)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5046") THEN
             SOLINT(I)=s(107)
!             print *, SOLINT(I)

      ELSE IF (STARR =="B3048") THEN
             SOLINT(I)=s(108)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3548") THEN
             SOLINT(I)=s(109)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4048") THEN
             SOLINT(I)=s(110)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4548") THEN
             SOLINT(I)=s(111)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5048") THEN
             SOLINT(I)=s(112)
!             print *, SOLINT(I)

      ELSE IF (STARR =="B3050") THEN
             SOLINT(I)=s(113)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3550") THEN
             SOLINT(I)=s(114)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4050") THEN
             SOLINT(I)=s(115)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4550") THEN
             SOLINT(I)=s(116)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5050") THEN
             SOLINT(I)=s(117)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3052") THEN
             SOLINT(I)=s(118)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3552") THEN
             SOLINT(I)=s(119)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4052") THEN
             SOLINT(I)=s(120)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4552") THEN
             SOLINT(I)=s(121)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5052") THEN
             SOLINT(I)=s(122)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3054") THEN
             SOLINT(I)=s(123)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3554") THEN
             SOLINT(I)=s(124)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4054") THEN
             SOLINT(I)=s(125)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4554") THEN
             SOLINT(I)=s(126)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5054") THEN
             SOLINT(I)=s(127)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3056") THEN
             SOLINT(I)=s(128)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3556") THEN
             SOLINT(I)=s(129)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4056") THEN
             SOLINT(I)=s(130)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4556") THEN
             SOLINT(I)=s(131)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5056") THEN
             SOLINT(I)=s(132)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3058") THEN
             SOLINT(I)=s(133)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3558") THEN
             SOLINT(I)=s(134)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4058") THEN
             SOLINT(I)=s(135)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4558") THEN
             SOLINT(I)=s(136)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5058") THEN
             SOLINT(I)=s(137)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3060") THEN
             SOLINT(I)=s(138)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3560") THEN
             SOLINT(I)=s(139)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4060") THEN
             SOLINT(I)=s(140)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4560") THEN
             SOLINT(I)=s(141)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5060") THEN
             SOLINT(I)=s(142)
!             print *, SOLINT(I)



      ELSE IF (STARR =="B3062") THEN
             SOLINT(I)=s(143)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3562") THEN
             SOLINT(I)=s(144)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4062") THEN
             SOLINT(I)=s(145)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4562") THEN
             SOLINT(I)=s(146)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5062") THEN
             SOLINT(I)=s(147)
!             print *, SOLINT(I)



      ELSE IF (STARR =="B3064") THEN
             SOLINT(I)=s(148)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3564") THEN
             SOLINT(I)=s(149)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4064") THEN
             SOLINT(I)=s(150)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4564") THEN
             SOLINT(I)=s(151)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5064") THEN
             SOLINT(I)=s(152)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3066") THEN
             SOLINT(I)=s(153)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3566") THEN
             SOLINT(I)=s(154)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4066") THEN
             SOLINT(I)=s(155)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4566") THEN
             SOLINT(I)=s(156)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5066") THEN
             SOLINT(I)=s(157)
!             print *, SOLINT(I)

      ELSE IF (STARR =="B3068") THEN
             SOLINT(I)=s(158)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3568") THEN
             SOLINT(I)=s(159)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4068") THEN
             SOLINT(I)=s(160)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4568") THEN
             SOLINT(I)=s(161)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5068") THEN
             SOLINT(I)=s(162)
!             print *, SOLINT(I)


      ELSE IF (STARR =="B3070") THEN
             SOLINT(I)=s(163)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3570") THEN
             SOLINT(I)=s(164)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4070") THEN
             SOLINT(I)=s(165)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4570") THEN
             SOLINT(I)=s(166)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5070") THEN
             SOLINT(I)=s(167)
!             print *, SOLINT(I)



      ELSE IF (STARR =="B3072") THEN
             SOLINT(I)=s(168)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B3572") THEN
             SOLINT(I)=s(169)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4072") THEN
             SOLINT(I)=s(170)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B4572") THEN
             SOLINT(I)=s(171)
!             print *, SOLINT(I)
      ELSE IF (STARR =="B5072") THEN
             SOLINT(I)=s(172)
!             print *, SOLINT(I)
       ELSE IF (STARR =="B4574") THEN
             SOLINT(I)=s(293)
          
       ELSE IF (STARR =="B4576") THEN
             SOLINT(I)=s(294)

       ELSE IF (STARR =="B4578") THEN
             SOLINT(I)=s(295)

       ELSE IF (STARR =="B4580") THEN
             SOLINT(I)=s(296)
!------------------------------------- +0.5 Metalliticities

      ELSE IF (STARR =="C3026") THEN
             SOLINT(I)=s(173)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3526") THEN
             SOLINT(I)=s(174)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4026") THEN
             SOLINT(I)=s(175)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4526") THEN
             SOLINT(I)=s(176)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5026") THEN
             SOLINT(I)=s(177)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3028") THEN
             SOLINT(I)=s(178)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3528") THEN
             SOLINT(I)=s(179)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4028") THEN
             SOLINT(I)=s(180)
!             print *, SOLINT(I)
          ELSE IF (STARR =="C4528") THEN
              SOLINT(I)=s(181)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5028") THEN
             SOLINT(I)=s(182)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3030") THEN
             SOLINT(I)=s(183)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3530") THEN
             SOLINT(I)=s(184)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4030") THEN
             SOLINT(I)=s(185)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4530") THEN
             SOLINT(I)=s(186)
!             print *, SOLINT(I)                
      ELSE IF (STARR =="C5030") THEN
             SOLINT(I)=s(187)
!             print *, SOLINT(I)                 


      ELSE IF (STARR =="C3032") THEN
             SOLINT(I)=s(188)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3532") THEN
             SOLINT(I)=s(189)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4032") THEN
             SOLINT(I)=s(190)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4532") THEN
             SOLINT(I)=s(191)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5032") THEN
             SOLINT(I)=s(192)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3034") THEN
             SOLINT(I)=s(193)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3534") THEN
             SOLINT(I)=s(194)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4034") THEN
             SOLINT(I)=s(195)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4534") THEN
             SOLINT(I)=s(196)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5034") THEN
             SOLINT(I)=s(197)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3036") THEN
             SOLINT(I)=s(198)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3536") THEN
             SOLINT(I)=s(199)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4036") THEN
             SOLINT(I)=s(200)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4536") THEN
             SOLINT(I)=s(201)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5036") THEN
             SOLINT(I)=s(202)
!             print *, SOLINT(I)      


      ELSE IF (STARR =="C3038") THEN
             SOLINT(I)=s(203)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3538") THEN
             SOLINT(I)=s(204)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4038") THEN
             SOLINT(I)=s(205)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4538") THEN
             SOLINT(I)=s(206)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5038") THEN
             SOLINT(I)=s(207)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3040") THEN
             SOLINT(I)=s(208)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3540") THEN
             SOLINT(I)=s(209)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4040") THEN
             SOLINT(I)=s(210)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4540") THEN
             SOLINT(I)=s(211)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5040") THEN
             SOLINT(I)=s(212)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3042") THEN
             SOLINT(I)=s(213)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3542") THEN
             SOLINT(I)=s(214)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4042") THEN
             SOLINT(I)=s(215)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4542") THEN
             SOLINT(I)=s(216)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5042") THEN
             SOLINT(I)=s(217)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3044") THEN
             SOLINT(I)=s(218)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3544") THEN
             SOLINT(I)=s(219)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4044") THEN
             SOLINT(I)=s(220)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4544") THEN
             SOLINT(I)=s(221)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5044") THEN
             SOLINT(I)=s(222)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3046") THEN
             SOLINT(I)=s(223)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3546") THEN
             SOLINT(I)=s(224)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4046") THEN
             SOLINT(I)=s(225)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4546") THEN
             SOLINT(I)=s(226)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5046") THEN
             SOLINT(I)=s(227)
!             print *, SOLINT(I)

      ELSE IF (STARR =="C3048") THEN
             SOLINT(I)=s(228)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3548") THEN
             SOLINT(I)=s(229)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4048") THEN
             SOLINT(I)=s(230)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4548") THEN
             SOLINT(I)=s(231)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5048") THEN
             SOLINT(I)=s(232)
!             print *, SOLINT(I)

      ELSE IF (STARR =="C3050") THEN
             SOLINT(I)=s(233)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3550") THEN
             SOLINT(I)=s(234)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4050") THEN
             SOLINT(I)=s(235)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4550") THEN
             SOLINT(I)=s(236)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5050") THEN
             SOLINT(I)=s(237)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3052") THEN
             SOLINT(I)=s(238)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3552") THEN
             SOLINT(I)=s(239)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4052") THEN
             SOLINT(I)=s(240)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4552") THEN
             SOLINT(I)=s(241)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5052") THEN
             SOLINT(I)=s(242)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3054") THEN
             SOLINT(I)=s(243)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3554") THEN
             SOLINT(I)=s(244)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4054") THEN
             SOLINT(I)=s(245)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4554") THEN
             SOLINT(I)=s(246)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5054") THEN
             SOLINT(I)=s(247)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3056") THEN
             SOLINT(I)=s(248)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3556") THEN
             SOLINT(I)=s(249)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4056") THEN
             SOLINT(I)=s(250)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4556") THEN
             SOLINT(I)=s(251)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5056") THEN
             SOLINT(I)=s(252)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3058") THEN
             SOLINT(I)=s(253)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3558") THEN
             SOLINT(I)=s(254)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4058") THEN
             SOLINT(I)=s(255)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4558") THEN
             SOLINT(I)=s(256)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5058") THEN
             SOLINT(I)=s(257)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3060") THEN
             SOLINT(I)=s(258)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3560") THEN
             SOLINT(I)=s(259)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4060") THEN
             SOLINT(I)=s(260)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4560") THEN
             SOLINT(I)=s(261)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5060") THEN
             SOLINT(I)=s(262)
!             print *, SOLINT(I)



      ELSE IF (STARR =="C3062") THEN
             SOLINT(I)=s(263)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3562") THEN
             SOLINT(I)=s(264)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4062") THEN
             SOLINT(I)=s(265)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4562") THEN
             SOLINT(I)=s(266)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5062") THEN
             SOLINT(I)=s(267)
!             print *, SOLINT(I)



      ELSE IF (STARR =="C3064") THEN
             SOLINT(I)=s(268)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3564") THEN
             SOLINT(I)=s(269)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4064") THEN
             SOLINT(I)=s(270)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4564") THEN
             SOLINT(I)=s(271)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5064") THEN
             SOLINT(I)=s(272)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3066") THEN
             SOLINT(I)=s(273)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3566") THEN
             SOLINT(I)=s(274)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4066") THEN
             SOLINT(I)=s(275)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4566") THEN
             SOLINT(I)=s(276)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5066") THEN
             SOLINT(I)=s(277)
!             print *, SOLINT(I)

      ELSE IF (STARR =="C3068") THEN
             SOLINT(I)=s(278)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3568") THEN
             SOLINT(I)=s(279)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4068") THEN
             SOLINT(I)=s(280)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4568") THEN
             SOLINT(I)=s(281)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5068") THEN
             SOLINT(I)=s(282)
!             print *, SOLINT(I)


      ELSE IF (STARR =="C3070") THEN
             SOLINT(I)=s(283)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3570") THEN
             SOLINT(I)=s(284)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4070") THEN
             SOLINT(I)=s(285)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4570") THEN
             SOLINT(I)=s(286)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5070") THEN
             SOLINT(I)=s(287)
!             print *, SOLINT(I)



      ELSE IF (STARR =="C3072") THEN
             SOLINT(I)=s(288)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C3572") THEN
             SOLINT(I)=s(289)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4072") THEN
             SOLINT(I)=s(290)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C4572") THEN
             SOLINT(I)=s(291)
!             print *, SOLINT(I)
      ELSE IF (STARR =="C5072") THEN
             SOLINT(I)=s(292)
!             print *, SOLINT(I)
      ELSE IF (STARR =="M5V") THEN
             SOLINT(I)=s(297)
             print *, SOLINT(I)

C    Here on out is Teal's muscles spectra additions... -- Teal
      ELSE
            READ(82,*) ll,x80,x81,x82,x83,x84,x85,x86,x87,x88,
     +      x89,x90,x,x
            IF(STARR=="G80") SOLINT(I) = x80
            IF(STARR=="G81") SOLINT(I) = x81
            IF(STARR=="G82") SOLINT(I) = x82
            IF(STARR=="G83") SOLINT(I) = x83
            IF(STARR=="G84") SOLINT(I) = x84
            IF(STARR=="G85") SOLINT(I) = x85
            IF(STARR=="G86") SOLINT(I) = x86
            IF(STARR=="G87") SOLINT(I) = x87
            IF(STARR=="G88") SOLINT(I) = x88
            IF(STARR=="G89") SOLINT(I) = x89
            IF(STARR=="G90") SOLINT(I) = x90
            IF(STARR=="G91") SOLINT(I) = x91

      ENDIF

        enddo
       rewind(10)
!      close(10)
!c             DO I=1,38
!c              print *, SOLINT(I)
!c         ENDDO
        TOTAL=0.
        DO I=1,NSOL
           TOTAL=SOLINT(I)+TOTAL
        ENDDO
c        PRINT *,'TOTALFLUX=',TOTAL

        END
