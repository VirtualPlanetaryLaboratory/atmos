
      SUBROUTINE RAYLEY(SIGR,AL2)
C
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NSOL=38, NGS=5)
      PARAMETER(NS=3, NS1=NS+1, NS4=NS+5) ! Adding parameter statement needed for FI(NS1,ND) 5/23/2011	  
C
       COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1)
C
C
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,Fwater, FO2NC, 
     2 FARNC, FI(NS1,ND) ! Don't need to pass FCH4NC nor FNC because already included in FI(NS1,ND) 5/3/2011
      COMMON/CONSS/C,BK,G,PI,SM,DM, DM2  ! Passing DM2 which is AMN2 in CONVEC. 5/3/2011
      DIMENSION A(4),B(4),DEL(4),SIG(4)  ! Changed dimension statement so that A,B,DEL, and SIG have 4 entries. SIGR is just a scalar. 5/3/2011
C
      REAL SIGR, AIRSCAT
      DATA A/29.06,26.63,43.9,27.92/
      DATA B/7.7,5.07,6.4,5.6/
      DATA DEL/.0305,.054,.0805,0.032/
	  
	  
C
      !DO 1135 I=1,NSOL                 Sending SIGRs one wavelength at a time so do not need these statements 5/3/2011
      !   AL2=ALAMBDA(I)*ALAMBDA(I)
	  
         AL4=AL2*AL2
		 
cRKK- THE EXPRESSION FOR SIG(J) IS FROM VARDAVAS & CARVER, 1984, 
cRKK- PLANETARY SPACE SCIENCE, VOL 32, P-1307.
C
         DO 1140 J=1,4
            PAREN=(1.E-5*A(J)*(1.+1.E-3*B(J)/AL2))**2
            SMDEL=(6.+3.*DEL(J))/(6.-7.*DEL(J))
            SIG(J)=4.577E-21*SMDEL*PAREN/AL4
 1140    CONTINUE
C	
	  ! SIG(1) is Nitrogen. SIG(2) is oxygen . SIG(3) is carbon dioxide. SIG(4) IS ARGON. These four components make up the rayleigh scattering 
      ! component due to air. The fifth term is the rayleigh scattering component due to water, which is assumed to be 73% that of air
	  !(Selsis et al.(2007) Icarus, 191, 453)
	  ! Unfortunately, empirical constants A,B, and DEL for methane don't exist yet so it is not included. It is such a tiny component,however
	  ! That it shouldn't affect the results much at all. 
         AIRSCAT = FARNC*SIG(4) +FCO2*SIG(3) + FO2NC*SIG(2) +
     &     FN2*SIG(1)
	     SIGR = AIRSCAT + 0.734*AIRSCAT*Fwater
 1135 CONTINUE 
C
      RETURN
      END  
