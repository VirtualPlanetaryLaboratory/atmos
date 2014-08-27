
c     SUBROUTINE INTERP(Temp1,p1,xkappa,kappa)
      SUBROUTINE INTERP(Temp1,p1,kappa)

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

C new common block, von Paris, 21/04/2006
	include 'CLIMA/INCLUDE/header.inc'
            COMMON/IRDATA/WEIGHT(8,3),xkappa(8,12,55,8,3), 
     & CIA(7,55), CPRW(ND,55)  
c-rr !3/23/11 put CIA matrix in IRDATA
 
      INTEGER LINE,i,j,lam,Tindex,pindex,Tp1,pp1,Temp
      REAL p,Temp1
c     REAL xkappa(8,12,55,8,3)
      REAL kappapp(55,8,3),kappapm(55,8,3)
      REAL kappa(55,8,3)
 
c-mm  These two loops will linearly interpolate kappa values to solve
c-mm  for any T and p.  kappapm interpolates xkappa for the two values
c-mm  having a pressure lower than p.  kappapp interpolates xkappa for
c-mm  the two values having a pressure higher than p.  kappa then
c-mm  interpolates temperature.
      do 9400 i=1,55
         do 9401 j=1,8
            do 9402 k=1,3
               kappa(i,j,k)=0.0
 9402       continue
 9401    continue
 9400 continue
        Temp=int(Temp1)
        Tindex=Temp/50
        p = p1*1000
        pindex=int(log10(p*1.0e7))
        Tp1=Tindex+1
        pp1=pindex+1
        do 9998 lam=9,48
           do 9995 i=1,8
              kappapm(lam,i,1)=((mod(Temp,50))/(50.))*
     2             xkappa(Tp1,pindex,lam,i,1)+
     3             (1.-((mod(Temp,50))/(50.)))*
     4             xkappa(Tindex,pindex,lam,i,1)
              kappapp(lam,i,1)=((mod(Temp,50))/(50.))*
     2             xkappa(Tp1,pp1,lam,i,1)+
     3             (1.-((mod(Temp,50))/(50.)))*
     4             xkappa(Tindex,pp1,lam,i,1)
              if (p.le.(1.0e0)) then
                 kappa(lam,i,1)=(1.-(abs(log10(p)-int(log10(p)))))*
     2             kappapp(lam,i,1)+abs(log10(p)-
     3             int(log10(p)))*kappapm(lam,i,1)
              else
                 kappa(lam,i,1)=(log10(p)-int(log10(p)))*
     2             kappapp(lam,i,1)+(1.-(log10(p)-int(log10(p))))*
     3             kappapm(lam,i,1)
              endif
 9995      continue
 9998   continue
  
        do 9949 lam=1,55
           do 9948 i=1,8
              kappapm(lam,i,2)=((mod(Temp,50))/50.)*
     2             xkappa(Tp1,pindex,lam,i,2)+
     3             (1.-((mod(Temp,50))/50.))*
     4             xkappa(Tindex,pindex,lam,i,2)
              kappapp(lam,i,2)=((mod(Temp,50))/50.)*
     2             xkappa(Tp1,pp1,lam,i,2)+
     3             (1.-((mod(Temp,50))/50.))*
     4             xkappa(Tindex,pp1,lam,i,2)
              if (p.le.(1.0e0)) then
                 kappa(lam,i,2)=(1-(abs(log10(p)-int(log10(p)))))*
     2             kappapp(lam,i,2)+abs(log10(p)-
     3             int(log10(p)))*kappapm(lam,i,2)
              else
                 kappa(lam,i,2)=(log10(p)-int(log10(p)))*
     2             kappapp(lam,i,2)+(1.-(log10(p)-int(log10(p))))*
     3             kappapm(lam,i,2)
              endif
 9948         continue
 9949   continue
 
        do 9607 lam=1,38
	if (Temp.le.150) then
		do 9620 i=1,6
		kappapm(lam,i,3) = xkappa(3,pindex,lam,i,3)
		kappapp(lam,i,3) = xkappa(3,pp1,lam,i,3)
		kappa(lam,i,3) = (1.-(abs(log10(p)-int(log10(p)))))*
     2  kappapp(lam,i,3)+abs(log10(p)-int(log10(p)))*kappapm(lam,i,3)
 9620  continue
        elseif (Temp.le.300) then
                do 9609 i=1,6
                        kappapm(lam,i,3)=((mod(Temp,150))/150.)*
     2  xkappa(6,pindex,lam,i,3)+(1.-((mod(Temp,150))/150.))*
     3  xkappa(3,pindex,lam,i,3)
                        kappapp(lam,i,3)=((mod(Temp,150))/150.)*
     2  xkappa(6,pp1,lam,i,3)+(1.-((mod(Temp,150))/150.))*
     3  xkappa(3,pp1,lam,i,3)
                        kappa(lam,i,3)=(1.-(abs(log10(p)-
     2  int(log10(p)))))*kappapp(lam,i,3)+
     3  abs(log10(p)-int(log10(p)))*kappapm(lam,i,3)
 9609   continue
        else
                do 9608 i=1,6
                        kappapm(lam,i,3)=((mod(Temp,300))/300.)*
     2  xkappa(8,pindex,lam,i,3)+(1.-((mod(Temp,300))/300.))*
     3  xkappa(6,pindex,lam,i,3)
                        kappapp(lam,i,3)=((mod(Temp,300))/300.)*
     2  xkappa(8,pp1,lam,i,3)+(1.-((mod(Temp,300))/300.))*
     3  xkappa(6,pp1,lam,i,3)
                        kappa(lam,i,3)=(1-(abs(log10(p)-
     2  int(log10(p)))))*kappapp(lam,i,3)+
     3  abs(log10(p)-int(log10(p)))*kappapm(lam,i,3)
 9608           continue
        endif
 9607   continue

        RETURN
        END
