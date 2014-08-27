 % This program reads inner edge output files and generates various plots.
clear all;clc

[TG0 SEFF PALB FH2O FTIR FTSO]=textread('waterloss_IHZ.dat', '%f %f %f %f %f %f', 39, 'headerlines',4)


TG02 = [220:10:1600];
PALB2 = spline(TG0,PALB,TG02);

%plot(TG02,PALB2,'rx-')
%hold on
%plot(TG0,PALB,'bo')
%title('RUNAWAY GREENHOUSE ALBEDO VS. SURFACE TEMP')
%xlabel('Surface Temperature (K)')
%ylabel('Planetary Albedo')
%axis([200 620 0.16 .36])
%hold off
%pause


%FTIR2 = interp1(TG0,FTIR,TG02);
%FTIR3 = spline(TG0,FTIR,TG02);


%TG02  = TG0;
%FTIR3 = supsmu(TG02,FTIR);

%plot(TG02,FTIR3,'go-')
%hold on
%plot(TG0,FTIR,'rx-')
%hold on
%plot(TG0,FTSO,'bx-')
%title('Net OUTGOING THERMAL AND INCOMING SOLAR FLUXES')
%xlabel('Surface Temperature (K)')
%ylabel('Flux (W/m^2)')

%axis([200 1400 200 320])
%hold off
%pause

  SEFF2 = splinefit(TG0,SEFF,28);
  SEFF3=PPVAL(SEFF2,TG02);
 plot(TG02,SEFF3,'rx-')
 hold on
 plot(TG0,SEFF,'bo');
 hold off
 
% xlabel('Surface Temperature (K)')
% ylabel('S_e_f_f')
axis([200 1500 0.5 1.8])
% 
% pause
% semilogy(SEFF, FH2O,'rx-')
% xlabel('S_e_f_f')
% ylabel('Volumetric H_2O Mixing Ratio')
% axis([0.5 1.3 1e-6 1])
% pause