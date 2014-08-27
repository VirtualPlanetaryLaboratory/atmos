% Create shaded region in the HZ plot
% Ravi kumar Kopparapu Oct 27 2012
% 
% Read the data file for observed planets 'habplanets.dat'
% Read flux data file for FGKM stars 'HZ_fluxes.dat'
clc;clear
[name mass teff mp aa lum flux ref]=textread('habplanets.dat','%s%f%f%f%f%f%f%s',9,'headerlines',1);
[steff,rv,rg,mg,maxg,em,AU2]= textread('files/HZ_fluxes.dat','%f%f%f%f%f%f%f',24,'headerlines',2);


% Plot 'Recent Venus' limit
h(1)=plot(rv,steff,'r');
hold on

% Plot 'Runaway Greenhouse' limit
h(2)=plot(rg,steff,'k');

% Plot 'Early Mars' limit
h(3)=plot(em,steff,'b');
axis([0 2.2 2600 7200])

% Store the values of moist greenhouse, maximum greenhouse & Teff
% Use fill command to shade the HZ (between moist greenhouse & maximum
% greenhouse)
maxg1 = [ 0.300 0.289 0.280 0.271 0.255 0.236  0.215  0.8423 0.8474 0.8530 0.8591 0.8660 0.8736 0.8821 0.8915 0.9019 0.9133 0.9256 0.9389 0.9530 0.9679 0.9833 0.9993 1.0156 1.0320 1.0483 1.0642 1.0795 1.0938 1.1069 1.1184 0.42];
steff1 = [5000  4800  4600  4400  4000  3400   2600   2600  2800  3000  3200  3400  3600  3800  4000  4200  4400  4600  4800  5000  5200  5400  5600  5800  6000  6200  6400  6600  6800  7000   7200  7200];
fill(maxg1,steff1,'g')

% Plot observed planets on the HZ plot
h(4)=plot(flux(1),teff(1),'^y','MarkerSize',10);
h(5)=plot(flux(2),teff(2),'dm','MarkerSize',10);
h(6)=plot(flux(3),teff(3),'+r','MarkerSize',10);
h(7)=plot(flux(4),teff(4),'*k','MarkerSize',10);
h(8)=plot(flux(5),teff(5),'ob','MarkerSize',10);
h(9)=plot(flux(6),teff(6),'sb','MarkerSize',10);
h(10)=plot(flux(7),teff(7),'pb','MarkerSize',10);
h(11)=plot(flux(8),teff(8),'pc','MarkerSize',10);

set(gca,'XDir','reverse')
Legend(h(4:11),'Gl 581d','Gl 581g','Gl 667Cc','Kepler 22b','Earth','HD 40307g','Tau Ceti e (?)','Tau Ceti f (?)')
xlabel('Effective flux incident on the planet (S/S_{o})')
ylabel('Stellar effective temperature T_{eff} (K)')
line([0.842 0.842],[2600 7200])
line([0.42 0.42],[2600 7200])


hold off



