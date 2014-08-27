clear  all;clc
% IR works in wavenumber space so spectra units W/m^2/cm^-1
format short
B = input('What is the file name (plus extension)?:','s')
A = load(B);
T = input('What is the surface temperature (in K)?:')

grids =[0.,40., 100., 160., 220., 280., 330., 380., 440., 495.,...
    545., 617., 667., 720., 800., 875., 940., 1000., 1065.,...
       1108., 1200., 1275., 1350., 1450., 1550., 1650., 1750., 1850.,...
       1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760.,...
       4030., 4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315.,...
       8850., 9350., 9650., 10400., 11220., 11870., 12790., 13300.,...
       14470., 15000.]; % Wavenumber (per centimeters)
   
 %for n =1:length(AV)-1
 %   AV(n) = 0.5*(AV(n)+AV(n+1));%bins need to be taken at midpoint (in microns)
 %end
   
AV=A(:,1) %bins need to be taken at midpoint (Wavenumber in (per centimeters))

Flux=A(:,2)/1000;   %Fluxes in W/m^2

Flux2(1)=Flux(1);
for i = 1:length(Flux)-1
    Flux2(i+1)= Flux(i+1)-Flux(i);% Fluxes given as cumulative sums in ir.f Needed to get individual fluxes!
    %Flux2(i)= Flux(i+1)-Flux(i)
end
%Flux2(length(Flux2)+1)=0.
%plot(AV, Flux2)
%xlabel('Wavenumber(per centimeters)') 
%ylabel('OLR(W/m^2)')
%pause


for j =1:length(AV)
    DWAV(j) = grids(j+1)-grids(j);  % Delta wavenumbers (in per centimeters)
end



%Fluxes in W*centimeter/m^2
Fluxbin=Flux2./DWAV;
TotalF1= sum(Flux2(1:29))
TotalF1= sum(Flux2)



% Blackbody emission curve for given surface temperature. Corresponds to FTIR with
% no atmospheric absorption

% in cgs units
c= 3E10 ; % cm/sec
h=6.63E-27;
k=1.38E-16;

v= [2000:-10:.1]; % wavenumber (inverse centimeters) (5 microns to 100,000 microns) This is why 
% wavenumber version is good for the IR.


NUM = 2*h*c^2.*v.^3;
DEN = exp((h*c.*v)/(k*T))-1 ;

PLANCK = .001*pi*NUM./DEN;  %divide by 1000 to convert from ergs/cm^2/sec/wavenumber/steradian to W/m^2/wavenumber/steradian
% Multiplied by pi to cancel steradians.


plot(AV, Fluxbin)
hold on
plot(v, PLANCK,'r')
axis([0 2000 0 0.4])
xlabel('Wavenumber(per cm)')
ylabel('OLR(W/m^2/cm^-1)')
%hold on
%%
clear  all;clc
% IR works in wavenumber space so spectra units W/m^2/cm^-1
format shortE
B = input('What is the file name (plus extension)?:','s')
A = load(B);
T = input('What is the surface temperature (in K)?:')

grids =[0.,40., 100., 160., 220., 280., 330., 380., 440., 495.,...
    545., 617., 667., 720., 800., 875., 940., 1000., 1065.,...
       1108., 1200., 1275., 1350., 1450., 1550., 1650., 1750., 1850.,...
       1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760.,...
       4030., 4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315.,...
       8850., 9350., 9650., 10400., 11220., 11870., 12790., 13300.,...
       14470., 15000.]; % Wavenumber (per centimeters)
   
AV=A(:,1) %bins need to be taken at midpoint (Wavenumber in (per centimeters))

Flux=A(:,2)/1000;   %Fluxes in W/m^2 
%Flux=B(:,2)/1000

Flux2(1)=Flux(1);
for i = 1:length(Flux)-1
    Flux2(i+1)= Flux(i+1)-Flux(i)% Fluxes given as cumulative sums in ir.f Needed to get individual fluxes!
    %Flux2(i)= Flux(i+1)-Flux(i)
end
%Flux2(length(Flux2)+1)=0.
%plot(AV, Flux2)
%xlabel('Wavenumber(per centimeters)') 
%ylabel('OLR(W/m^2)')
%pause

for j =1:length(AV)
    DWAV(j) = grids(j+1)-grids(j);  % Delta wavenumbers(in per centimeters)
end



%Fluxes in W*centimeter/m^2
Fluxbin=Flux2./DWAV;
TotalF2= sum(Flux2(1:29))
TotalF2= sum(Flux2)

% Blackbody emission curve for given surface temperature. Corresponds to FTIR with
% no atmospheric absorption

% in cgs units
c= 3E10 ; % cm/sec
h=6.63E-27;
k=1.38E-16;

v= [2000:-10:.1]; % wavenumber (inverse centimeters) (5 microns to 100,000 microns) This is why 
% wavenumber version is good for the IR.


NUM = 2*h*c^2.*v.^3;
DEN = exp((h*c.*v)/(k*T))-1 ;

PLANCK = .001*pi*NUM./DEN;  %divide by 1000 to convert from ergs/cm^2/sec/wavenumber/steradian to W/m^2/wavenumber/steradian
% Multiplied by pi to cancel steradians.


plot(AV, Fluxbin)
hold on
plot(v, PLANCK,'r')
axis([0 2000 0 0.4])
xlabel('Wavenumber(per cm)')
ylabel('OLR(W/m^2/cm^-1)')