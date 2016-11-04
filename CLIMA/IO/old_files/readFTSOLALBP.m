clear all;clc
%This reads in FUPSOL in first column, FDNSOl in second. It uses them to
% calculate the individual planetary albedos components for every bin and the planetary albedo

format shortE
B = input('What is the file name (plus extension)?:','s')
SOLCON = input('What is the SOLCON value (S/So) for your planet?:')
Di= 0.5 % diurnal averaging
% Assuming 60 degree zenith angle...

A = load(B)
 
%solar grid points (in microns)
grids = [.2376 .2750 .2850 .3071 .3292 .3412 .3900 .4500 .5400 .5495 .5666 .6050 .6250 .6667 .6910 .7520 .7840 .8420 .8910 .9620 1.0360 1.0700 1.1300 1.2030 1.3070 1.4310 1.5650 1.6880 1.8620 2.0200 2.2030 2.4810 2.6600 2.920 3.2390 3.5770 4.0100 4.172 4.545];

for n =1:length(grids)-1
    grids2(n) = 0.5*(grids(n)+grids(n+1));%bins need to be taken at midpoint (in microns)
end

AV=grids2; % wavelength midpoints in microns

% The cumulative upward and downward fluxes, respectively.
FUPSOLBIN = SOLCON*Di*A(:,1)
FDNSOLBIN = SOLCON*Di*A(:,2)


FUPSOL(1)=FUPSOLBIN(1)
FDNSOL(1)=FDNSOLBIN(1)

for i = 1:length(FUPSOLBIN)-1
    FUPSOL(i+1)= (FUPSOLBIN(i+1)-FUPSOLBIN(i))% individual upward fluxes leaving planet in every interval calculated from cumulative sum (in W/m^2).
    FDNSOL(i+1)= (FDNSOLBIN(i+1)-FDNSOLBIN(i))% individual downward fluxes entering planet in every interval calculated from cumulative sum (in W/m^2).
end

ALBP = FUPSOL./FDNSOL;
Planetary_Albedo=sum(FUPSOL)/sum(FDNSOL);
Downwelling_flux= sum(FDNSOL)/1000; % in W/m^2
Upwelling_flux=sum(FUPSOL)/1000; % in W/m^2
fprintf('The Planetary Albedo is %f \n',Planetary_Albedo)
fprintf('The downwelling flux at the top is is %f W/m^2 \n',Downwelling_flux)
fprintf('The upwelling flux at the top is %f W/m^2 \n',Upwelling_flux)



plot(AV, ALBP)
title('Planetary Albedo vs. Wavelength')
xlabel('Wavelength(microns)')
ylabel('Planetary Albedo')

