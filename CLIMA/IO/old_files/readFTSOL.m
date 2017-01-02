clear all;clc
% Solar works in wavlength space because as wavelength increases, fluxes
% increase.

format shortE
B = input('What is the file name (plus extension)?:','s')
%A = load('FTSOL10Gliese.dat');

A = load(B)
SOLCON = input('What is the SOLCON value (S/So) for your planet?:')
Di= 0.5 % diurnal averaging
% Assuming 60 degree zenith angle...
 
%solar grid points (in microns)
grids = [.2376 .2750 .2850 .3071 .3292 .3412 .3900 .4500 .5400 .5495 .5666 .6050 .6250 .6667 .6910 .7520 .7840 .8420 .8910 .9620 1.0360 1.0700 1.1300 1.2030 1.3070 1.4310 1.5650 1.6880 1.8620 2.0200 2.2030 2.4810 2.6600 2.920 3.2390 3.5770 4.0100 4.172 4.545];

for n =1:length(grids)-1
    grids2(n) = 0.5*(grids(n)+grids(n+1));%bins need to be taken at midpoint (in microns)
end

 AV=grids2; % wavelength midpoints in microns


Fluxbinstar = SOLCON*A(:,2)/1000  % Binned Stellar fluxes in W/m^2. Re-normalized to flux of current star at planet's distance(not 1360 W/m^2)
Fluxbin=SOLCON*Di*(A(:,1)/1000);   %Binned cumulative net incident fluxes in W/m^2 (takes absorption coefficients and rayleigh scattering both into account)
%Renormalized to current star flux at planet's distance and taking diurnal averaging into
    %account.

Flux2bin(1)=Fluxbin(1)
for i = 1:length(Fluxbin)-1
    Flux2bin(i+1)= (Fluxbin(i+1)-Fluxbin(i))% individual net incident fluxes received by planet in every interval calculated from cumulative sum (in W/m^2).
    
end

 %DWAV(1) = AV(1)
 
 for j =1:length(AV)
     DWAV(j) = abs(grids(j+1)-grids(j)) % Delta wavelengths (in microns)
 end
 
plot(AV, Fluxbinstar)  % plots stellar binned spectrum at midpoints
title('Binned Stellar Spectrum')
xlabel('Wavelength(microns)')
ylabel('Flux(W/m^2)')
%pause

%globally averaged binned stellar flux is Fluxbinstar/4 
GFluxbinstar= Fluxbinstar/4

plot(AV, GFluxbinstar)
title('Binned Globally-Averaged Stellar Spectrum')
xlabel('Wavelength(microns)')
ylabel('Flux(W/m^2)')
%pause

%Converting binned fluxes into spectrum fluxes (in W/m^2/micron)
Flux=Flux2bin./DWAV;  % net incident flux spectrum
GFluxstar = GFluxbinstar./DWAV'; % globally averaged stellar flux

% Total fluxes (add all bins together)
Totalfluxstar = sum(Fluxbinstar); %  total stellar flux at planetary distance before reaching atmosphere(in W/m^2)
TotalGfluxstar=sum(GFluxbinstar); %  total globally-averaged stellar flux (in W/m^2)
TotalFluxplanet= sum(Flux2bin); %  total net incident fluxes received at surface(in W/m^2)
fprintf('Total stellar flux at planet distance is %f W/m^2 \n',Totalfluxstar)
fprintf('Total globally-averaged stellar flux is %f W/m^2 \n',TotalGfluxstar)
fprintf('Total net incident flux is %f W/m^2:',TotalFluxplanet)


%plot (AV, GFluxstar, AV, (GFluxstar-Flux')) % plots globally-averaged star spectrum vs. globally averaged star-planet
%xlabel('Wavelength(microns)')
%ylabel('ISW(W/m^2/micron)')
%legend('Star', 'Star -Planet')
%pause

%To smooth out spectra by averaging and weighting the neighboring points
%alpha=0.5
%GFluxstar(2:38-1) = alpha*GFluxstar(2:38-1) + (1-alpha)*0.5*(GFluxstar(1:38-2)+GFluxstar(3:38))
%Flux(2:38-1) = alpha*Flux(2:38-1) + (1-alpha)*0.5*(Flux(1:38-2)+Flux(3:38))

plot(AV,GFluxstar,AV, Flux)  % plots globally-averaged spectrum vs. planet spectrum only
legend('Globally-Averaged Stellar Spectrum','Net Incident Flux Spectrum')
xlabel('Wavelength(microns)')
ylabel('Net Incoming Solar Flux(W/m^2/micron)')
hold on


%format longG
C(:,1) = (1:38);
C(:,2) = AV;
C(:,3) = Flux;
diary on
        disp('  Interval     *Wavelength( microns )     Net incident flux(W/m^2/micron)')
for j = 1:38
        fprintf('%4.0f              %5.4f                   %g',C(j,1), C(j,2), C(j,3))
        fprintf('\n')
end
disp('*Note that these solar flux data are plotted against the midpoints of their intervals')
diary off


