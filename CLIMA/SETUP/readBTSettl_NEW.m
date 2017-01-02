clear all ; clc 

s=[3 3 4 4 5];
ss= [0 5 0 5 0];

for T = 26:2:72
       for j = 1:5
             
file=sprintf('lte0%d-%d.%d-0.0a+0.0.BT-Settl.7', T, s(j), ss(j));


A = load(file);
A= sortrows(A);  % sorts first column (wavelength) in ascending order with associated fluxes

WL = A(:,1)/10000;  % Wavelength in angstroms converted to microns

F= 10000.*(10.^(A(:,2)-8.d0)); % converting from log10(f_lambda) to ergs/cm^2/sec/micron
grids = [.2376 .2750 .2850 .3071 .3292 .3412 .3900 .4500 .5400 .5495 .5666 .6050 .6250 .6667 .6910 .7520 .7840 .8420 .8910 .9620 1.0360 1.0700 1.1300 1.2030 1.3070 1.4310 1.5650 1.6880 1.8620 2.0200 2.2030 2.4810 2.6600 2.920 3.2390 3.5770 4.0100 4.172 4.545];


indices = find(WL>=grids(1) & WL<=grids(length(grids)));  % find indices of all wavelengths/fluxes within our grid
WL=WL(indices);
F =F(indices);

% Plots function
%subplot(1,2,1)
%plot(WL, F) 
%pause

%interpolate the whole thing
DWL=.0001;   % each wavelength point separated by just one angstrom
first = WL(1);  % first wavelength point in our grid
last = WL(length(WL)); % last wavelength point in our grid
WLi= [first:DWL:last];
F=interp1(WL,F,WLi);

%Plots interpolated function
%subplot(1,2,2)
%plot(WLi,F)
%xlabel('microns')
%ylabel('ergs/cm^2/sec/micron')
%pause

% Getting delta lamba and F vectors needed to get dimensions of 
% infinitesimal fluxes that need to be added together for each bin
for n = 1:length(WLi)-1
Flux(n) = 0.5*(F(n)+ F(n+1)); % midpoint flux approximation
end


Flux = Flux*DWL;   % infinitesimal fluxes now in ergs/cm^2/sec
Flux(length(Flux)+1)= 0 ;  % because need to increment Flux vector by 1 to equal number of points 

for n = 1: length(grids)-1  %length(grid)-1 
index=(find(WLi>=grids(n) & WLi<=grids(n+1)));  % find indices of fluxes within a grid square
Bin(n) = sum(Flux(index)); % summing all the fluxes within given bin
%empty=find(Bin==isempty(Bin)) 
end

FatEarth = sum(Bin);

Normalized = Bin.*(1/(sum(Bin)/1.36e6));
% sum(Normalized)

diary on
file
fprintf(' %e\n', Normalized)
%fprintf(' %i', T, s(j), ss(j))
fprintf('\n')
diary off
% For plotting at midpoint
%for n =1:length(grids)-1
%    grids2(n) = 0.5*(grids(n)+grids(n+1));
%end

%plot(grids2,Normalized)
%xlabel('microns')
%ylabel('ergs/cm^2/sec')
         end
 end

%%
clear all ; clc 
A=load('3400logg5Met0.txt');
A= sortrows(A);  % sorts first column (wavelength) in ascending order with associated fluxes

WL = A(:,1)/10000;  % Wavelength in angstroms converted to microns

F= 10000.*(10.^(A(:,2)-8.d0)); % converting from log10(f_lambda) to ergs/cm^2/sec/micron

grids = [.2376 .2750 .2850 .3071 .3292 .3412 .3900 .4500 .5400 .5495 .5666 .6050 .6250 .6667 .6910 .7520 .7840 .8420 .8910 .9620 1.0360 1.0700 1.1300 1.2030 1.3070 1.4310 1.5650 1.6880 1.8620 2.0200 2.2030 2.4810 2.6600 2.920 3.2390 3.5770 4.0100 4.172 4.545];


indices = find(WL>=grids(1) & WL<=grids(length(grids)));  % find indices of all wavelengths/fluxes within our grid
WL=WL(indices);
F =F(indices);
plot(WL, F)
hold on

%%

A=load('3400logg5Met-4.txt');
A= sortrows(A);  % sorts first column (wavelength) in ascending order with associated fluxes

WL = A(:,1)/10000;  % Wavelength in angstroms converted to microns

F= 10000.*(10.^(A(:,2)-8.d0)); % converting from log10(f_lambda) to ergs/cm^2/sec/micron

grids = [.2376 .2750 .2850 .3071 .3292 .3412 .3900 .4500 .5400 .5495 .5666 .6050 .6250 .6667 .6910 .7520 .7840 .8420 .8910 .9620 1.0360 1.0700 1.1300 1.2030 1.3070 1.4310 1.5650 1.6880 1.8620 2.0200 2.2030 2.4810 2.6600 2.920 3.2390 3.5770 4.0100 4.172 4.545];


indices = find(WL>=grids(1) & WL<=grids(length(grids)));  % find indices of all wavelengths/fluxes within our grid
WL=WL(indices);
F =F(indices);
plot(WL, F,'r')

%% 
clear all;clc
grids2=[2.563000e-01 2.800000e-01 2.960500e-01 3.181500e-01 3.352000e-01 3.656000e-01 4.200000e-01 4.950000e-01 5.447500e-01 5.580500e-01 5.858000e-01 6.150000e-01 6.458500e-01 6.788500e-01 7.215000e-01 7.680000e-01 8.130000e-01 8.665000e-01 9.265000e-01 9.990000e-01 1.053000e+00 1.100000e+00 1.166500e+00 1.255000e+00 1.369000e+00 1.498000e+00 1.626500e+00 1.775000e+00 1.941000e+00 2.111500e+00 2.342000e+00 2.570500e+00 2.790000e+00 3.079500e+00 3.408000e+00 3.793500e+00 4.091000e+00 4.358500e+00]

%2800logg4.5Met0.5 (high Metallicity)

HighMet2800=[7.364316e-02 6.749493e-02 2.640066e+00 6.438083e+01 7.549063e+01 3.201181e+02 1.634924e+03 4.795518e+03 6.914937e+02 1.052291e+03 2.577947e+03 1.718818e+03 4.788281e+03 2.397839e+03 1.597098e+04 1.178389e+04 3.853944e+04 2.687344e+04 6.656962e+04 8.253973e+04 4.119029e+04 7.364068e+04 8.650522e+04 1.217559e+05 1.075981e+05 9.390919e+04 9.225593e+04 9.722589e+04 6.188772e+04 6.467791e+04 8.024212e+04 2.952551e+04 3.357818e+04 3.133111e+04 2.859355e+04 3.025333e+04 8.749188e+03 1.468317e+04]

%2800logg4.5Met-4 (low Metallicity)

LowMet2800=[8.379089e-01 4.907273e-02 6.006292e-01 1.761347e+02 1.175437e+02 2.978984e+02 3.767707e+03 3.322496e+04 4.888461e+03 1.258441e+04 3.505887e+04 2.343035e+04 5.407834e+04 2.535147e+04 9.431371e+04 5.622176e+04 1.099124e+05 9.283264e+04 1.324599e+05 1.113799e+05 3.997596e+04 6.325428e+04 7.096127e+04 8.528426e+04 7.564808e+04 5.091766e+04 2.922210e+04 2.639200e+04 1.625400e+04 1.393494e+04 1.570388e+04 9.482192e+03 1.424503e+04 1.695024e+04 1.516966e+04 1.448857e+04 4.212624e+03 7.805287e+03] 

%3800logg4.5Met0.5 (high Metallicity)

HighMet3800=[9.626573e+00 9.088990e+00 6.478489e+01 3.978287e+02 4.555437e+02 2.413908e+03 1.133756e+04 3.589075e+04 4.483690e+03 8.575352e+03 2.256306e+04 1.149089e+04 2.462724e+04 1.407443e+04 5.338286e+04 3.333690e+04 6.517030e+04 5.136523e+04 7.721256e+04 7.613387e+04 3.392053e+04 5.693625e+04 6.512572e+04 8.604621e+04 9.184156e+04 8.987019e+04 7.632114e+04 8.705548e+04 5.726479e+04 5.218667e+04 5.398688e+04 2.171347e+04 2.627096e+04 2.300495e+04 1.770585e+04 1.579060e+04 4.452092e+03 7.511157e+03]

%3800logg4.5Met-4 (low Metallicity)
LowMet3800=[1.619540e+02 6.127073e+01 3.347856e+02 1.382215e+03 9.905589e+02 5.918369e+03 1.975076e+04 5.892420e+04 8.002075e+03 1.591732e+04 3.939569e+04 2.200188e+04 4.814221e+04 2.911760e+04 7.519963e+04 4.015840e+04 7.218808e+04 5.956509e+04 8.307954e+04 8.100501e+04 3.511590e+04 5.854126e+04 6.583149e+04 8.343693e+04 8.381722e+04 7.108952e+04 4.976755e+04 5.289813e+04 3.605748e+04 3.185915e+04 3.512146e+04 1.728096e+04 2.024039e+04 1.891710e+04 1.469924e+04 1.321578e+04 3.791110e+03 7.022724e+03]

plot(grids2, HighMet2800,'r')
%legend('2800(High Met)')
pause
hold on
plot(grids2, LowMet2800,'g')
%legend('2800(High Met)', '2800(Low Met)')
%pause
%plot(grids2, HighMet3800,'bl')
%legend('2800(High Met)', '2800(Low Met)', '3800(High Met)')
%pause
%plot(grids2, LowMet3800,'k')
xlabel('microns')
ylabel('ergs/cm^2/sec')
%legend('2800(High Met)', '2800(Low Met)', '3800(High Met)', '3800(Low Met)')


%%

clear all ; clc 


A=load('2800logg4.5Met0.5.txt');
A= sortrows(A);  % sorts first column (wavelength) in ascending order with associated fluxes

WL = A(:,1)/10000;  % Wavelength in angstroms converted to microns

F= 10000.*(10.^(A(:,2)-8.d0)); % converting from log10(f_lambda) to ergs/cm^2/sec/micron
grids = [.2376 .2750 .2850 .3071 .3292 .3412 .3900 .4500 .5400 .5495 .5666 .6050 .6250 .6667 .6910 .7520 .7840 .8420 .8910 .9620 1.0360 1.0700 1.1300 1.2030 1.3070 1.4310 1.5650 1.6880 1.8620 2.0200 2.2030 2.4810 2.6600 2.920 3.2390 3.5770 4.0100 4.172 4.545];


indices = find(WL>=grids(1) & WL<=grids(length(grids)));  % find indices of all wavelengths/fluxes within our grid
WL=WL(indices);
F =F(indices);
plot(WL, F) 
hold on



A=load('2800logg4.5Met-4.txt');
A= sortrows(A);  % sorts first column (wavelength) in ascending order with associated fluxes

WL = A(:,1)/10000;  % Wavelength in angstroms converted to microns

F= 10000.*(10.^(A(:,2)-8.d0)); % converting from log10(f_lambda) to ergs/cm^2/sec/micron
grids = [.2376 .2750 .2850 .3071 .3292 .3412 .3900 .4500 .5400 .5495 .5666 .6050 .6250 .6667 .6910 .7520 .7840 .8420 .8910 .9620 1.0360 1.0700 1.1300 1.2030 1.3070 1.4310 1.5650 1.6880 1.8620 2.0200 2.2030 2.4810 2.6600 2.920 3.2390 3.5770 4.0100 4.172 4.545];


indices = find(WL>=grids(1) & WL<=grids(length(grids)));  % find indices of all wavelengths/fluxes within our grid
WL=WL(indices);
F =F(indices);
%subplot(1,2,1)
plot(WL, F) 
xlabel('microns')
ylabel('ergs/cm^2/sec/micron')