clear all;clc

[J P ALT T CONVEC DT TOLD FH2O FSAVE TCOOL THEAT] = textread('clima_allout1bar.tab','%f %f %f %f %f %f %f %f %f %f %f',101,'headerlines',939)
n=101  %Number of layers down starting from top of atmosphere
% Plots T vs Altitude
plot(T(1:n), ALT(1:n))
xlabel('Temperature(K)')
ylabel('Altitude(km)')


hold on

[J2 P2 ALT2 T2 CONVEC2 DT2 TOLD2 FH2O2 FSAVE2 TCOOL2 THEAT2] = textread('clima_allout1e-2met.tab','%f %f %f %f %f %f %f %f %f %f %f',101,'headerlines',739)

% Plots T vs Altitude
plot(T2(1:n), ALT2(1:n))
xlabel('Temperature(K)')
ylabel('Altitude(km)')