% Checking Temperature profiles...
clear all;clc

[ALT PRESS TEMP H2O_mix_ratio] = textread('profile.dat', '%f %f %f %f', 101, 'headerlines',1)
plot(TEMP,PRESS)
set(gca, 'YDir', 'reverse')

axis([160 250 0 0.5])
xlabel('Temperature (K)')
ylabel('Pressure (bars)')

%%
% Upward Fluxes
clear all;clc

[J PF ALT FTOTAL FTIR FDNIR FUPIR FTSOL FDNSOL FUPSOL DIVF] = textread('clima_allout.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 132)

FUPSOL = FUPSOL/1000 % convert to W/m^2
FUPIR  = FUPIR/1000 % convert to W/m^2

%Ours
%subplot 121
semilogy(FUPSOL, PF)
set(gca, 'YDir', 'reverse')
%axis([15 30 .1E-5 1])
xlabel('Upward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
%hold on
pause

%subplot 122
semilogy(FUPIR, PF)
set(gca, 'YDir', 'reverse')
%axis([110 230 .1E-5 1])
xlabel('Upward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on


%%
% Downward Fluxes
clear all;clc

[J PF ALT FTOTAL FTIR FDNIR FUPIR FTSOL FDNSOL FUPSOL DIVF] = textread('clima_allout.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 132)

FDNSOL = FDNSOL/1000 % convert to W/m^2
FDNIR  = FDNIR/1000 % convert to W/m^2

%Ours
%subplot 121
semilogy(FDNSOL, PF)
set(gca, 'YDir', 'reverse')
%axis([80 115 .1E-5 1])
xlabel('Downward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
%hold on
pause

%subplot 122
semilogy(FDNIR, PF)
set(gca, 'YDir', 'reverse')
axis([-5 200 .1E-5 2])
xlabel('Downward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

%% Net Fluxes

[J PF ALT FTOTAL FTIR FDNIR FUPIR FTSOL FDNSOL FUPSOL DIVF] = textread('clima_allout.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 132)

FTSOL = FTSOL/1000 % convert to W/m^2
FTIR  = abs(FTIR/1000) % convert to W/m^2

%Ours
%subplot 121
semilogy(FTSOL, PF)
set(gca, 'YDir', 'reverse')
%axis([65 85 .1E-5 1])
xlabel('Net Incoming Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
%hold on
pause

%subplot 122
semilogy(FTIR, PF)
set(gca, 'YDir', 'reverse')
%axis([20 130 .1E-5 1])
xlabel('Net Outgoing IR Flux (W/m^2)')
ylabel('Pressure (bars)')


%pause
%FT= FTSOL - FTIR
%semilogy(FT, PF)
%set(gca, 'YDir', 'reverse')
%xlabel('Net Flux (W/m^2)')
%ylabel('Pressure (bars)')




