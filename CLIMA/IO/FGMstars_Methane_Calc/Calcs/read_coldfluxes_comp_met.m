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

[JF PFF ALTF FTOTALF FTIRF FDNIRF FUPIRF FTSOLF FDNSOLF FUPSOLF DIVFF] = textread('F7200/COLD_EARLYMARS/clima_allout1e-2met.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 856)
[JM PFM ALTM FTOTALM FTIRM FDNIRM FUPIRM FTSOLM FDNSOLM FUPSOLM DIVFM] = textread('M3600/M3600COLD/clima_allout1e-2met.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 2356)
[JS PFS ALTS FTOTALS FTIRS FDNIRS FUPIRS FTSOLS FDNSOLS FUPSOLS DIVFS] = textread('COLDSUN_EarlyMars/clima_allout1barwithmethane.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 1156)

FUPSOLF = FUPSOLF/1000 % convert to W/m^2
FUPIRF  = FUPIRF/1000 % convert to W/m^2

FUPSOLM = FUPSOLM/1000 % convert to W/m^2
FUPIRM  = FUPIRM/1000 % convert to W/m^2

FUPSOLS = FUPSOLS/1000 % convert to W/m^2
FUPIRS  = FUPIRS/1000 % convert to W/m^2

%Ours
%subplot 121
semilogy(FUPSOLF, PFF, 'k')
set(gca, 'YDir', 'reverse')
%axis([15 30 .1E-5 1])
xlabel('Upward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FUPSOLM, PFM, 'r')
set(gca, 'YDir', 'reverse')
%axis([15 30 .1E-5 1])
xlabel('Upward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FUPSOLS, PFS)
set(gca, 'YDir', 'reverse')
%axis([15 30 .1E-5 1])
xlabel('Upward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

pause
hold off

%subplot 122
semilogy(FUPIRF, PFF, 'k')
set(gca, 'YDir', 'reverse')
%axis([110 230 .1E-5 1])
xlabel('Upward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FUPIRM, PFM, 'r')
set(gca, 'YDir', 'reverse')
%axis([110 230 .1E-5 1])
xlabel('Upward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FUPIRS, PFS)
set(gca, 'YDir', 'reverse')
%axis([110 230 .1E-5 1])
xlabel('Upward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on


%% 
% Downward Fluxes
clear all;clc

[JF PFF ALTF FTOTALF FTIRF FDNIRF FUPIRF FTSOLF FDNSOLF FUPSOLF DIVFF] = textread('F7200/COLD_EARLYMARS/clima_allout1e-2met.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 856)
[JM PFM ALTM FTOTALM FTIRM FDNIRM FUPIRM FTSOLM FDNSOLM FUPSOLM DIVFM] = textread('M3600/M3600COLD/clima_allout1e-2met.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 2356)
[JS PFS ALTS FTOTALS FTIRS FDNIRS FUPIRS FTSOLS FDNSOLS FUPSOLS DIVFS] = textread('COLDSUN_EarlyMars/clima_allout1barwithmethane.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 1156)


FDNSOLF = FDNSOLF/1000 % convert to W/m^2
FDNIRF  = FDNIRF/1000 % convert to W/m^2

FDNSOLM = FDNSOLM/1000 % convert to W/m^2
FDNIRM  = FDNIRM/1000 % convert to W/m^2

FDNSOLS = FDNSOLS/1000 % convert to W/m^2
FDNIRS  = FDNIRS/1000 % convert to W/m^2


%Ours
%subplot 121
semilogy(FDNSOLF, PFF, 'k')
set(gca, 'YDir', 'reverse')
%axis([80 115 .1E-5 1])
xlabel('Downward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FDNSOLM, PFM, 'r')
set(gca, 'YDir', 'reverse')
%axis([80 115 .1E-5 1])
xlabel('Downward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FDNSOLS, PFS)
set(gca, 'YDir', 'reverse')
%axis([80 115 .1E-5 1])
xlabel('Downward Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

pause
hold off

%subplot 122
semilogy(FDNIRF, PFF, 'k')
set(gca, 'YDir', 'reverse')
axis([-5 200 .1E-5 2])
xlabel('Downward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FDNIRM, PFM, 'r')
set(gca, 'YDir', 'reverse')
axis([-5 200 .1E-5 2])
xlabel('Downward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FDNIRS, PFS)
set(gca, 'YDir', 'reverse')
axis([-5 200 .1E-5 2])
xlabel('Downward IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on




%% Net Fluxes
clear all;clc

[JF PFF ALTF FTOTALF FTIRF FDNIRF FUPIRF FTSOLF FDNSOLF FUPSOLF DIVFF] = textread('F7200/COLD_EARLYMARS/clima_allout1e-2met.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 856)
[JM PFM ALTM FTOTALM FTIRM FDNIRM FUPIRM FTSOLM FDNSOLM FUPSOLM DIVFM] = textread('M3600/M3600COLD/clima_allout1e-2met.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 2356)
[JS PFS ALTS FTOTALS FTIRS FDNIRS FUPIRS FTSOLS FDNSOLS FUPSOLS DIVFS] = textread('COLDSUN_EarlyMars/clima_allout1barwithmethane.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 1156)

FTSOLF = FTSOLF/1000 % convert to W/m^2
FTIRF  = abs(FTIRF/1000) % convert to W/m^2

FTSOLM = FTSOLM/1000 % convert to W/m^2
FTIRM  = abs(FTIRM/1000) % convert to W/m^2

FTSOLS = FTSOLS/1000 % convert to W/m^2
FTIRS  = abs(FTIRS/1000) % convert to W/m^2



%Ours
%subplot 121
semilogy(FTSOLF, PFF, 'k')
set(gca, 'YDir', 'reverse')
%axis([65 85 .1E-5 1])
xlabel('Net Incoming Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FTSOLM, PFM, 'r')
set(gca, 'YDir', 'reverse')
%axis([65 85 .1E-5 1])
xlabel('Net Incoming Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FTSOLS, PFS)
set(gca, 'YDir', 'reverse')
%axis([65 85 .1E-5 1])
xlabel('Net Incoming Solar Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

pause
hold off



%subplot 122
semilogy(FTIRF, PFF, 'k')
set(gca, 'YDir', 'reverse')
%axis([20 130 .1E-5 1])
xlabel('Net Outgoing IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FTIRM, PFM, 'r')
set(gca, 'YDir', 'reverse')
%axis([20 130 .1E-5 1])
xlabel('Net Outgoing IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on

semilogy(FTIRS, PFS)
set(gca, 'YDir', 'reverse')
%axis([20 130 .1E-5 1])
xlabel('Net Outgoing IR Flux (W/m^2)')
ylabel('Pressure (bars)')
hold on



%pause
%FT= FTSOL - FTIR
%semilogy(FT, PF)
%set(gca, 'YDir', 'reverse')
%xlabel('Net Flux (W/m^2)')
%ylabel('Pressure (bars)')




