
clear all ;clc
[ALT P T FH2O O3 THEAT TCOOL PSATCO2 FCO2]=textread('clima_last.tab','%f %f %f %f %f %f %f %f %f', 'headerlines', 1)

plot(THEAT, ALT,'g-x')
hold on
plot(TCOOL,ALT,'k-x')
xlabel('Heating rate(Kelvin/day)')
ylabel('Altitude (km)')
legend('Heating rate', 'Cooling rate')

