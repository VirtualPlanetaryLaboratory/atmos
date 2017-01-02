% This program reads in Temperature(K), Pressure(bars), and Altitude(bars)
% and computes the temperature profiles as a function of pressure and
% altitude.
% This version overlays the CO2 saturation vapor pressure curve
% This file reads in the clima_last.tab file
clear all;clc

[ALT PTOT T FH2O O3 THEAT TCOOL PSATCO2 FCO2 ] = textread('clima_last3e-2.tab','%f %f %f %f %f %f %f %f %f',101,'headerlines',1)


PCO2 = FCO2.*PTOT 

for j=1:length(T)
    if T(j)<216.56
PSL(j) = 6.760956 - 1284.07/(T(j) - 4.718) + 1.256E-4*(T(j) - 143.15);
    else
PSL(j) = 3.128082 - 867.2124/T(j) + 1.865612E-2*T(j) - 7.248820E-5*T(j)^2 + 9.3E-8*T(j)^3;
    end
end
PATM = 10.^PSL;
PSCO2 = 1.013*PATM;

n=101  %Number of layers down starting from top of atmosphere
% Plots T vs Altitude
plot(T(1:n), ALT(1:n))
xlabel('Temperature(K)')
ylabel('Altitude(km)')

hold on
%pause

% Plots T vs. PCO2 in logspace  
%plot(T(1:n),log10(PCO2(1:n)),T(1:n),log10(PSCO2(1:n)),'r-.') 
%set(gca, 'YDir', 'reverse')
%xlabel('Temperature(K)')
%ylabel('CO2 partial pressure(bars)')
%legend('Partial pressure curve', 'CO2 saturation vapor pressure curve')



%%

% This program reads in Temperature(K), Pressure(bars), and Altitude(bars)
% and computes the temperature profiles as a function of pressure and
% altitude.
% This version overlays the water vapor saturation pressure curve

clear all;clc

[ALT PTOT T FH2O O3 THEAT TCOOL PSATCO2 FCO2] = textread('clima_last.tab','%f %f %f %f %f %f %f %f %f',101,'headerlines',1)

% Constants(in cgs units)
R= 1.9872
RV=R/18
SUBL= 677
TOP = 273.15
POP= 6.103E-3

% TTAB vector
TTAB(1) = 1e-2
TTAB(2:75) = [5:5:370]
TTAB(76)= 373.97

%PVAP vector
PVAP = [6.116560e-003 , 8.725050e-003 , 1.228000e-002 , 1.705380e-002 , 2.338590e-002 ,...
    3.168740e-002 , 4.245150e-002 , 5.626370e-002 , 7.381260e-002 , 9.589780e-002 ,...
    1.234460e-001 , 1.575200e-001 , 1.993280e-001 , 2.502340e-001 , 3.117710e-001 ,...
    3.856460e-001 , 4.737520e-001 , 5.781750e-001 , 7.012010e-001 , 8.453250e-001 ,...
    1.013250e+000 , 1.207910e+000 , 1.432440e+000 , 1.690230e+000 , 1.984860e+000 ,...
    2.320170e+000 , 2.700230e+000 , 3.129320e+000 , 3.611970e+000 , 4.152930e+000 ,...
    4.757190e+000 , 5.429950e+000 , 6.176650e+000 , 7.002950e+000 , 7.914720e+000 ,...
    8.918060e+000 , 1.001930e+001 , 1.122490e+001 , 1.254160e+001 , 1.397640e+001 ,...
    1.553640e+001 , 1.722890e+001 , 1.906150e+001 , 2.104200e+001 , 2.317810e+001 ,...
    2.547820e+001 , 2.795040e+001 , 3.060350e+001 , 3.344590e+001 , 3.648690e+001 ,...
    3.973540e+001 , 4.320080e+001 , 4.689290e+001 , 5.082140e+001 , 5.499650e+001 ,...
    5.942850e+001 , 6.412830e+001 , 6.910680e+001 , 7.437560e+001 , 7.994640e+001 ,...
    8.583150e+001 , 9.204390e+001 , 9.859700e+001 , 1.055050e+002 , 1.127830e+002 ,...
    1.204470e+002 , 1.285140e+002 , 1.370030e+002 , 1.459330e+002 , 1.553280e+002 ,...
    1.652120e+002 , 1.756140e+002 , 1.865680e+002 , 1.981180e+002 , 2.103270e+002 , 2.205820e+002];

PH2O = FH2O.*PTOT

for j = 1:length(T)
    
    if (T(j)<=273.16)
        HL=SUBL
        PSH2O(j) = POP*exp(-HL/RV*(1/T(j)-1/TOP))
    elseif (T(j)>273.16)&(T(j)<=646.96)
        TC(j) = T(j)-273.15;
        N(j) = round(TC(j)/5)+1;
        FR(j) = (TC(j) - TTAB(N(j)))/(TTAB(N(j)+1)-TTAB(N(j)));
        PSH2O(j) = FR(j)*PVAP(N(j)+1) + (1 -FR(j))*PVAP(N(j));
    elseif (T(j)>646.96)
        PSH2O= 1E99
    end
end

% Plots T vs. Altitude
n=101   %Number of layers down starting from top of atmosphere 
plot(T(1:n), ALT(1:n))
xlabel('Temperature(K)')
ylabel('Altitude(km)')
pause

% Plots T vs. PH2O in logspace

plot(T(1:n),log10(PH2O(1:n)),T(1:n),log10(PSH2O(1:n)),'r-.') 
set(gca, 'YDir', 'reverse')
xlabel('Temperature(K)')
ylabel('H2O partial pressure(bars)')
legend('Partial pressure curve', 'H2O saturation vapor pressure curve')

%% Plots Heating rates vs. altitude
clear all;clc

[ALT PTOT T FH2O O3 THEAT TCOOL PSATCO2 FCO2 ] = textread('clima_last1e-2met.tab','%f %f %f %f %f %f %f %f %f',101,'headerlines',1)

% Plots heating rate vs. Altitude
n=101   %Number of layers down starting from top of atmosphere 
plot(THEAT(1:n), ALT(1:n),'--')
hold on
plot([0,70],[ALT(69),ALT(69)],'--')  % x vector, then y vector
xlabel('Heating rate (K/day)')
ylabel('Altitude(km)')

hold on
%pause