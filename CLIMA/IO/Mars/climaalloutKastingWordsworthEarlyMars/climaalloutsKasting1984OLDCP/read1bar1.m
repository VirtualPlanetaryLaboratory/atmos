[J P ALT T CONVEC DT TOLD FH20 FSAVE TCOOL THEAT]=textread('clima_allout1bar.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 488);
[J P ALT T CONVEC DT TOLD FH20 FSAVE TCOOL THEAT]=textread('clima_allout1bar.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 488);




plot(T(6:101),ALT(6:101))
title('Kasting et al. (1984) at 1 bar')
xlabel('T(K)')
ylabel('ALTITUDE(km)')
xlegend('Old Cps', 'New Cps')