[J P ALT T CONVEC DT TOLD FH20 FSAVE TCOOL THEAT]=textread('clima_allout1barnewcp.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 488);
[Jold Pold ALTold Told CONVECold DTold TOLDold FH20old FSAVEold TCOOLold THEATold]=textread('clima_allout1baroldcp.tab','%f %f %f %f %f %f %f %f %f %f %f', 101, 'headerlines', 638);




plot(T(6:101),ALT(6:101), Told(6:101), ALT(6:101))
title('Kasting et al. (1984) at 1 bar')
xlabel('T(K)')
ylabel('ALTITUDE(km)')
legend('New Cps', 'Old Cps')
pause


plot(T(80:101),ALT(80:101), Told(80:101), ALT(80:101))
title('Kasting et al. (1984) at 1 bar')
xlabel('T(K)')
ylabel('ALTITUDE(km)')
legend('New Cps', 'Old Cps')