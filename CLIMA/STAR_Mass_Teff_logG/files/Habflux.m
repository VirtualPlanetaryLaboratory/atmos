clc
clear



seffsun=[1.76320E+00   1.11420E+00   1.00710E+00   3.43800E-01   2.98100E-01];
a=[2.62400E-04   1.66020E-04   1.49880E-04   5.89420E-05   5.11050E-05];
b=[2.34510E-08   1.79010E-08   1.33950E-08   1.65580E-09   1.43550E-09];
c=[-1.34390E-11  -9.06330E-12  -7.67660E-12  -3.00450E-12  -2.60500E-12];
d=[-2.76200E-15  -2.01540E-15  -1.57770E-15  -5.29830E-16  -4.59380E-16];

teff = [3498 3498 3600 5518 4720 5780 5214];
tstar = teff-5780;

seffobserve = seffsun(3) + a(3)*tstar + b(3)*tstar.^2 + c(3)*tstar.^3 + d(3)*tstar.^4;
plot(seffobserve,teff,'o')
hold on

teff = 2600:200:7200;
 
   tstar = teff-5780;
   seffihz = seffsun(3) + a(3)*tstar + b(3)*tstar.^2 + c(3)*tstar.^3 + d(3)*tstar.^4;
   seffohz = seffsun(4) + a(4)*tstar + b(4)*tstar.^2 + c(4)*tstar.^3 + d(4)*tstar.^4;
   teff = teff + 200;


plot(seffihz,teff,'r')
plot(seffohz,teff,'b')
hold off