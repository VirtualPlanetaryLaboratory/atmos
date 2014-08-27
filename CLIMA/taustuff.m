clear all;clc

% for K1.eq.1 and standard atmosphere for 25, 51, 75, and 101 layers
format long
[TAU25 PATH25 INT25] =   textread('tau25', '%f %f %f', 55, 'headerlines', 1)
[TAU51 PATH51 INT51] =   textread('tau51', '%f %f %f', 55, 'headerlines', 1)
[TAU75 PATH75 INT75]=    textread('tau75', '%f %f %f', 55, 'headerlines', 1)
[TAU101 PATH101 INT101]=    textread('tau101', '%f %f %f', 55, 'headerlines', 1)

semilogy(INT25, TAU25,'k')
hold on

semilogy(INT51, TAU51)
hold on


semilogy(INT75, TAU75,'r')
hold on


semilogy(INT101, TAU101,'g')
legend('25 layers','51 layers', '75 layers', '101 layers')

%%
clear all;clc

% for K1.eq.1 and IUP=1 profile for 25, 51, 75, and 101 layers
format long
[TAU25 PATH25 INT25] =   textread('tau25', '%f %f %f', 55, 'headerlines', 1)
[TAU51 PATH51 INT51] =   textread('tau51', '%f %f %f', 55, 'headerlines', 1)
[TAU75 PATH75 INT75]=    textread('tau75', '%f %f %f', 55, 'headerlines', 1)
[TAU101 PATH101 INT101]=    textread('tau101', '%f %f %f', 55, 'headerlines', 1)

semilogy(INT25, TAU25,'k')
hold on

semilogy(INT51, TAU51)
hold on


semilogy(INT75, TAU75,'r')
hold on


semilogy(INT101, TAU101,'g')
legend('25 layers','51 layers', '75 layers', '101 layers')