% script to find the analytic prediction for steady-state thermal boundary
% layer thickness from eqn 6.329 in turcotte and schubert

% Mousumi Roy, Sept 13, 2014

clear all

Ra_c = 657.5;
TTop = 273;
TBottom = 1300 + 273;

g       = 9.8; %N/kg
rho0    = 3300; %kg/m^3
alpha   = 2.5e-5; % thermal expansion, SI
kappa_0 = 1.E-6; %m^2/s
mu      = [1.e19,1.e20,1.e21,1.e22,1.e23,1.e24];
DT      = (TBottom-TTop); %200; %C

y_L = ((Ra_c*2*mu*kappa_0)/(rho0*alpha*g*DT)).^(1/3);
y_L = y_L*1e-3  %in km
fac = y_L(2)/y_L(1)
fac = y_L(3)/y_L(2)
mu=6.5e+21;
Z = 1e-3*((Ra_c*2*mu*kappa_0)/(rho0*alpha*g*DT)).^(1/3)
mu=9.3e+21;
Z = 1e-3*((Ra_c*2*mu*kappa_0)/(rho0*alpha*g*DT)).^(1/3)
mu=4.7e+21;
Z = 1e-3*((Ra_c*2*mu*kappa_0)/(rho0*alpha*g*DT)).^(1/3)

mu=1e+18;
Z = 1e-3*((Ra_c*2*mu*kappa_0)/(rho0*alpha*g*DT)).^(1/3)
mu=2.2e+18;
Z = 1e-3*((Ra_c*2*mu*kappa_0)/(rho0*alpha*g*DT)).^(1/3)
mu=3.1e+18;
Z = 1e-3*((Ra_c*2*mu*kappa_0)/(rho0*alpha*g*DT)).^(1/3)
