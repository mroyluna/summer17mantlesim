% calculate Ra number
% MR -- Feb 2016 - for JGR revisions v2

clear all

ymax = 400;

g      = 9.8 %N/kg
rho0   = 3300 %kg/m^3
L      = 100; % in km
alpha  = 2.5e-5 % thermal expansion, SI
TT    = 1573;
kappa = 1e-6 %m^2/s

h = 1000; % in km
a = 100; % in km

%Ep    = 0.014057;  %normalized activation energy
%adiabat = 0.25; % K / km
adiabat = 0;
dTemp = 1300+adiabat*(h-L)
%dTemp = 1300+(1500-1300)/(h-L)
%put in parameters for wet diffusion creep from Hirth and Kohlstedt 2003
R = 8.3145; % universal gas constant J/K/mol
E = 375.e3; % activation energy in J/mol 
V = 1.e-6; % activation volume in m^3/mol
Ep = E/(R*dTemp)

C = rho0*g*V*h*1e3/(R*dTemp)
K = 2e-13;

mu_a = [1.e18,1.e19,1.e20,1.e21,5.e21];%,1.e23];

Ra_a = rho0*alpha*g*dTemp*h*h*h*1e9./(kappa*mu_a)