%test solomatov scaling
%MR -- 10/28/15
TTop = 273;
TBottom = 1300 + 273;

g       = 9.8; %N/kg
rho0    = 3300; %kg/m^3
alpha   = 2.5e-5; % thermal expansion, SI
kappa_0 = 1.E-6; %m^2/s
mu      = [1.e18,1.e20,1.e21];
DT      = (TBottom-TTop); %200; %C

h=400e3;
b=12.7;

mu_int=8e18;

Rafac   = rho0*alpha*g*DT*(h^3)/(kappa_0);
Ra = Rafac./mu

Ra_int = Rafac/mu_int

zlinf = h*(b^(4/3))*(Ra_int^(-1/3))*1e-3

Ra_c=20.9*b^4
mu_c = Rafac/Ra_c

mu=1e19
vscale  = rho0*alpha*g*DT*(h^2)/mu

Pe = 1.66e-9*h/kappa_0

tc = (0.01*h*h/kappa_0)/3e7