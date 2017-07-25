%code to test effect of Delta_T on viscosty and density in BLNK
% MR -- April 3, 2015

    g      = 9.8 %N/kg
    rho0   = 3300; %kg/m^3
    L      = 100; % in km
    alpha  = 2.5e-5; % thermal expansion, SI
    TT    = 1573;
    kappa = 1e-6; %m^2/s

    h = 1000; % in km
    a = 100; % in km

    %Ep    = 0.014057;  %normalized activation energy
    adiabat = 0.25; % K / km
    %adiabat = 0;
    dTemp = 1300+adiabat*(h-L)
    mu_a  = 1e23;
    Tb    = 800 + 273;
    b     = 12.7;
    Ep    = b/dTemp
%b     = Ep*dTemp;
%Ep    = b/dTemp
%cc    = 0;
    cc    = log(128)
%Ep     = 0.014057;
%b     = 6.9;
    z=[1:h]'; % in km
    zl = find(z.*(z<=L)); % in km
    za = find(z.*(z>L)); % in km
    ind_l = find(z==L); % index for the LAB

    Tl = 273 + (Tb-273)*zl/L;
    Ta = TT + adiabat*(za - a);
    T = [Tl; Ta];
    %plot(T-273,z,'linewidth',[1.1]); hold on
    unitarr  = ones(length(T),1);
    dT       = (T-TT*unitarr);
    rho      = rho0*(unitarr - alpha*dT);
    
    
    mu   = mu_a*exp(-Ep*(T - TT) + cc*(z/h));
    
    figure(1);clf
    
    semilogx(mu,z)
    set(gca,'ydir','reverse','xlim', [1.e17,1.e30],'ylim',[0 200],'fontsize',[16]); grid on; 
    xlabel('\eta, Pa s')

    hold on
    semilogx(mu(ind_l-1:ind_l+1), z(ind_l-1:ind_l+1),'ro')
    
    delmu = mu(ind_l)/mu(ind_l+1)
    