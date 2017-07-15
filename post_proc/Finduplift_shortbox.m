% script to find uplift from T, pstar, mu, velocity 
% output files 

% Step 1: using output written in .vtu from Fenics code
% utilizing bash script by VTKtoMatlab.sh by S. Gold
% utilizing bash script by VTKtoMatlab_vec.sh by M. Roy
% based upon StreamlineTplot.m that reads in and analyzes T and pstar from
% the Fenics output.

% Mousumi Roy, March 4, 2014

clear all

colormap(gray)
% output interval
%
int = 5; %output interval
lw  = 1;  %default linewidth
%prange = [1.1, 1.5]*1.e9;
%gradprange = [-1000 1000];

%things to calculate flexural response
pm = 3300.0;
ps = 2500.0;
g  = 9.8;
E  = 50.0e9; % in N/m^2
v2 = 0.25 * 0.25;
%     Elastic plate thickness
Te = 35.0e3; % in m
%     Flexural rigidity of the lithospheric plate
Drig = (E *(Te^3.) ) / (12.0*(1-v2)); %in Nm
%     Alpha is the flexural parameter
alphaf = (4* Drig / (g*(pm-ps)))^(0.25) %in m

flag = 1; % flag to save/print figures

% first choose folder where data are:
Tb   = 1300;
muscale = 1.e+18;
ms   = ['1e+18'];%num2str(muscale);
mstr = ['mu\=' ms '/'];
mstrc= ['mu=' ms '/'];
Tbstr= ['Tb\=' num2str(Tb)];
Tbstrc= ['Tb=' num2str(Tb)];
locroot = ['tanhstep_smallbox_400km_no_adiabat/']

%locroot = ['cylinder_100/BLNK89/b12.7clog128/']
%locroot = ['test_steady_state/BLNK_LAB.5/longrun/']
%locroot = ['cylinder_50/BLNK/']/   
loc  = [locroot mstr Tbstr '/t6t']
loc2 = [locroot mstr Tbstr '/velocity']
loc3 = [locroot mstr Tbstr '/gradp']
loc4 = [locroot mstr Tbstr '/mu']
loc5 = [locroot mstr Tbstr '/pstar']
loc_cd = [locroot mstrc Tbstrc ]
%loc2 = ['HK03/' mstr Tbstr '/pstar']

%from MultipleRuns.py or codes like it in the folder(s) 
% above, we establish the Temp scale
%temp_values = [27.+273, Tb+273, 1300.+273, 1500.+273]
%dTemp = temp_values[3] - temp_values[0]

Tscale = 1500-27;
Tval   = 1300+273;
T_0    = 1300+273;
h      = 4e2; % box dimension in km
hscale = 4e5; % box scale in m
Hcrust = h - 30;  % crustal thickness in km; surface is at Yr = h km
tcont  = [300:100:1600];
pscale = 1192135725.0; % pressure scale in Pa from 
                        %MultipleRuns.py
rho_0 = 3300.;  % SI
rhof  = 3295;
alpha = 2.5e-5; % thermal expansion, SI
g     = 9.81;   % SI
mu_a  = muscale;
kappa_0 = 1.E-6;
%Ra    = rho_0*alpha*g*Tscale*(hscale^3)/(kappa_0*mu_a)
%vscale= rho_0*alpha*g*Tscale*(hscale^2)/mu_a;
Ra      = rho_0*alpha*g*Tscale*(hscale^3)/(kappa_0*muscale)
Rafac   = rho_0*alpha*g*Tscale*(hscale^3)/(kappa_0);
vscale  = rho_0*alpha*g*Tscale*(hscale^2)/muscale;
testpscale = muscale*vscale/hscale
pscale

% for streamline calculation, use the following from 
%paraview:  dimgradP = u*1192135725.0/1e6
% wvel = -(dimgradP-150*9.8*jHat)*1e-15/1e-2
% all SI units
drho = 150;
 
kovermu = 1e-15/1e-2;
startx = [1:25:990]'; %note here units must be in km as 
                       %displayed in box
starty = 500*ones(length(startx),1);
starty = 180*ones(length(startx),1);

figure(1); clf
figure(3); clf

for ii = 1:int:9
    clear dat val T X Y Z U V DPDX DPDY
    eval(['! sh VTKtoMatlab.sh ' loc '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = dat(:,3);
    
    clear dat
    eval(['! sh VTKtoMatlab_vec.sh ' loc2 '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    vx  = dat(:,3);
    vy  = dat(:,4);
    
    clear dat
    eval(['! sh VTKtoMatlab_vec.sh ' loc3 '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    dpdx  = dat(:,3);
    dpdy  = dat(:,4);

    eval(['! sh VTKtoMatlab.sh ' loc4 '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    mu = dat(:,3);
    
    eval(['! sh VTKtoMatlab.sh ' loc5 '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    pstar = dat(:,3);
    
    %F = TriScatteredInterp(x,y,val);
    %profile = pscale*F(posx,posy);
    
    %restructure the data into a matrix using the arrangement of points 
    % in the arrays x and y;
    % to do this, we first see that the arrays x and y are in order going
    % from left to right in the x-direction (constant y) starting at the
    % lower left corner of the mesh and the mesh is uniform with spacing 
    % dx and dy
    
    if ii==1 
        nn   = length(x);
        maxx = max(x);
        minx = min(x);
        delx   = diff(x);
        bigchange = minx - maxx;
        ind1 = 0;
        xind(1) = 1;
        yind(1) = 1;
    
%         ncol = find(delx == bigchange, 1);
%         nrow = sum((diff(y) ~= 0))+1;
%         dx   = h*(delx(1));
%         dy   = h*(y(ncol+1) - y(1));
%         xvec = [0:dx:max(x)];
%         yvec = [0:dy:max(y)];
        
        for i=2:nn-1
            xind(i) = xind(i-1) + 1;
            yind(i) = yind(i-1);
            if delx(i-1) == bigchange
                xind(i) = 1;
                yind(i) = yind(i-1)+1;
            end
        end
        xind(nn) = xind(nn-1) + 1;
        yind(nn) = yind(nn-1);
    end
    
    % now that the indexing is done, make the T array that contains the
    % temperature values for each uploaded timestep 
    
    for i=1:nn
        T(yind(i),xind(i)) = val(i); 
        X(yind(i),xind(i)) = x(i);
        Y(yind(i),xind(i)) = y(i); 
        U(yind(i),xind(i)) = vx(i); 
        V(yind(i),xind(i)) = vy(i); 
        DPDX(yind(i),xind(i)) = dpdx(i); 
        DPDY(yind(i),xind(i)) = dpdy(i);
        MU(yind(i),xind(i))   = mu(i);
        PSTAR(yind(i),xind(i))= pstar(i);
    end
%scale to real dimensional values
    Xr = 1000*X;  % in km
    Yr = 1000*Y;% in km
    Tr = Tscale*T;
    Ur = U*vscale;
    Vr = V*vscale;
    DPDY = pscale*DPDY/hscale;
    DPDX = pscale*DPDX/hscale;
    WY = -(DPDY + rhof*g)*kovermu;
    WX = -(DPDX)*kovermu;
    MU    = MU*muscale;
    PSTAR = pscale*PSTAR;
    
    if ii==1 
        Tr_init  = Tr;
        [ny, nx] = size(Tr_init);
        dy       = Yr(2,1)-Yr(1,1); %in km
        dx       = Xr(1,2)-Xr(1,1); %in km
        Ymantle  = (Yr < Hcrust) & (Yr >= 200 );% & (Tr <= Tval);
    end
% components of the stress tensor:
    [gradVrx gradVry]  = gradient(Vr,dy*1e3,dy*1e3);
    [gradUrx gradUry]  = gradient(Ur,dy*1e3,dy*1e3);
    
    Sxx = 2*MU.*gradUrx;
    Syy = 2*MU.*gradVry;
    Sxy = MU.*(gradVrx + gradUry);
    
    % test
%     divV   = gradUrx+gradVry; 
%     mesh(divV); %very small (order 4e-16 /s or smaller)
%     
    [dSxxdx dSxxdy] = gradient(Sxx, dy*1e3,dy*1e3);
    [dSyydx dSyydy] = gradient(Syy, dy*1e3,dy*1e3);
    [dSxydx dSxydy] = gradient(Sxy, dy*1e3,dy*1e3);
    
    %predicted pressure gradients from the calculated stress tensor
    gPx  = dSxxdx + dSxydy;
    gPy  = dSyydy + dSxydx;
    [psx psy] = gradient(PSTAR, dy*1e3,dy*1e3);
    

    test = Sxx + Syy; %should be small
    %[psx psy] = gradient(PSTAR, dy*1e3,dy*1e3);
    %figure(3); clf; mesh(test); pause
   
    WY = -(-psy - drho*g)*kovermu;
    WX = -(psx)*kovermu;

    % at each timestep, integrate T-Tinit to find the isostatic rock 
    % uplift (add up contributions below crust only, assuming crust density
    % unchanged and the alpha value is for mantle only):
    delT     = Tr - Tr_init;
    rhoarr   = rho_0*(1 - alpha*(Tr-T_0));
    integral = cumsum(delT.*Ymantle,1);
    ru_iso   = alpha*dy*1e3*integral(ny,:);
    ru_iso   = ru_iso - ru_iso(1);
    
    Plith    = g*cumsum(rhoarr,1,'reverse')*dy*1e3;
    
    Ptot = PSTAR + Plith;
    %pause
    % now find the dynamic rock uplift by considering tau_zz at the base of
    % the plate... see notes in magma_fluids_notes_3_13_14.pdf
    % first, find d u_y / dy
    
    %dVrdy    = -diff(Vr,1,1)*0.001/dy;
    %ru_array = (2*MU(1:ny-1,:).*dVrdy - PSTAR(1:ny-1,:))/fac;
    fac      = rho_0*g;
    %ru_test  = -filter(1,[1 -0.999], dSyydy*1e3*dy/fac);%
    ru_test  = cumsum(dSyydy*dy*1e3,1)/fac;
    %ru_test2 = -filter(1,[1 -0.999], dSxydx*1e3*dy/fac);%
    ru_test2 = cumsum(dSxydx*dy*1e3,1)/fac; 
    % need -ve sign as integrating down along -y-hat
    % note that the result does not depend on whether we integrate up or
    % down! the pressure field is conservative... same answer
    ru_test3  = cumsum(dSxxdx*dx*1e3,2)/fac;
    ru_test4  = cumsum(dSxydy*dx*1e3,2)/fac; 
    
    %pause;
    
    CC   = contour(Xr,Yr,Tr, [Tval, Tval]);
    cval = CC(1,1);  % contour value
    cpts = CC(2,1);  % number of points in contour -- this picks out the number 
    % for the first contour in the set; if number of points is less
    % than the size of CC, then there are multiple contours being
    % plotted... so we need to worry about the LAB contour... may have to
    % choose this by inspection (usually the one we want has the yvalues greater 
    % than 700 km?)  --- IN THIS CASE for shortbox ABOVE 350 KM?? 
    % first, find how all the contour lines that are being picked out:
    
     contind = find(CC(1,:)==cval); %indices of the beginnings of the contours
     ncont   = numel(contind); 
     maxy    = max(CC(2, contind+1));
     allinds = find(CC(2,:)==maxy);
     cind = allinds(1,1)-1;
     cval = CC(1,cind);
     cpts = CC(2,cind);
    
    cloc = [CC(1, cind+1:cind+cpts)' CC(2, cind+1:cind+cpts)']; % x,y locations of points
    Yrn  = Yr(1:ny-1,:);
    Xrn  = Xr(1:ny-1,:);
    
%     %figure(2); clf
%     cloc = [CC(1, 2:cpts)' CC(2, 2:cpts)']; % x,y locations of points
%     Yrn  = Yr(1:ny-1,:);
%     Xrn  = Xr(1:ny-1,:);
%     
    Y_isotherm = interp1(cloc(:,1),cloc(:,2),Xr(1,:));
    
    %FF     = scatteredInterpolant(Yr(:),Xr(:),ru_array(:));
    %FF     = scatteredInterpolant(Yrn(:),Xrn(:),ru_array(:));
    %FF      = scatteredInterpolant(Xr(:),Yr(:),ru_test(:));
    %ru_dyn1 = FF(cloc(:,1),cloc(:,2));
    %FF2     = scatteredInterpolant(Xr(:),Yr(:),ru_test2(:));
    %ru_dyn2 = FF2(cloc(:,1),cloc(:,2));
    % first interpolate to find rock uplift at contour location of LAB
    ru_dyn_a_cont  = interp2(Xr, Yr, ru_test,cloc(:,1),cloc(:,2));
    ru_dyn_b_cont  = interp2(Xr, Yr, ru_test2,cloc(:,1),cloc(:,2));
    % now interpolate the 1D function above onto the required Xr locations
    ru_dyn_a  = interp1(cloc(:,1), ru_dyn_a_cont, Xr(1,:));
    ru_dyn_b  = interp1(cloc(:,1), ru_dyn_b_cont, Xr(1,:));
    
    ru_dyn  = ru_dyn_a + ru_dyn_b;
    ru_dyn  = ru_dyn - ru_dyn(1,1);
    
    ru_dyn_a = ru_dyn_a - ru_dyn_a(1,1);
    ru_dyn_b = ru_dyn_b - ru_dyn_b(1,1);
    
    figure(2); clf; 
        plot(Xr(1,:),ru_dyn_b,'--','linewidth',[1.2],'color',[0 0 0.8]); hold on; pause(0.5);  
        plot(Xr(1,:),ru_dyn_a,'-','linewidth',[1.2],'color',[1 0 0]); hold on; pause(0.5);
        plot(Xr(1,:),ru_dyn,'-','linewidth',[1.2],'color',[0 0 0]); hold on; pause(0.5);
        legend('integrated shear stress gradient', 'vertical normal stress', 'total')
        set(gca,'fontsize',[14],'fontname','Helvetica')
    
    filename = ['DynTopo_' num2str(ii) ];
        WD1 = cd;
        cd(loc_cd)
        %eval(['print -dpng ' filename])
        savefig(filename);
        cd(WD1)
    %convolution to get flexural response for each part of the rock uplift
    % use Drig, the flexural rigidity and alphaf, the flexural parameter.
    % for 2D, the Green's function is a sum of sin and cosine mutliplied by exp:
    w_iso = zeros(1,length(Xr(1,:)));
    w_dyn = w_iso;
    norm  = w_iso;
    for jj = 1:length(Xr(1,:))
        xova  = 1e3*abs(Xr(1,:) - Xr(1,jj))/alphaf; %Xr in km
        Aiso  = ru_iso(1,jj);
        w_iso = w_iso + Aiso.*exp(-xova).*(sin(xova)+cos(xova));
        Adyn  = ru_dyn(1,jj);
        w_dyn = w_dyn + Adyn.*exp(-xova).*(sin(xova)+cos(xova));
        norm  = norm + exp(-xova).*(sin(xova)+cos(xova));
    end
    w_iso = w_iso./norm;
    w_dyn = w_dyn./norm;
    
    %output
    figure(1);clf
    %[cs h] = 
    contourf(Xr,Yr,Tr,tcont);hold on
    contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
    hanLAB=contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
%     surf(Xr,Yr,Tr);
%     view(2);  shading interp
    %clabel(cs, h, tcont,'fontsize',[12]);
    
    hold on
    plot(startx,starty,'wo'); 
    han = streamline(Xr,Yr,WX,WY,startx,starty); 
    set(han,'color','r','linewidth',[1.25]);
    quiver(Xr(1:2:end,1:2:end),Yr(1:2:end,1:2:end),U(1:2:end,1:2:end),V(1:2:end,1:2:end),3,'k'); 
        
    title([loc(1:4) ', t = ' num2str(ii) ' m.y.']);
    set(gca,'fontname','Helvetica','fontsize',[14])
    xlabel('km')
    box on
    %colorbar
    hold off
%     if flag == 1
%         filename = ['Streamplot_t_' num2str(ii) ];
%         WD1 = cd;
%         cd(loc_cd)
%         eval(['print -dpng ' filename])
%         cd(WD1)
%     end
    pause(0.5);  %pause
    
    figure(3); clf
    %subplot(411)
    %plot(Xr(1,:),ru_iso,'o'); 
    hold on
    plot(Xr(1,:),ru_iso,'--','linewidth',[1.2],'color',[0.8 0 0]); pause(0.5);  
    plot(Xr(1,:),w_iso,'-','linewidth',[2],'color',[1 0 0]); pause(0.5); %pause
    %subplot(412)
    %plot(Xr(1,:),ru_dyn,'ro'); hold on
    plot(Xr(1,:),ru_dyn,'--','linewidth',[1.2],'color',[0 0 0.8]); pause(0.5);  
    plot(Xr(1,:),w_dyn,'-','linewidth',[2],'color',[0 0 1]); pause(0.5); %pause
    %subplot(413)
    %plot(Xr(1,:),ru_dyn+ru_iso,'mo'); hold on
    plot(Xr(1,:),ru_dyn+ru_iso,'--','linewidth',[3],'color',[0.6 0.6 0.6]); pause(0.5);  
    plot(Xr(1,:),w_dyn+w_iso,'-','linewidth',[3],'color',[0 0 0]); pause(0.5);  %pause
    xlabel('km');ylabel('Rock Uplift, m');
    legend('Airy Isostasy','Flexure','Dynamic (no flexure)', 'Dynamic (with flexure)',...
        'Total (no flexure)', 'Total (with flexure)','Location','NorthEast');
    set(gca,'fontsize',[14],'fontname','Helvetica')
    
    filename = ['RockUplift_' num2str(ii) ];
        WD1 = cd;
        cd(loc_cd)
        %eval(['print -dpng ' filename])
        savefig(filename);cd(WD1)
        
    figure(4); clf
    hh = subplot(511);
    plot(Xr(1,:),w_dyn+w_iso,'-','linewidth',[3],'color',[0 0 0]); pause(0.5);  hold on
    plot([0 1000],[0 0],'--','linewidth',[1],'color',[0 0 0]);
    box off;
    
    pos = get(hh,'Position');
    new_h = axes('Position',pos);
    linkaxes([hh new_h],'y');
    pos(3) = eps;
    set(new_h,'Position',pos,'XTick',[],'XTickLabel',[]);
    set(gca,'YaxisLocation', 'right');%pause
    set(hh,'Visible','off');

    hold on
    %pause
    filename = ['TotalRockUplift_' num2str(ii) ];
        WD1 = cd;
        cd(loc_cd)
        savefig(filename);%eval(['print -dpng ' filename])
        cd(WD1)
        
    figure(1)
end
%pause
for ii = 10:int:99
    clear dat val profile
    eval(['! sh VTKtoMatlab.sh ' loc '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = dat(:,3);
    
    clear dat
    eval(['! sh VTKtoMatlab_vec.sh ' loc2 '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    vx  = dat(:,3);
    vy  = dat(:,4);
    
    clear dat
    eval(['! sh VTKtoMatlab_vec.sh ' loc3 '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    dpdx  = dat(:,3);
    dpdy  = dat(:,4);
    
    eval(['! sh VTKtoMatlab.sh ' loc4 '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    mu = dat(:,3);
    
    eval(['! sh VTKtoMatlab.sh ' loc5 '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    pstar = dat(:,3);
    
     for i=1:nn
        T(yind(i),xind(i)) = val(i); 
        X(yind(i),xind(i)) = x(i);
        Y(yind(i),xind(i)) = y(i); 
        U(yind(i),xind(i)) = vx(i); 
        V(yind(i),xind(i)) = vy(i); 
        DPDX(yind(i),xind(i)) = dpdx(i); 
        DPDY(yind(i),xind(i)) = dpdy(i);
        MU(yind(i),xind(i))   = mu(i);
        PSTAR(yind(i),xind(i))= pstar(i);
    end
   
 
%scale to real dimensional values
    Xr = 1000*X;  % in km
    Yr = 1000*Y;  % in km
    Tr = Tscale*T;
    Ur = U*vscale;
    Vr = V*vscale;
    DPDY = pscale*DPDY/hscale;
    DPDX = pscale*DPDX/hscale;
    WY = -(DPDY - rhof*g)*kovermu;
    WX = -(DPDX)*kovermu;
    MU    = MU*muscale;
    PSTAR = pscale*PSTAR;
    
    % components of the stress tensor:
    [gradVrx gradVry]  = gradient(Vr,dy*1e3,dy*1e3);
    [gradUrx gradUry]  = gradient(Ur,dy*1e3,dy*1e3);
    
    Sxx = 2*MU.*gradUrx;
    Syy = 2*MU.*gradVry;
    Sxy = MU.*(gradVrx + gradUry);
    
    % this is old and wrong
    %P   = 0.5*(Sxx+Syy); [px py] = gradient(P, dy*1e3,dy*1e3);
    
    [dSxxdx dSxxdy] = gradient(Sxx, dy*1e3,dy*1e3);
    [dSyydx dSyydy] = gradient(Syy, dy*1e3,dy*1e3);
    [dSxydx dSxydy] = gradient(Sxy, dy*1e3,dy*1e3);
    
    %predicted pressure gradients from the calculated stress tensor
    gPx  = dSxxdx + dSxydy;
    gPy  = dSyydy + dSxydx;
    [psx psy] = gradient(PSTAR, dy*1e3,dy*1e3);
    
    %test = Sxx + Syy;
    %[psx psy] = gradient(PSTAR, dy*1e3,dy*1e3);
    %figure(3); clf; mesh(test); pause(1)
   
    WY = -(-psy - drho*g)*kovermu;
    WX = -(psx)*kovermu;
    
    % at each timestep, integrate T-Tinit to find the isostatic rock 
    % uplift (add up contributions below crust only, assuming crust density
    % unchanged and the alpha value is for mantle only):
    delT   = Tr - Tr_init;
    rhoarr   = rho_0*(1 - alpha*(Tr-T_0));
    integral = cumsum(delT.*Ymantle,1);
    ru_iso = alpha*dy*1e3*integral(ny,:);
    ru_iso = ru_iso - ru_iso(1);
    
    % now find the dynamic rock uplift by considering tau_zz at the base of
    % the plate... see notes in magma_fluids_notes_3_13_14.pdf
    % first, find d u_y / dy
    dVrdy    = -diff(Vr,1,1)*0.001/dy;%dy is in km, Vr in m/s (order 1e-10)
    fac      = rho_0*g;
    ru_array = (2*MU(1:ny-1,:).*dVrdy - PSTAR(1:ny-1,:))/fac;
    %ru_test  = -(Syy)/fac;
    %ru_test2 = -cumsum(gPy*dy*1e3,1)/fac; % need -ve sign as integrating down along -y-hat
    ru_test  = -filter(1,[1 -0.999], dSyydy*1e3*dy/fac);%cumsum(dSyydy*dy*1e3,1)/fac;
    ru_test2 = -filter(1,[1 -0.999], dSxydx*1e3*dy/fac);%cumsum(dSxydx*dy*1e3,1)/fac; 
    % need -ve sign as integrating down along -y-hat
    % note that the result does not depend on whether we integrate up or
    % down! the pressure field is conservative... same answer
    
    CC   =contour(Xr,Yr,Tr, [Tval, Tval]);
    cval = CC(1,1);  % contour value
    cpts = CC(2,1);  % number of points in contour -- this picks out the number 
    % for the first contour in the set; if number of points is less
    % than the size of CC, then there are multiple contours being
    % plotted... so we need to worry about the LAB contour... may have to
    % choose this by inspection (usually the one we want has the yvalues greater 
    % than 700 km?)  --- IN THIS CASE for shortbox ABOVE 350 KM?? 
    % first, find how all the contour lines that are being picked out:
    
     contind = find(CC(1,:)==cval); %indices of the beginnings of the contours
     ncont   = numel(contind); 
     maxy    = max(CC(2, contind+1));
     allinds = find(CC(2,:)==maxy);
     cind = allinds(1,1)-1;
     cval = CC(1,cind);
     cpts = CC(2,cind);
    
    cloc = [CC(1, cind+1:cind+cpts)' CC(2, cind+1:cind+cpts)']; % x,y locations of points
    Yrn  = Yr(1:ny-1,:);
    Xrn  = Xr(1:ny-1,:);
    
    %FF     = scatteredInterpolant(Yr(:),Xr(:),ru_array(:));
    %FF     = scatteredInterpolant(Yrn(:),Xrn(:),ru_array(:));
%     FF      = scatteredInterpolant(Xr(:),Yr(:),ru_test(:));
%     ru_dyn1 = FF(cloc(:,1),cloc(:,2));
%     FF2     = scatteredInterpolant(Xr(:),Yr(:),ru_test2(:));
%     ru_dyn2 = FF2(cloc(:,1),cloc(:,2));
%     ru_dyn_a  = interp1(cloc(:,1),ru_dyn1,Xr(1,:));
%     ru_dyn_b  = interp1(cloc(:,1),ru_dyn2,Xr(1,:));
%     ru_dyn  = ru_dyn_b;% + ru_dyn_b;
%     ru_dyn  = ru_dyn - ru_dyn(1,1);
    % first interpolate to find rock uplift at contour location of LAB
    ru_dyn_a_cont  = interp2(Xr, Yr, ru_test,cloc(:,1),cloc(:,2));
    ru_dyn_b_cont  = interp2(Xr, Yr, ru_test2,cloc(:,1),cloc(:,2));
    % now interpolate the 1D function above onto the required Xr locations
    ru_dyn_a  = interp1(cloc(:,1), ru_dyn_a_cont, Xr(1,:));
    ru_dyn_b  = interp1(cloc(:,1), ru_dyn_b_cont, Xr(1,:));
   
    ru_dyn  = ru_dyn_a + ru_dyn_b;
    ru_dyn  = ru_dyn - ru_dyn(1,1);
    
    
    ru_dyn_a = ru_dyn_a - ru_dyn_a(1,1);
    ru_dyn_b = ru_dyn_b - ru_dyn_b(1,1);
    figure(2); clf; 
        plot(Xr(1,:),ru_dyn_b,'--','linewidth',[1.2],'color',[0 0 0.8]); hold on; pause(0.5);  
        plot(Xr(1,:),ru_dyn_a,'-','linewidth',[1.2],'color',[1 0 0]); hold on; pause(0.5);
        plot(Xr(1,:),ru_dyn,'-','linewidth',[1.2],'color',[0 0 0]); hold on; pause(0.5);
        legend('integrated shear stress gradient', 'vertical normal stress', 'total')
        set(gca,'fontsize',[14],'fontname','Helvetica')
    
    filename = ['DynTopo_' num2str(ii) ];
        WD1 = cd;
        cd(loc_cd)
        savefig(filename);%eval(['print -dpng ' filename])
        cd(WD1)
        
    %convolution to get flexural response for each part of the rock uplift
    % use Drig, the flexural rigidity and alphaf, the flexural parameter.
    % for 2D, the Green's function is a sum of sin and cosine mutliplied by exp:
    w_iso = zeros(1,length(Xr(1,:)));
    w_dyn = w_iso;
    norm  = w_iso;
   
    for jj = 1:length(Xr(1,:))
        xova  = 1e3*abs(Xr(1,:) - Xr(1,jj))/alphaf; %Xr in km, alphaf in m
        Aiso  = ru_iso(1,jj);
        w_iso = w_iso + Aiso.*exp(-xova).*(sin(xova)+cos(xova));
        Adyn  = ru_dyn(1,jj);
        w_dyn = w_dyn + Adyn.*exp(-xova).*(sin(xova)+cos(xova));
        norm  = norm + exp(-xova).*(sin(xova)+cos(xova));
    end
    w_iso = w_iso./norm;
    w_dyn = w_dyn./norm;
    
 %output
    figure(1);clf
    %[cs h] = 
    contourf(Xr,Yr,Tr,tcont);hold on
    contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
    %view(2);  shading interp
    %clabel(cs, h, tcont,'fontsize',[12]);
    
    hold on
    plot(startx,starty,'wo'); 
    han = streamline(Xr,Yr,WX,WY,startx,starty); 
    set(han,'color','r','linewidth',[1.25]);
    quiver(Xr(1:2:end,1:2:end),Yr(1:2:end,1:2:end),U(1:2:end,1:2:end),V(1:2:end,1:2:end),3,'k'); 
        
    title([loc(1:4) ', t = ' num2str(ii) ' m.y.']);
    set(gca,'fontname','Helvetica','fontsize',[14])
    xlabel('km')
    box on
    ax = gca;
    set(gca, 'ylim', [0 400]);
    ax.Clipping = 'off';
    %colorbar
    hold off
    
%      if flag == 1
%         filename = ['Streamplot_t_' num2str(ii) ];
%         WD1 = cd;
%         cd(loc_cd)
%         eval(['print -dpng ' filename])
%         cd(WD1)
%     end
%     pause(0.5);      
    
    figure(3); clf
    %subplot(411)
    %plot(Xr(1,:),ru_iso,'o'); 
    hold on
    plot(Xr(1,:),ru_iso,'--','linewidth',[1.2],'color',[0.8 0 0]); pause(0.5);  
    plot(Xr(1,:),w_iso,'-','linewidth',[2],'color',[1 0 0]); pause(0.5); %pause
    %subplot(412)
    %plot(Xr(1,:),ru_dyn,'ro'); hold on
    plot(Xr(1,:),ru_dyn,'--','linewidth',[1.2],'color',[0 0 0.8]); pause(0.5);  
    plot(Xr(1,:),w_dyn,'-','linewidth',[2],'color',[0 0 1]); pause(0.5); %pause
    %subplot(413)
    %plot(Xr(1,:),ru_dyn+ru_iso,'mo'); hold on
    plot(Xr(1,:),ru_dyn+ru_iso,'--','linewidth',[3],'color',[0.6 0.6 0.6]); pause(0.5);  
    plot(Xr(1,:),w_dyn+w_iso,'-','linewidth',[3],'color',[0 0 0]); pause(0.5);  %pause
    xlabel('km');ylabel('Rock Uplift, m');
    legend('Airy Isostasy','Flexure','Dynamic (no flexure)', 'Dynamic (with flexure)',...
        'Total (no flexure)', 'Total (with flexure)','Location','NorthEast');
    set(gca,'fontsize',[14],'fontname','Helvetica')
    
    filename = ['RockUplift_' num2str(ii) ];
        WD1 = cd;
        cd(loc_cd)
        savefig(filename);%eval(['print -dpng ' filename])
        cd(WD1)
        
    figure(4); clf
    hh = subplot(511);
    plot(Xr(1,:),w_dyn+w_iso,'-','linewidth',[3],'color',[0 0 0]); pause(0.5);  hold on
    plot([0 1000],[0 0],'--','linewidth',[1],'color',[0 0 0]);
    box off;
    
    pos = get(hh,'Position');
    new_h = axes('Position',pos);
    linkaxes([hh new_h],'y');
    pos(3) = eps;
    set(new_h,'Position',pos,'XTick',[],'XTickLabel',[]);
    set(gca,'YaxisLocation', 'right');%pause
    set(hh,'Visible','off');

    hold on
    %pause
    filename = ['TotalRockUplift_' num2str(ii) ];
        WD1 = cd;
        cd(loc_cd)
        savefig(filename);%eval(['print -dpng ' filename])
        cd(WD1)
    
    figure(1)
end

