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

flag = 0; % flag to save/print figures

% first choose folder where data are:
Tb   = 1300;
muscale = 1.e+21;
ms   = ['1e+21'];%num2str(muscale);
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
starty = 200*ones(length(startx),1);

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
    
    test = Sxx + Syy; %should be small
    [psx psy] = gradient(PSTAR, dy*1e3,dy*1e3);
    %figure(3); clf; mesh(test); pause
   
    WY = -(gPy - drho*g)*kovermu;
    WX = -(gPx)*kovermu;

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
    
    
end

fgPy = flipud(gPy);
fgPx = flipud(gPx);
%pause
subx = gPx(3:end-3,3:end-3); %choose region away from edges
suby = gPy(3:end-3,3:end-3);
subXr = Xr(3:end-3,3:end-3);
subYr = Yr(3:end-3,3:end-3);

qx = cumtrapz(subx,2)*dx*1e3;
qy = cumtrapz(suby,1)*dy*1e3;

 subqx = qx;
 subqy = qy;
 y_add = qy(:,1)*ones(1,length(qx(1,:)));
 x_add = ones(length(qy(:,1)),1)*qx(1,:);
%pause
% path integral along top edge (x_add) and left edge (y_add) of model is added in
 subqx = qx + y_add - x_add; %- x_add;% - y_add;
 slope = ((subqx(:,end)-subqx(:,1))/(max(max(subXr))))*ones(1,length(qx(:,1)));
 trend = ones(1,length(qx(:,1)))*subqx(:,1) + slope.*subXr;
 subqx = subqx - trend;
 %pause;

colormap(jet)
figure(2);clf
mesh(subXr,subYr,subqx);
colormap(jet)
figure(3);clf
mesh(subXr, subYr, subqy);

%pause;
figure(1)
zout = [328];
zout = [328,96,24];
%zout=[24:8:380]
for i = 1:length(zout)
    clear irow icol ind
    [irow] = find(subYr(:,1)> zout(i)-1 & subYr(:,1) < zout(i)+1);
    inddd = irow(1);
%     datx = subXr(inddd,:);
%     daty = subqx(inddd,:);
%     [coeff] = polyfit(datx,daty,1);
%     trend = coeff(1)*datx+coeff(2);
%     detrendy = daty - trend;
    %plot(subXr(inddd-1,:),subqx(inddd-1,:),'.','markersize',[12]);hold on
    %plot(datx,trend,'m-');
    plotno = ['3' '1' num2str(i)];
    subplot(plotno)
    plot(subXr(inddd,:),subqx(inddd,:)-mean(subqx(inddd,:)),'.','markersize',[12]);hold on
    %plot(subXr(inddd,:),trend(inddd,:),'m.','markersize',[12]);hold on
    %plot(subXr(inddd+1,:),subqx(inddd+1,:),'.','markersize',[12]);hold on
    %plot(subXr(inddd-1,:),subqy(inddd-1,:),'--')
    plot(subXr(inddd,:),subqy(inddd,:)-mean(subqy(inddd,:)),'--')
    %plot(subXr(inddd+1,:),subqy(inddd+1,:),'--')
    tt = ['Depth z=',num2str(400-zout(i)), ' km'];
    title(tt)
    pause(0.5)
    %clf
end
%savefig('ProofPressure1e18.fig')
pause(2); clf
clear qx qy subx suby

%start with PSTAR, see if you can recover path independence 

subx = psx(3:end-3,3:end-3); %choose region away from edges
suby = psy(3:end-3,3:end-3);
subXr = Xr(3:end-3,3:end-3);
subYr = Yr(3:end-3,3:end-3);

qx = cumtrapz(subx,2)*dx*1e3;
qy = cumtrapz(suby,1)*dy*1e3;

 subqx = qx;
 subqy = qy;
 y_add = qy(:,1)*ones(1,length(qx(1,:)));
 x_add = ones(length(qy(:,1)),1)*qx(1,:);
%pause
% path integral along top edge (x_add) and left edge (y_add) of model is added in
 subqx = qx + y_add - x_add; %- x_add;% - y_add;
 slope = ((subqx(:,end)-subqx(:,1))/(max(max(subXr))))*ones(1,length(qx(:,1)));
 trend = ones(1,length(qx(:,1)))*subqx(:,1) + slope.*subXr;
 subqx = subqx - trend;
 %pause;

colormap(jet)
figure(2);clf
mesh(subXr,subYr,subqx);
colormap(jet)
figure(3);clf
mesh(subXr, subYr, subqy);

%pause;
figure(1)
zout = [328];
zout = [328,96,24];
%zout=[24:8:380]
for i = 1:length(zout)
    clear irow icol ind
    [irow] = find(subYr(:,1)> zout(i)-1 & subYr(:,1) < zout(i)+1);
    inddd = irow(1);
%     datx = subXr(inddd,:);
%     daty = subqx(inddd,:);
%     [coeff] = polyfit(datx,daty,1);
%     trend = coeff(1)*datx+coeff(2);
%     detrendy = daty - trend;
    %plot(subXr(inddd-1,:),subqx(inddd-1,:),'.','markersize',[12]);hold on
    %plot(datx,trend,'m-');
    plotno = ['3' '1' num2str(i)];
    subplot(plotno)
    plot(subXr(inddd,:),subqx(inddd,:)-mean(subqx(inddd,:)),'.','markersize',[12]);hold on
    %plot(subXr(inddd,:),trend(inddd,:),'m.','markersize',[12]);hold on
    %plot(subXr(inddd+1,:),subqx(inddd+1,:),'.','markersize',[12]);hold on
    %plot(subXr(inddd-1,:),subqy(inddd-1,:),'--')
    plot(subXr(inddd,:),subqy(inddd,:)-mean(subqy(inddd,:)),'--')
    %plot(subXr(inddd+1,:),subqy(inddd+1,:),'--')
    tt = ['Depth z=',num2str(400-zout(i)), ' km'];
    title(tt)
    pause(0.5)
    %clf
end
difx = gPx(3:end-3,3:end-3)-subx;
dify = gPy(3:end-3,3:end-3)-suby;
figure(3);clf;mesh(subXr, subYr, difx/(max(max(psx))))
figure(2);clf;mesh(subXr, subYr, dify/(max(max(psy))))
