% script to plot color plots of T, with velocity arrows, and streamlines in Matlab

% Step 1: using output written in .vtu from Fenics code
% utilizing bash script by VTKtoMatlab.sh by S. Gold
% utilizing bash script by VTKtoMatlab_vec.sh by M. Roy

% Mousumi Roy, Feb 7, 2014; revised Oct 10, 2015

clear all

colormap(gray)

%time at which we want the output plotted
tout = 9
timeint = 1; %if writing output every 2 my instead of 1 my, need this to be 2 (for longruns)

%index at which we want the output plotted
outind = tout/timeint;

% output interval
%
int = 5; %output interval
lw  = 1;  %default linewidth
%prange = [1.1, 1.5]*1.e9;
%gradprange = [-1000 1000];

flag = 0; % flag to save/print figures
long = 0; % if longer run, then do last for-loop

%Tbs = [800 1000 1300];
Tbs = [1300];
% first choose folder where data are:%
%locroot = ['test_steady_state/BLNK_LAB.9/longrun/']
%locroot = ['test_steady_state/BLNK_LAB.9/shortbox/']
%locroot = ['test_steady_state/BLNK_LAB.7/shortbox/']
%locroot = ['tanhstep/BLNK/w0.2h0.05sc0.05/']
%locroot = ['cylinder_100/BLNK89/b12.7clog128/']
%locroot = ['cylinder_100/HK03/']
locroot = ['tanhstep_smallbox_400km_no_adiabat/']
%locroot = ['longruns_no_adiabat/']

mstr = ['mu\=5e+18/'];
mstrc= ['mu=5e+18/'];

%Define constants
rho_0 = 3300.;  % SI
alpha = 2.5e-5; % thermal expansion, SI
g     = 9.81;   % SI
kappa_0 = 1.E-6;

%from MultipleRuns.py or codes like it in the folder(s) above, we 
%establish the Temp scale
%temp_values = [27.+273, Tb+273, 1300.+273, 1500.+273]
%dTemp = temp_values[3] - temp_values[0]
Tscale = 1500-27;
Tval   = 1573;
h      = 1e3; % box dimension in km

hscale = 1e6; % box scale in m
tcont  = [300:100:1600]; %in K
pscale = 1192135725.0; % pressure scale from MultipleRuns.py in Pa
muscale = str2num(mstrc(4:8)); % in Pa s
Ra      = rho_0*alpha*g*Tscale*(hscale^3)/(kappa_0*muscale)
Rafac   = rho_0*alpha*g*Tscale*(hscale^3)/(kappa_0);
vscale  = rho_0*alpha*g*Tscale*(hscale^2)/muscale;

% for streamline calculation, use the following from paraview:
% dimgradP = u*1192135725.0/1e6
% wvel = -(dimgradP-150*9.8*jHat)*1e-15/1e-2
% all SI units
drhomax = 500;
drhomin = 100;
rhomelt = 2800; %kg/m^3
drho    = 400;%[drhomin:drhomax];
kovermu = 1e-15/1e-2;
startx  = [1:20:990]'; %note here units must be in km as displayed in box
starty  = 180*ones(length(startx),1);

figure(2);clf;
figure(1);clf;colormap(gray)
numt = 1;
    
ii=1
combTracers = [];

%     function trackStream(newStreams, n, lower, upper)
%         numStreams = length(newStreams);
%         xEndPoints = zeros(1,numStreams);
%         for index = 1:numStreams
%             xData = get(newStreams(index), 'XData');
%             xEndPoints(index) = xData(end);
%         end
%         figure(10); clf; hold on;
%         subplot(2,1,1); hold on; title('Bin Dist at Current step');
%         hist(xEndPoints, n); xlim([lower,upper]);
%         
%         combTracers = [combTracers; histc(xEndPoints, linspace(lower,upper,n))];
%         subplot(2,1,2); hold on; title('Cumulative distribution');
%         area(combTracers'); xlim([1,n]);
%     end
%for 
    Tbcount = 1:length(Tbs)
    Tb   = Tbs(Tbcount);
    Tbstr= ['Tb\=' num2str(Tb)];
    Tbstrc= ['Tb=' num2str(Tb)];
    %locroot = ['tanhstep/BLNK/w0.2h0.05sc0.05/']
    %locroot = ['test_steady_state/BLNK_LAB.6/']
    loc  = [locroot mstr Tbstr '/t6t']
    loc2 = [locroot mstr Tbstr '/velocity']
    loc3 = [locroot mstr Tbstr '/gradp']
    loc4 = [locroot mstr Tbstr '/mu']
    loc_cd = [locroot mstrc Tbstrc ]
    %loc2 = ['HK03/' mstr Tbstr '/pstar']
%%
    for ii = 1%:int:9
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

        clear dat
        eval(['! sh VTKtoMatlab.sh ' loc4 '00000' num2str(ii) '.vtu'])
        dat = load('PythonSoln');
        mu = dat(:,3);
        
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
            T(yind(i),xind(i))  = val(i);
            MU(yind(i),xind(i)) = mu(i);
            X(yind(i),xind(i)) = x(i);
            Y(yind(i),xind(i)) = y(i); 
            U(yind(i),xind(i)) = vx(i); 
            V(yind(i),xind(i)) = vy(i); 
            DPDX(yind(i),xind(i)) = dpdx(i); 
            DPDY(yind(i),xind(i)) = dpdy(i); 
        end
    %scale to real dimensional values
        Xr = h*X; 
        Yr = (h*Y);
        Tr = Tscale*T;
        rho  = rho_0*(1 - alpha*(Tr - Tval));
        drho = rho - rhomelt;
        DPDY = pscale*DPDY/hscale;
        DPDX = pscale*DPDX/hscale;
        WY = -(DPDY - drho*g)*kovermu;
        WX = -(DPDX)*kovermu;
        MU = muscale*MU;
        U  = U*vscale;
        V  = V*vscale;
        Vmeltx = WX + U;
        Vmelty = WY + V;
        %output
        figure(1);
        %clf
        %[cs h] = 
         contourf(Xr,Yr,Tr,tcont);hold on
         contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
        figure(2);
         [clith] = contour(Xr,Yr,Tr, [Tval, Tval]);hold on
         zll    = mean(clith(2,2:clith(2,1)));
         plot([min(min(Xr)), max(max(Xr))],[zll,zll],'k--');
        
   %     %make an array containting the mean depth of the Tval contour as a fn
   %     %of time
          % find mean viscosity in the convecting interior  - used to find Ra_i
          mu_int         = mean(MU(Yr<zll));
          Ra_int         = Rafac/mu_int
          meanz(numt,:)  = [ii*timeint, zll, Ra_int];
          numt          = numt + 1; 
           
    %     surf(Xr,Yr,Tr);
    %     view(2);  shading interp
    %     clabel(cs, h, tcont,'fontsize',[12]);
    
        figure(1)% now overlay streamlines and velocity vectors
        hold on
        plot(startx,starty,'wo'); 
        han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
        set(han,'color','r','linewidth',[1.25]);
        quiver(Xr,Yr,U,V,3,'k'); 
        %pause
        
        title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14])
        xlabel('km');ylabel('km'); 
        box on
        %colorbar
        hold off
        
%         if flag == 1
%             filename = ['Streamplot_t_' num2str(ii) ];
%             WD1 = cd;
%             cd(loc_cd)
%             eval(['print -dpdf ' filename])
%             cd(WD1)
%         end
        
        %find pressure-gradients along average isotherm depth, zll
        close = min(min(abs(Yr-zll)));
        xpos  = Xr(abs(Yr-zll)==close);
        p_x   = DPDX(abs(Yr-zll)==close); %profile of dpdx along zll
        %output
    
        %output to file for each timestep
%             filename = ['dpdx_zll_t_' num2str(ii)];
%             dat = [xpos p_x];
%             WD1 = cd;
%             cd(loc_cd)
%             eval(['save ' filename ' dat -ascii'])
%             cd(WD1)
            
        pause(0.5);  %pause
        
       
    end
%%
%     if outind <10 
        ii=outind;
        clear dat val profile
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
        
        clear dat
        eval(['! sh VTKtoMatlab.sh ' loc4 '00000' num2str(ii) '.vtu'])
        dat = load('PythonSoln');
        mu = dat(:,3);
        
         for i=1:nn
            T(yind(i),xind(i)) = val(i); 
            MU(yind(i),xind(i)) = mu(i);
            X(yind(i),xind(i)) = x(i);
            Y(yind(i),xind(i)) = y(i); 
            U(yind(i),xind(i)) = vx(i); 
            V(yind(i),xind(i)) = vy(i); 
            DPDX(yind(i),xind(i)) = dpdx(i); 
            DPDY(yind(i),xind(i)) = dpdy(i); 
        end


    %scale to real dimensional values
        Xr = h*X; 
        Yr = h*Y;
        Tr = Tscale*T;
        DPDY = pscale*DPDY/hscale;
        DPDX = pscale*DPDX/hscale;
        rho  = rho_0*(1 - alpha*(Tr - Tval));
        drho = rho - rhomelt;
        WY = -(DPDY - drho*g)*kovermu;
        WX = -(DPDX)*kovermu;
        MU = muscale*MU;
        U  = U*vscale;
        V  = V*vscale;
        Vmeltx = WX + U;
        Vmetly = WY + V;
     %output
        figure(1);clf
        %[cs h] = 
        subplot('Position',[0.1 0.25 0.8 0.7])
        contourf(Xr,Yr,Tr,tcont);hold on
        contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
%          figure(2);
%          [clith] = contour(Xr,Yr,Tr, [Tval, Tval]);hold on
%           zll    = mean(clith(2,2:clith(2,1)));
%           plot([min(min(Xr)), max(max(Xr))],[zll,zll],'k--');
%           
    %     %make an array containting the mean depth of the Tval contour as a fn
    %     %of time
          % find mean viscosity in the convecting interior  - used to find Ra_i
          mu_int         = mean(MU(Yr<zll));
          Ra_int         = Rafac/mu_int
          meanz(numt,:)  = [ii*timeint, zll, Ra_int];
          numt          = numt + 1; 
        %view(2);  shading interp
        %clabel(cs, h, tcont,'fontsize',[12]);

        plot(startx,starty,'wo'); hold on
        han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
        set(han,'color','r','linewidth',[1.25]);
        quiver(Xr,Yr,U,V,3,'k'); 

        %title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14],'xlim',[0 1000],'ylim',[200 1000])
        ylabel('km'); 
        box on
        %colorbar
        hold off
        
        
        
%         if flag == 1
%             filename = ['Streamplot_t_' num2str(ii) ];
%             WD1 = cd;
%             cd(loc_cd)
%             eval(['print -dpdf ' filename])
%             cd(WD1)
%         end
        pause(0.5);      
        %find pressure-gradients along average isotherm depth, zll
        close = min(min(abs(Yr-zll)));
        xpos  = Xr(abs(Yr-zll)==close);
        p_x   = DPDX(abs(Yr-zll)==close);
        figure(1);
        subplot('Position',[0.1 0.1 0.8 0.1])
        plot(xpos,p_x/(rho_0*g),'k-');
        set(gca,'fontname','Helvetica','fontsize',[14],'ylim',[-0.005 0.005])
        xlabel('km'); %ylabel('Pa/m')
%         %output to file
%             filename = ['dpdx_zll_t_' num2str(ii)];
%             dat = [xpos p_x];
%             WD1 = cd;
%             cd(loc_cd)
%             eval(['save ' filename ' dat -ascii'])
%             cd(WD1)
%         figure(1)

        
        %make a figure for paper at a chosen time...
            
%             disp('write output?')
%             pause
%             WD1 = cd;
%             cd(loc_cd)
% %                figure(1); title('')
% %                print -dpdf no_convection_20my_melt.pdf
% %                savefig('no_convection_20my_melt.fig')
% %                figure(4)
% %                print -dpdf no_convection_20my_dpdx.pdf
% %                savefig('no_convection_20my_dpdx.fig')
%                figure(1); title('')
%                print -dpdf with_convection_8my_melt.pdf
%                savefig('with_convection_8my_melt.fig')
% %                figure(4)
% %                print -dpdf with_convection_98my_dpdx.pdf
% %                savefig('with_convection_98my_dpdx.fig')
%             cd(WD1)
%             
      %end

    %pause
%%
    %for ii = 10:int:99
    if outind >=9 & outind <= 99
        ii=outind;
        clear dat val profile
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
        
        clear dat
        eval(['! sh VTKtoMatlab.sh ' loc4 '00000' num2str(ii) '.vtu'])
        dat = load('PythonSoln');
        mu = dat(:,3);
        
       
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
         
        
                
         for i=1:nn
            T(yind(i),xind(i)) = val(i); 
            MU(yind(i),xind(i)) = mu(i);
            X(yind(i),xind(i)) = x(i);
            Y(yind(i),xind(i)) = y(i); 
            U(yind(i),xind(i)) = vx(i); 
            V(yind(i),xind(i)) = vy(i); 
            DPDX(yind(i),xind(i)) = dpdx(i); 
            DPDY(yind(i),xind(i)) = dpdy(i); 
        end


    %scale to real dimensional values
        Xr = h*X; 
        Yr = h*Y;
        Tr = Tscale*T;
        DPDY = pscale*DPDY/hscale;
        DPDX = pscale*DPDX/hscale;
        rho  = rho_0*(1 - alpha*(Tr - Tval));
        drho = rho - rhomelt;
        WY = -(DPDY - drho*g)*kovermu;
        WX = -(DPDX)*kovermu;
        MU = muscale*MU;
        U  = U*vscale;
        V  = V*vscale;
        Vmeltx = WX + U;
        Vmetly = WY + V;
     %output
         
        
        figure(1);clf
        %[cs h] = 
        %subplot('Position',[0.1 0.25 0.8 0.7])
        %colormap('jet')
        contourf(Xr,Yr,Tr,tcont);hold on
        hanLAB=contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
%          figure(2);
%          [clith] = contour(Xr,Yr,Tr, [Tval, Tval]);hold on
%           zll    = mean(clith(2,2:clith(2,1)));
%           plot([min(min(Xr)), max(max(Xr))],[zll,zll],'k--');
%           
    %     %make an array containting the mean depth of the Tval contour as a fn
    %     %of time
          % find mean viscosity in the convecting interior  - used to find Ra_i
          mu_int         = mean(MU(Yr<zll));
          Ra_int         = Rafac/mu_int
          meanz(numt,:)  = [ii*timeint, zll, Ra_int];
          numt          = numt + 1; 
        %view(2);  shading interp
        %clabel(cs, h, tcont,'fontsize',[12]);

        %plot(startx,starty,'wo'); 
        hold on
        han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
        set(han,'color','r','linewidth',[1.25]);
        quiver(Xr(1:2:end,1:2:end),Yr(1:2:end,1:2:end),U(1:2:end,1:2:end),V(1:2:end,1:2:end),3,'k'); 

        %title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14],'xlim',[0 1000],'ylim',[0 1000])
        ylabel('km'); 
        box on
        %colorbar
        hold off

        figure(6);clf;colormap(gray)%plot gradpx
            subplot(311)
            contourf(Xr,Yr,DPDX/(rho_0*g),10);hold on
            plot(clith(1,2:end),clith(2,2:end),'y-','linewidth',[2])
            strh = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
            set(strh,'color','r','linewidth',[1.25]);
            set(gca,'fontname','Helvetica','fontsize',[14],'xlim',[0 1000],'ylim',[150 400])
            ch=colorbar
            caxis([-0.1 0.1])
        %pause
 
%         figure(4);clf;
% %         subplot(211);
%           colormap(gray)
% %         contourf(Xr,Yr,Tr,tcont);hold on
% %         contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
% %         plot(startx,starty,'wo'); hold on
% %         han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
% %         set(han,'color','r','linewidth',[1.25]);
% %         quiver(Xr,Yr,U,V,3,'k'); 
% %         set(gca,'fontname','Helvetica','fontsize',[14],'xlim',[0 1000],'ylim',[600 1000])
% %         ylabel('km'); colorbar
% %         box on
%         %colorbar
%         
%         subplot(211);%colormap(copper)
%         contourf(Xr,Yr,DPDX./(rho_0*g));hold on;colorbar
%         plot(startx,starty,'wo'); hold on
%         %contour(Xr,Yr,Tr, [Tval, Tval],'y','linewidth',[2]);
%         plot(hanLAB(1,2:end),hanLAB(2,2:end),'y','linewidth',[2]);
%         han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
%         set(han,'color','r','linewidth',[1.25]);
%         title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
%         set(gca,'fontname','Helvetica','fontsize',[14],'xlim',[0 1000],'ylim',[600 1000])
%         xlabel('km');ylabel('km'); 
%         box on
%     
%         if flag == 1
%             filename = ['Streamplot_t_rev' num2str(ii) ];
%             WD1 = cd;
%             cd(loc_cd)
%             eval(['print -dpng ' filename])
%             cd(WD1)
%         end
%         pause(0.5);      
%         %find pressure-gradients along average isotherm depth, zll
%         close = min(min(abs(Yr-zll)));
%         xpos  = Xr(abs(Yr-zll)==close);
%         p_x   = DPDX(abs(Yr-zll)==close);
%         figure(1);
% %         subplot('Position',[0.1 0.1 0.8 0.1])
% %         plot(xpos,p_x/(rho_0*g),'k-');
% %         set(gca,'fontname','Helvetica','fontsize',[14],'ylim',[-0.005 0.005])
% %         xlabel('km'); %ylabel('Pa/m')
% % %         %output to file
% %             filename = ['dpdx_zll_t_' num2str(ii)];
% %             dat = [xpos p_x];
% %             WD1 = cd;
% %             cd(loc_cd)
% %             eval(['save ' filename ' dat -ascii'])
% %             cd(WD1)
% %         figure(1)
% 
%         figure(5);subplot(211)
%         
%         %make a figure for paper at a chosen time...
%             
%             disp('write output?')
%             pause
             WD1 = cd;
             cd(loc_cd)
%                 figure(1); title('')
%                 filename = ['Streamplot_t_rev' num2str(ii) ];
%                 eval(['print -dpng ' filename])
% %                print -dpdf no_convection_20my_melt.pdf
% %                savefig('no_convection_20my_melt.fig')
% %                figure(4)
% %                print -dpdf no_convection_20my_dpdx.pdf
% %                savefig('no_convection_20my_dpdx.fig')
% %                figure(1); title('')
%                 print -dpdf with_convection_10my_760_Tb1300_melt_rev.pdf
%                 savefig('with_convection_10my_760_Tb1300_melt_rev.fig')
%                figure(4)
%                 %print -dpdf dpdx_stream_10my_rev.pdf
%                 %savefig('dpdx_stream_10my_rev.fig')
%                figure(6)
%                filename = ['dpdx_contours_t_' num2str(ii) ];
%                  filename = ['Streamplot_t_' num2str(ii) ];
% %                figure(1); title('')
%                 eval(['print -dpng ' filename])
% %                %print -dpdf with_convection_398my_melt.pdf
               %savefig('with_convection_398my_melt.fig')
%                figure(4)
%                print -dpdf with_convection_1164my_dpdx.pdf
%                savefig('with_convection_1164my_dpdx.fig')
            
            cd(WD1)
            
      end

%% comment here 
    
    if long==1
    %if we have long runs, uncomment the following for loop
      if outind <= 800 & outind >= 100
        ii = outind;
        clear dat val profile
        eval(['! sh VTKtoMatlab.sh ' loc '000' num2str(ii) '.vtu'])
        dat = load('PythonSoln');
        x = dat(:,1);
        y = dat(:,2);
        val = dat(:,3);

        clear dat
        eval(['! sh VTKtoMatlab_vec.sh ' loc2 '000' num2str(ii) '.vtu'])
        dat = load('PythonSoln');
        vx  = dat(:,3);
        vy  = dat(:,4);

        clear dat
        eval(['! sh VTKtoMatlab_vec.sh ' loc3 '000' num2str(ii) '.vtu'])
        dat = load('PythonSoln');
        dpdx  = dat(:,3);
        dpdy  = dat(:,4);

        clear dat
        eval(['! sh VTKtoMatlab.sh ' loc4 '000' num2str(ii) '.vtu'])
        dat = load('PythonSoln');
        mu = dat(:,3);
        
         for i=1:nn
            T(yind(i),xind(i)) = val(i); 
            MU(yind(i),xind(i)) = mu(i);
            X(yind(i),xind(i)) = x(i);
            Y(yind(i),xind(i)) = y(i); 
            U(yind(i),xind(i)) = vx(i); 
            V(yind(i),xind(i)) = vy(i); 
            DPDX(yind(i),xind(i)) = dpdx(i); 
            DPDY(yind(i),xind(i)) = dpdy(i); 
        end


    %scale to real dimensional values
        Xr = h*X; 
        Yr = h*Y;
        Tr = Tscale*T;
        DPDY = pscale*DPDY/hscale;
        DPDX = pscale*DPDX/hscale;
        rho  = rho_0*(1 - alpha*(Tr - Tval));
        drho = rho - rhomelt;
        WY = -(DPDY - drho*g)*kovermu;
        WX = -(DPDX)*kovermu;
        MU = muscale*MU;
        U  = U*vscale;
        V  = V*vscale;
        Vmeltx = WX + U;
        Vmetly = WY + V;
     %output
        figure(1);clf
        %[cs h] = 
        %subplot('Position',[0.1 0.25 0.8 0.7])
        contourf(Xr,Yr,Tr,tcont);hold on
        contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
%          figure(2);
%          [clith] = contour(Xr,Yr,Tr, [Tval, Tval]);hold on
%           zll    = mean(clith(2,2:clith(2,1)));
%           plot([min(min(Xr)), max(max(Xr))],[zll,zll],'k--');
%           
    %     %make an array containting the mean depth of the Tval contour as a fn
    %     %of time
          % find mean viscosity in the convecting interior  - used to find Ra_i
          mu_int         = mean(MU(Yr<zll));
          Ra_int         = Rafac/mu_int
          meanz(numt,:)  = [ii*timeint, zll, Ra_int];
          numt          = numt + 1; 
        %view(2);  shading interp
        %clabel(cs, h, tcont,'fontsize',[12]);

        plot(startx,starty,'wo'); hold on
        han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
        set(han,'color','r','linewidth',[1.25]);
        %quiver(Xr,Yr,U,V,3,'k'); 
        quiver(Xr(1:2:end,1:2:end),Yr(1:2:end,1:2:end),U(1:2:end,1:2:end),V(1:2:end,1:2:end),3,'k'); 

        %title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14],'xlim',[0 1000],'ylim',[0 1000])
        ylabel('km'); 
        box on
        %colorbar
        hold off

         if flag == 1
            figure(1)
            filename = ['Streamplot_t_' num2str(ii) ];
            WD1 = cd;
            cd(loc_cd)
            %eval(['print -dpdf ' filename])
            cd(WD1)
        end
        pause(0.5);      
        %find pressure-gradients along average isotherm depth, zll
        close = min(min(abs(Yr-zll)));
        xpos  = Xr(abs(Yr-zll)==close);
        p_x   = DPDX(abs(Yr-zll)==close);
        figure(1);
% %         subplot('Position',[0.1 0.1 0.8 0.1])
% %         plot(xpos,p_x/(rho_0*g),'k-');
% %         set(gca,'fontname','Helvetica','fontsize',[14],'ylim',[-0.005 0.005])
% %         xlabel('km'); %ylabel('Pa/m')
%         %output to file
%             filename = ['dpdx_zll_t_' num2str(ii)];
%             dat = [xpos p_x];
%             WD1 = cd;
%             cd(loc_cd)
%             eval(['save ' filename ' dat -ascii'])
%             cd(WD1)
%         figure(1)
        
        %make a figure for paper at a chosen time...
            
            %disp('write output?')
            %pause
            WD1 = cd;
            cd(loc_cd);
                filename = ['dpdx_contours_t_' num2str(ii) ];
%                filename = ['Streamplot_t_' num2str(ii) ];
%                figure(1); title('')
                set(ch,'ylim',[-0.1 0.1])
                %eval(['print -dpng ' filename])
%                %print -dpdf with_convection_398my_melt.pdf
               %savefig('with_convection_398my_melt.fig')
%                figure(4)
%                print -dpdf with_convection_1164my_dpdx.pdf
%                savefig('with_convection_1164my_dpdx.fig')
            cd(WD1)
      end
    end
    %%
    
  