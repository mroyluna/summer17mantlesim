function StreamlineTplot
% script to plot color plots of T, with velocity arrows, and streamlines in Matlab

% Step 1: using output written in .vtu from Fenics code
% utilizing bash script by VTKtoMatlab.sh by S. Gold
% utilizing bash script by VTKtoMatlab_vec.sh by M. Roy

% Mousumi Roy, Feb 7, 2014 - revised Oct 2015

clear all
set(0,'DefaultFigureWindowStyle','docked')
colormap(gray)

figure(6); clf

interval = 2; %output interval
lw  = 1;  %default linewidth
%prange = [1.1, 1.5]*1.e9;
%gradprange = [-1000 1000];

flag = 1; % flag to save/print figures
long = 1; % if longer run, then do last for-loop
timeint = 2; %if writing output every 2 my instead of 1 my, need this to be 2 if longrun
endind = 500;

%Tbs = [800 1000 1300];
Tbs = [1300];
% first choose folder where data are:%
%locroot = '/home/alex/Desktop/shortbox/'
%locroot = ['test_steady_state/BLNK_LAB.9/shortbox/']
%locroot = ['test_steady_state/BLNK_LAB.7/shortbox/']
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
Tscale = 1305-27;
Tval   = 1563;
h      = 1e3; % box dimension in km

hscale = 1e6; % box scale in m
tcont  = [100:100:1600]; %in K
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
drho    = mean([drhomin:drhomax]);
kovermu = 1e-13/1;
streamint = 20;
startx  = [1:streamint:990]'; %note here units must be in km as displayed in box
nstream = length(startx)*0.5;
starty  = 200*ones(length(startx),1);
prevDist = [];

figure(2);clf;
figure(1);clf;
colormap(gray)
numt = 1;

zeroString   = '000000';
combTracers  = [];
stackTracers = [];

    function trackStream(newStreams, n, lower, upper, ind)
        numStreams = length(newStreams);
        xEndPoints = zeros(1,numStreams);
        for index = 1:numStreams
            xData = get(newStreams(index), 'XData');
            xEndPoints(index) = xData(end);
        end
        set(0,'DefaultFigureWindowStyle','docked')
        figure(5); clf; hold on;
        subplot(2,1,1); hold on; title('Bin Dist at Current step');
        hist(xEndPoints, n); xlim([lower,upper]);
        
        combTracers = [combTracers; histc(xEndPoints, linspace(lower,upper,n))];
        subplot(2,1,2); hold on; title('Cumulative distribution');
        area(combTracers'); xlim([1,n]);
        
        if ind == 31
            figure(6); 
            colormap(hot);shading faceted
            
            subplot(3,2,1);  %title('Cumulative distribution');
            xrange = linspace(lower,upper,n);
            plotx  = xrange(3:end-3);
            indsout = [1, 3, 5, 7];%, 9, 11]; %if going to 50 my, keep 9, 11
            ploty  = combTracers(indsout,3:end-3);
            area(plotx,ploty'); xlim([lower,upper]);
            set(gca,'fontname','Helvetica','fontsize',[14]);
            set(gca,'xlim',[0 1000],'yTickLabel',' ','box','off')
            
            clear ploty;
            ploty  = combTracers(:,3:end-3);
            subplot(3,2,2);
            m1 = mean(ploty(1:2,:));
            m2 = mean(ploty(3:6,:));
            m3 = mean(ploty(7:end,:));
            %normalize to 1
            m1 = m1/max(m3);m2 = m2/max(m3); m3 = m3/max(m3);
            plot(plotx, movingmean(m1',15),'color',[0 0 0],'linewidth',[2]); hold on
            plot(plotx, movingmean(m2',15),'color',[1 0 0],'linewidth',[2]); 
%             plot(plotx, movingmean(m3',15),'color',[0 0.2 0.8],'linewidth',[2]); 
            xlim([lower,upper]);
            legend('0-10 my', '10-30 my','30-50 my', 'location','EastOutside')
            set(gca,'fontname','Helvetica','fontsize',[14]);
            set(gca,'xlim',[0 1000],'ylim',[0 2])
            
        end
    end

for Tb = Tbs
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
    for ind = 1:interval:endind
        %pause
        clear dat val T X Y Z U V DPDX DPDY
        numZeros = int64(length(zeroString)-length(num2str(ind)));
        fileSuffix = [zeroString(1:numZeros) num2str(ind)];
        
        eval(['! sh VTKtoMatlab.sh ' loc fileSuffix '.vtu'])
        dat = load('PythonSoln');
        x = dat(:,1);
        y = dat(:,2);
        val = dat(:,3);
        
        clear dat
        eval(['! sh VTKtoMatlab_vec.sh ' loc2 fileSuffix '.vtu'])
        dat = load('PythonSoln');
        vx  = dat(:,3);
        vy  = dat(:,4);
        
        clear dat
        eval(['! sh VTKtoMatlab_vec.sh ' loc3 fileSuffix '.vtu'])
        dat = load('PythonSoln');
        dpdx  = dat(:,3);
        dpdy  = dat(:,4);
        
        clear dat
        eval(['! sh VTKtoMatlab.sh ' loc4 fileSuffix '.vtu'])
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
        
        if ind==1
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
        divmelt = divergence(Xr, Yr, Vmeltx, Vmelty);
        divsolid = divergence(Xr, Yr, U, V);
        %         figure(5);clf
        %         subplot(121); mesh(Xr, Yr, log10(abs(divmelt)));
        %         subplot(122); mesh(Xr, Yr, log10(abs(divsolid)));
        %         figure(6); clf;
        %         quiver(Xr,Yr,Vmeltx,Vmelty,3,'k');set(gca,'ylim',[700 1000]);hold on
        %         %contour(Xr,Yr, sqrt(Vmeltx.*Vmeltx + Vmelty.*Vmelty),100);
        %         plot(startx,starty,'wo');
        %         [verts] = stream2(Xr,Yr,Vmeltx,Vmelty,startx,starty);
        %         sl = streamline(verts);
        %         set(sl,'Color','r');
        %         axis tight manual off;
        %         ax = gca;
        %         ax.Position = [0,0,1,1];
        %         %streamparticles(verts,5,'ParticleAlignment','on')
        %         streamparticles(verts,5,'ParticleAlignment','on')%,'animate',[1]);pause
        
        %         figure(7);clf;subplot(311)
        %         plot(Xr(end,:),Vmelty(end,:),'k-');hold on
        %         plot(Xr(end,:),Vmelty(1,:),'r-');hold on
        %         set(gca,'fontname','Helvetica','fontsize',[14])
        %         xlabel('km'); ylabel('m/s')
        %pause;
        
        %output
        figure(1); clf
        %[cs h] =
        dx = Xr(1,2)-Xr(1,1);
        dy = Yr(2,1)-Yr(1,1);
        Xn = [min(min(Xr)):dx*1.5:max(max(Xr))];
        Yn = [min(min(Yr)):dy*1.5:max(max(Yr))];
        [Xf, Yf] = meshgrid(Xn,Yn);
        %sTr = interp2(Xr, Yr, Tr, Xf, Yf,'cubic');
        
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
        meanz(numt,:)  = [ind*timeint, zll, Ra_int];
        numt          = numt + 1;
        
        %     surf(Xr,Yr,Tr);
        %     view(2);  shading interp
        %     clabel(cs, h, tcont,'fontsize',[12]);
        
        figure(1)% now overlay streamlines and velocity vectors
        hold on
        plot(startx,starty,'wo');
        han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty);
        
        set(han,'color','r','linewidth',[1.25]);
        quiver(Xr(1:2:end,1:2:end),Yr(1:2:end,1:2:end),U(1:2:end,1:2:end),V(1:2:end,1:2:end),3,'k');
        %pause
        
        title([loc(1:4) ', t = ' num2str(timeint*ind) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14])
        xlabel('km');ylabel('km');
        box on
        ax = gca;
        set(gca, 'ylim', [0 400]);
        ax.Clipping = 'off';
        
        trackStream(han, nstream, min(min(Xr)), max(max(Xr)), ind);
        %colorbar
        hold off
        
        
%         figure(10)
%         prevDist = trackStream(combDist, han, 20,min(min(Xr)),max(max(Xr)));
        
        %find pressure-gradients along average isotherm depth, zll
        close = min(min(abs(Yr-zll)));
        xpos  = Xr(abs(Yr-zll)==close);
        p_x   = DPDX(abs(Yr-zll)==close); %profile of dpdx along zll
        %output
        figure(4);clf;subplot(311)
        plot(xpos,p_x,'k-');
        set(gca,'fontname','Helvetica','fontsize',[14],'ylim',[-150 150])
        xlabel('km'); ylabel('Pa/m')
        %output to file for each timestep
        filename = ['dpdx_zll_t_' num2str(ind)];
        dat = [xpos p_x];
        WD1 = cd;
        cd(loc_cd)
        eval(['save ' filename ' dat -ascii'])
        cd(WD1)
        
        pause(0.5);  %pause
        
        if flag == 1
            WD1 = cd
            cd(loc_cd)
            figure(1); title('');
            chan=colorbar
            set(chan,'ydir','reverse')
            brighten(1,.05)
            caxis([27 1350])
            filename = ['Streamplot_t_' num2str(ind) ];
            eval(['print -dpdf ' filename])
            cd(WD1)
        end
    end
    
    %%
     if flag == 1
            WD1 = cd;
            cd(loc_cd)
            figure(1); title('');
            colorbar
            chan=colorbar
            set(chan,'ydir','reverse')
            caxis([27 1350])
            brighten(1,.05)
            filename = ['Streamplot_t_' num2str(ind) ];
            eval(['print -dpdf ' filename])
            figure(6); 
            subplot(321);
            get(gca,'ylim');hold on;
            ylim = [0 30];
            plot([250 250], ylim, 'k--','linewidth',[2.0]);
            plot([750 750], ylim, 'k--' ,'linewidth',[2.0]); 
            set(gca,'xlim',[0 1000],'ylim',ylim);
            filename = ['Meltout' num2str(ind) '_' mstrc(4:8)];
            eval(['print -dpdf ' filename])
            cd(WD1)
    end
end
end