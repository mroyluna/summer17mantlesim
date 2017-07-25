% script to plot color plots of T, with velocity arrows, and streamlines in Matlab

% Step 1: using output written in .vtu from Fenics code
% utilizing bash script by VTKtoMatlab.sh by S. Gold
% utilizing bash script by VTKtoMatlab_vec.sh by M. Roy

% Mousumi Roy, Feb 7, 2014; revised Oct 10, 2015

clear all

colormap(gray)

%time at which we want the output plotted
tout = 199
timeint = 2; %if writing output every 2 my instead of 1 my, need this to be 2 (for longruns)

%index at which we want the output plotted
outind = 380;%tout/timeint;

% output interval
%
int = 5; %output interval
lw  = 1;  %default linewidth
%prange = [1.1, 1.5]*1.e9;
%gradprange = [-1000 1000];

flag = 0; % flag to save/print figures
long = 1; % if longer run, then do last for-loop

%Tbs = [800 1000 1300];
Tbs = [1300];
% first choose folder where data are:%
%locroot = ['test_steady_state/BLNK_LAB.9/longrun/']
%locroot = ['test_steady_state/BLNK_LAB.9/shortbox/']
%locroot = ['test_steady_state/BLNK_LAB.7/shortbox/']
%locroot = ['tanhstep/BLNK/w0.2h0.05sc0.05/']
%locroot = ['cylinder_100/BLNK89/b12.7clog128/']
%locroot = ['cylinder_100/HK03/']
%locroot = ['tanhstep_smallbox_400km_no_adiabat/']
locroot = ['longruns_no_adiabat/test2/']

mstr = ['mu\=1e+18/'];
mstrc= ['mu=1e+18/'];
muscale = str2num(mstrc(4:8)); % in Pa s

%Define constants
rho_0 = 3300.;  % SI
rhof  = 3295;
alpha = 2.5e-5; % thermal expansion, SI
g     = 9.81;   % SI
kappa_0 = 1.E-6;
mu_a  = muscale;

%from MultipleRuns.py or codes like it in the folder(s) above, we 
%establish the Temp scale
%temp_values = [27.+273, Tb+273, 1300.+273, 1500.+273]
%dTemp = temp_values[3] - temp_values[0]
Tscale = 1305-27;
Tval   = 1270+273;
tcont  = [300:100:1600];

h      = 4e2; % box dimension in km

hscale = 4e5; % box scale in m
tcont  = [300:100:1600]; %in K
pscale = 1192135725.0; % pressure scale from MultipleRuns.py in Pa
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
    loc  = [locroot mstrc Tbstrc '/temp']
    loc2 = [locroot mstrc Tbstrc '/velocity']
    loc3 = [locroot mstrc Tbstrc '/gradp']
    loc4 = [locroot mstrc Tbstrc '/mu']
    loc_cd = [locroot mstrc Tbstrc ]
    %loc2 = ['HK03/' mstr Tbstr '/pstar']
    
    %if long==1
    %if we have long runs, uncomment the following for loop
    %  if outind <= 800 & outind >= 100
        ii = outind*2;
        clear dat val profile
        name=[loc '_' num2str(ii) '.csv'];
        dat = csvread(name,1,0);
        x = dat(:,2);
        y = dat(:,3);
        val = dat(:,1);

        clear dat
        name=[loc2 '_' num2str(ii) '.csv'];
        dat = csvread(name,1,0);
        vx  = dat(:,1);
        vy  = dat(:,2);

        clear dat
        name=[loc3 '_' num2str(ii) '.csv'];
        dat = csvread(name,1,0);
        dpdx  = dat(:,1);
        dpdy  = dat(:,2);

        clear dat
        name=[loc4 '_' num2str(ii) '.csv'];
        dat = csvread(name,1,0);
        mu = dat(:,1);
        
        uy     = unique(y);
        ux     = unique(x);
%         bigmat = [x y val vx vy dpdx dpdy mu];
%         newmat = sortrows(bigmat, 2);
%         for i=1:size(uy);
%             rrs    = newmat(:,2)==uy(i);
%             choose = rrs*ones(1,8);
%             subarr = newmat.*choose;
%             extracted = subarr(find(rrs),:);
%             sorted = sortrows(extracted,1);
%             if i==1 
%                 stack = sorted; 
%             else 
%                 stack = [stack;sorted];
%             end
%         end
%         
%         x = stack(:,1);
%         y = stack(:,2);
%         val = stack(:,3);
%         vx = stack(:,4);
%         vy = stack(:,5);
%         dpdx = stack(:,6);
%         dpdy = stack(:,7);
%         mu  = stack(:,8);
%         
        
%             nn   = length(x);
%             maxx = max(x);
%             minx = min(x);
%             delx   = diff(x);
%             bigchange = minx - maxx;
%             ind1 = 0;
%             xind(1) = 1;
%             yind(1) = 1;

    %         ncol = find(delx == bigchange, 1);
    %         nrow = sum((diff(y) ~= 0))+1;
    %         dx   = h*(delx(1));
    %         dy   = h*(y(ncol+1) - y(1));
    %         xvec = [0:dx:max(x)];
    %         yvec = [0:dy:max(y)];

%             for i=2:nn-1
%                 xind(i) = xind(i-1) + 1;
%                 yind(i) = yind(i-1);
%                 if delx(i-1) == bigchange
%                     xind(i) = 1;
%                     yind(i) = yind(i-1)+1;
%                 end
%             end
%             xind(nn) = xind(nn-1) + 1;
%             yind(nn) = yind(nn-1);

        for i=1:size(uy)
             for j=1:size(ux)
                r1 = (x==ux(j));
                r2 = (y==uy(i));
                rp = r1.*r2;
                inds = find(rp);
                ind  = inds(1);
                T(i,j)  = val(ind); 
                MU(i,j) = mu(ind);
                X(i,j)  = x(ind);
                Y(i,j)  = y(ind); 
                U(i,j)  = vx(ind); 
                V(i,j)  = vy(ind); 
                DPDX(i,j) = dpdx(ind); 
                DPDY(i,j) = dpdy(ind); 
             end
        end


    %scale to real dimensional values
        Xr = 1000*X;  % in km
        Yr = 1000*Y;
        Tr = 273+1300*(T-.2347);
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
        Vmelty = WY + V;
        
     %output
        figure(1);clf
        %[cs h] = 
        %subplot('Position',[0.1 0.25 0.8 0.7])
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
        %view(2);  shading interp
        %clabel(cs, h, tcont,'fontsize',[12]);

        
        figure(1)% now overlay streamlines and velocity vectors
        hold on
        plot(startx,starty,'wo'); hold on
        han = streamline(Xr,Yr,Vmeltx,Vmelty,startx,starty); 
        %han = streamline(Xr,Yr,WX,WY,startx,starty); 
        set(han,'color','r','linewidth',[1.25]);
        %quiver(Xr,Yr,U,V,3,'k'); 
        quiver(Xr(1:2:end,1:2:end),Yr(1:2:end,1:2:end),U(1:2:end,1:2:end),V(1:2:end,1:2:end),3,'k'); 

        %title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14],'xlim',[0 1000],'ylim',[0 400])
        ylabel('km'); 
        xlabel('km')
        box on
        ax = gca;
        ax.Clipping = 'off';
        %colorbar
        hold off
         if flag == 1
            filename = ['Streamplot_t_' num2str(ii) ];
            WD1 = cd;
            cd(loc_cd)
            figure(1); title('');
            %chan=colorbar
            %set(chan,'ydir','reverse')
            brighten(1,.05)
            caxis([27 1350])
            eval(['print -dpdf ' filename])
            cd(WD1)
        end
        %pause(0.5);      
        %find pressure-gradients along average isotherm depth, zll
%         close = min(min(abs(Yr-zll)));
%         xpos  = Xr(abs(Yr-zll)==close);
%         p_x   = DPDX(abs(Yr-zll)==close);
%         figure(1);
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
            
%             disp('write output?')
%             pause
%             WD1 = cd;
%             cd(loc_cd);
%                filename = ['Streamplot_t_' num2str(ii) ];
%                figure(1); title('')
%                eval(['print -dpng ' filename])
%                %print -dpdf with_convection_398my_melt.pdf
%                %savefig('with_convection_398my_melt.fig')
% %                figure(4)
% %                print -dpdf with_convection_1164my_dpdx.pdf
% %                savefig('with_convection_1164my_dpdx.fig')
%             cd(WD1)
      %end
    %end
    %%
    
  