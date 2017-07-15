% script to plot color plots of T, with velocity arrows, and streamlines in Matlab

% Step 1: using output written in .vtu from Fenics code
% utilizing bash script by VTKtoMatlab.sh by S. Gold
% utilizing bash script by VTKtoMatlab_vec.sh by M. Roy

% Mousumi Roy, Feb 7, 2014

clear all

colormap(gray)
% output interval
%
int = 1; %output interval
lw  = 1;  %default linewidth
%prange = [1.1, 1.5]*1.e9;
%gradprange = [-1000 1000];

flag = 0; % flag to save/print figures
long = 1; % if longer run, then do last for-loop
timeint = 2; %if writing output every 2 my instead of 1 my, need this
%Tbs = [800 1000 1300];
Tbs = [1300];
% first choose folder where data are:%
locroot = ['test_steady_state/BLNK_LAB.9/longrun/']
%locroot = ['test_steady_state/BLNK_LAB.9/shortbox/']
%locroot = ['test_steady_state/BLNK_LAB.7/shortbox/']
%locroot = ['tanhstep/BLNK/w0.2h0.1sc0.05/shortbox/']
mstr = ['mu\=1e+19/'];
mstrc= ['mu=1e+19/'];

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
Tval   = 1500;
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
drho    = 500;%[drhomin:drhomax];
kovermu = 1e-15/1e-2;
startx  = [1:20:990]'; %note here units must be in km as displayed in box
starty  = 750*ones(length(startx),1);

figure(2);clf;
figure(1);clf;colormap(gray)
numt = 1;

prevDist = []
    
for Tbcount = 1:length(Tbs)
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
        DPDY = pscale*DPDY/hscale;
        DPDX = pscale*DPDX/hscale;
        WY = -(DPDY - drho*g)*kovermu;
        WX = -(DPDX)*kovermu;
        MU = muscale*MU;
        
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
        han = streamline(Xr,Yr,WX+U*vscale,WY+V*vscale,startx,starty); 
        prevDist = trackStream(prevDist, han, 20,min(min(Xr)),max(max(Xr)));
        
        
        set(han,'color','r','linewidth',[1.25]);
        quiver(Xr,Yr,U,V,3,'k'); 
        %pause
        
        title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14])
        xlabel('km');ylabel('km'); 
        box on
        %colorbar
        hold off
        
        if flag == 1
            filename = ['Streamplot_t_' num2str(ii) ];
            WD1 = cd;
            cd(loc_cd)
            eval(['print -dpdf ' filename])
            cd(WD1)
        end
        
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
            filename = ['dpdx_zll_t_' num2str(ii)];
            dat = [xpos p_x];
            WD1 = cd;
            cd(loc_cd)
            eval(['save ' filename ' dat -ascii'])
            cd(WD1)
            
        pause(0.5);  %pause
        
       
    end
    %pause
    %%
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
        
        clear dat
        eval(['! sh VTKtoMatlab.sh ' loc4 '0000' num2str(ii) '.vtu'])
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
        WY = -(DPDY - drho*g)*kovermu;
        WX = -(DPDX)*kovermu;
        MU = muscale*MU;
        
     %output
        figure(1);clf
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
          numt           = numt + 1; 
          
        %view(2);  shading interp
        %clabel(cs, h, tcont,'fontsize',[12]);

        figure(1); hold on
        plot(startx,starty,'wo'); 
        han = streamline(Xr,Yr,WX+U*vscale,WY+V*vscale,startx,starty); 
        set(han,'color','r','linewidth',[1.25]);
        quiver(Xr,Yr,U,V,3,'k'); 

        title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14])
        xlabel('km');ylabel('km'); 
        box on
        %colorbar
        hold off
       
         if flag == 1
            filename = ['Streamplot_t_' num2str(ii) ];
            WD1 = cd;
            cd(loc_cd)
            eval(['print -dpdf ' filename])
            cd(WD1)
        end
        %find pressure-gradients along average isotherm depth, zll
        close = min(min(abs(Yr-zll)));
        xpos  = Xr(abs(Yr-zll)==close);
        p_x   = DPDX(abs(Yr-zll)==close);
        figure(4);clf;subplot(411)
        plot(xpos,p_x,'k-');
        set(gca,'fontname','Helvetica','fontsize',[14],'ylim',[-150 150])
        xlabel('km'); ylabel('Pa/m')
         %output to file
            filename = ['dpdx_zll_t_' num2str(ii)];
            dat = [xpos p_x];
            WD1 = cd;
            cd(loc_cd)
            eval(['save ' filename ' dat -ascii'])
            cd(WD1)
        
        pause(0.5);  
        
        %make a figure for paper at a chosen time...
        if ii == 0.0001
        %if ii == 10
            colormap(gray)
            figure(3)
            %[cs h] = 
            contourf(Xr,Yr,Tr,tcont);hold on
            contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
            quiver(Xr,Yr,U,V,3,'k'); 
            title(' ')
            ylabel('km')
            xlabel('km')
            set(gca,'fontname','Helvetica','fontsize',[14]); 
            box on
            
            pause
            WD1 = cd;
            cd(loc_cd)
               print -dpdf no_convection_20my.pdf
               savefig('no_convection_20my.fig')
               figure(1); title('')
               print -dpdf no_convection_20my_melt.pdf
               savefig('no_convection_20my_melt.fig')
               figure(4)
               print -dpdf no_convection_20my_dpdx.pdf
               savefig('no_convection_20my_dpdx.fig')
            cd(WD1)
            
        end
        if ii == 0.0001
        %if ii == 30
            colormap(gray)
            figure(3)
            %[cs h] = 
            contourf(Xr,Yr,Tr,tcont);hold on
            contour(Xr,Yr,Tr, [Tval, Tval],'k','linewidth',[2]);
            quiver(Xr,Yr,U,V,3,'k'); 
            title(' ')
            ylabel('km')
            xlabel('km')
            set(gca,'fontname','Helvetica','fontsize',[14]); 
            box on
            pause
            WD1 = cd;
            cd(loc_cd)
               print -dpdf with_convection_60my.pdf
               savefig('with_convection_60my.fig')
               figure(1); title('')
               print -dpdf with_convection_60my_melt.pdf
               savefig('with_convection_60my_melt.fig')
               figure(4)
               print -dpdf with_convection_60my_dpdx.pdf
               savefig('with_convection_60my_dpdx.fig')
            cd(WD1)
            
        end
    end
    %% comment here 
    
    if long==1
    %if we have long runs, uncomment the following for loop
      for ii = 100:int:800
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
        WY = -(DPDY - drho*g)*kovermu;
        WX = -(DPDX)*kovermu;
        MU = muscale*MU;
        
     %output
        figure(1);clf
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
        %view(2);  shading interp
        %clabel(cs, h, tcont,'fontsize',[12]);

        figure(1);hold on
        plot(startx,starty,'wo'); 
        han = streamline(Xr,Yr,WX+U*vscale,WY+V*vscale,startx,starty); 
        set(han,'color','r','linewidth',[1.25]);
        quiver(Xr,Yr,U,V,3,'k'); 

        title([loc(1:4) ', t = ' num2str(timeint*ii) ' m.y.']);
        set(gca,'fontname','Helvetica','fontsize',[14])
        xlabel('km'); ylabel('km'); 
        box on
        %colorbar
        hold off

         if flag == 1
            filename = ['Streamplot_t_' num2str(ii) ];
            WD1 = cd;
            cd(loc_cd)
            eval(['print -dpdf ' filename])
            cd(WD1)
        end
        pause(0.5);      
        %find pressure-gradients along average isotherm depth, zll
        close = min(min(abs(Yr-zll)));
        xpos  = Xr(abs(Yr-zll)==close);
        p_x   = DPDX(abs(Yr-zll)==close);
        figure(4);clf;subplot(411)
        plot(xpos,p_x,'k-');
        set(gca,'fontname','Helvetica','fontsize',[14],'ylim',[-150 150])
        xlabel('km'); ylabel('Pa/m')
        %output to file
            filename = ['dpdx_zll_t_' num2str(ii)];
            dat = [xpos p_x];
            WD1 = cd;
            cd(loc_cd)
            eval(['save ' filename ' dat -ascii'])
            cd(WD1)
        figure(1)
      end
    end
    
    
    WD1 = cd;
    cd(loc_cd)
      filename = ['meanz_' num2str(Tval-273) '.dat'];
      save(filename, 'meanz', '-ascii')
    cd(WD1)
end
%%