% script to plot profiles of dynamic pstar and gradp in Matlab

% Step 1: using output written in .vtu from Fenics code
% utilizing bash script by VTKtoMatlab.sh by S. Gold
% utilizing bash script by VTKtoMatlab_vec.sh by M. Roy


% Mousumi Roy, Dec 6, 2013

clear all
% output interval
%
int = 10; %output interval
lw  = 1;  %default linewidth
prange = [1.1, 1.5]*1.e9;
gradprange = [-1000 1000];
flag = 0; % flag to save/print figure

% first choose folder where data are:

Tb   = 1300;
mstr = ['mu\=1e+19/'];
Tbstr= ['Tb\=' num2str(Tb)];
loc  = ['BLNK89/b12.7clog128/' mstr Tbstr '/pstar']
loc2 = ['HK03/' mstr Tbstr '/pstar']

pscale = 1192135725.0;
h      = 1e6; % box dimension in m

dep    = 50; %in km
depstr = ['_' num2str(dep) 'km_'];
depf   = (h - dep*1e3)/h;

filename = ['Tb' num2str(Tb) 'mu' mstr(5:9) 'dep' num2str(dep) 'km']


%vertical profile at center
posx = 0.5*ones(100,1);
posy = [0.01:0.01:1]';

%horizontal profile at 200 km depth
posy = depf*ones(100,1);
posx = [0.01:0.01:1]';


figure(1);clf

for ii = 1:int:9
    clear dat val profile
    eval(['! sh VTKtoMatlab.sh ' loc '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy);
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    subplot(221)
    plot(posx*1e3,profile,'r-','linewidth',lw);hold on
    title([loc(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(222)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
    %pause(0.5)
    
end

for ii = 10:int:99
    clear dat val profile
    eval(['! sh VTKtoMatlab.sh ' loc '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy);
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    
    subplot(221)
    plot(posx*1e3,profile,'r-','linewidth',lw); hold on
    title([loc(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(222)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
    %pause(0.5)
end
plot(posx*1e3,profile,'k-','linewidth',[1.5]); hold on
set(gca,'fontname','Helvetica','fontsize',[14],'ylim',prange)
xlabel('km')

for ii = 1:int:9
    clear dat val profile
    eval(['! sh VTKtoMatlab.sh ' loc2 '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy);
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    
    subplot(223)
    plot(posx*1e3,profile,'b-','linewidth',lw);hold on
    title([loc2(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(224)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
%     
    %pause(0.5)
    
end

for ii = 10:int:99
    clear dat val profile
    eval(['! sh VTKtoMatlab.sh ' loc2 '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy);
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    
    subplot(223)
    plot(posx*1e3,profile,'b-','linewidth',lw); hold on
    title([loc2(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(224)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
%     %pause(0.5)
end
%
plot(posx*1e3,profile,'k-','linewidth',[1.5]); hold on
set(gca,'fontname','Helvetica','fontsize',[14],'ylim',prange)
xlabel('km')

% --------------- Now repeat for gradp ---------------------
%
loc  = ['BLNK89/b12.7clog128/' mstr Tbstr '/gradp']
loc2 = ['HK03/' mstr Tbstr '/gradp']

for ii = 1:int:9
    clear dat val profile
    eval(['! sh VTKtoMatlab_vec.sh ' loc '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = sqrt(dat(:,3).*dat(:,3) + dat(:,4).*dat(:,4));
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy)/h;
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    
    subplot(222)
    plot(posx*1e3,profile,'r-','linewidth',lw);hold on
    title([loc(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(222)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
    %pause(0.5)
    
end

for ii = 10:int:99
    clear dat val profile
    eval(['! sh VTKtoMatlab_vec.sh ' loc '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = sqrt(dat(:,3).*dat(:,3) + dat(:,4).*dat(:,4));
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy)/h;
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    
    subplot(222)
    plot(posx*1e3,profile,'r-','linewidth',lw); hold on
    title([loc(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(222)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
    %pause(0.5)
end
plot(posx*1e3,profile,'k-','linewidth',[1.5]); hold on
set(gca,'fontname','Helvetica','fontsize',[14],'ylim',gradprange)
xlabel('km')

for ii = 1:int:9
    clear dat val profile
    eval(['! sh VTKtoMatlab_vec.sh ' loc2 '00000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = sqrt(dat(:,3).*dat(:,3) + dat(:,4).*dat(:,4));
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy)/h;
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    
    subplot(224)
    plot(posx*1e3,profile,'b-','linewidth',lw);hold on
    title([loc2(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(224)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
%     
    %pause(0.5)
    
end

for ii = 10:int:99
    clear dat val profile
    eval(['! sh VTKtoMatlab_vec.sh ' loc2 '0000' num2str(ii) '.vtu'])
    dat = load('PythonSoln');
    x = dat(:,1);
    y = dat(:,2);
    val = sqrt(dat(:,3).*dat(:,3) + dat(:,4).*dat(:,4));
    val = dat(:,3);
    F = TriScatteredInterp(x,y,val);
    profile = pscale*F(posx,posy)/h;
    
%     subplot(211)
%     plot3(x,y,val,'o'); 
    
    subplot(224)
    plot(posx*1e3,profile,'b-','linewidth',lw); hold on
    title([loc2(1:4) depstr ', t = ' num2str(ii) ' m.y.']);
%     subplot(224)
%     plot(posx*1e3,gradient(profile/h),'-'); hold on
%     title([mstr depstr num2str(Tb) ' t = ' num2str(ii) ' m.y.']);
%     %pause(0.5)
end
plot(posx*1e3,profile,'k-','linewidth',[1.5]); hold on
set(gca,'fontname','Helvetica','fontsize',[14],'ylim',gradprange)
xlabel('km')


if flag == 1
    eval(['print -dpng ' filename])
end

%---------------------------------------------------------------
% OLD STUFF - DO NOT ERASE
% using output written in .csv from Paraview (after analysis_step.txt)
% utilizing bash script by VTKtoMatlab.sh by S. Gold

% Mousumi Roy, Dec 3, 2013

% first choose folder where data are:
% clear all
% 
% figure(1);clf
% 
% mstr = ['mu\=1e+19/'];
% depstr = ['_200km_'];
% Tbstr= ['Tb\=800'];
% loc  = ['BLNK89/b12.7clog128/' mstr Tbstr '/Pstar' depstr 'depth.']
% loc2 = ['HK03/' mstr Tbstr '/Pstar' depstr 'depth.']
% 
% for ii = 1:100
%     clear dat val profile
%     eval(['! sh readscalars.sh ' loc num2str(ii-1) '.csv'])
%     dat = load('datarr');
%     x = dat(:,2);
%     y = dat(:,3);
%     val = dat(:,1);
%     %F = TriScatteredInterp(x,y,val)
%     %profile = F(posx,posy);
%     
%     subplot(231)
%     %plot3(x,y,val,'o'); 
%     title(['t = ' num2str(ii) ' m.y.'])
%     %subplot(223)
%     plot(x,val,'-');hold on
%     %subplot(224)
%     %plot(posy,gradient(profile),'o')
%     
%     %
%     
% end
% 
% pause(0.5)
% 
% 
% for ii = 1:100
%     clear dat val profile
%     eval(['! sh readscalars.sh ' loc2 num2str(ii-1) '.csv'])
%     dat = load('datarr');
%     x = dat(:,2);
%     y = dat(:,3);
%     val = dat(:,1);
%     %F = TriScatteredInterp(x,y,val)
%     %profile = F(posx,posy);
%     
%     subplot(234)
%     %plot3(x,y,val,'o'); 
%     title(['t = ' num2str(ii) ' m.y.'])
%     %subplot(223)
%     plot(x,val,'-');hold on
%     %subplot(224)
%     %plot(posy,gradient(profile),'o')
%     
%     %pause(0.5)
%     
% end
% 
% pause(0.5)
% 
% Tbstr= ['Tb\=1000'];
% loc  = ['BLNK89/b12.7clog128/' mstr Tbstr '/Pstar' depstr 'depth.']
% loc2 = ['HK03/' mstr Tbstr '/Pstar' depstr 'depth.']
% 
% for ii = 1:100
%     clear dat val profile
%     eval(['! sh readscalars.sh ' loc num2str(ii-1) '.csv'])
%     dat = load('datarr');
%     x = dat(:,2);
%     y = dat(:,3);
%     val = dat(:,1);
%     %F = TriScatteredInterp(x,y,val)
%     %profile = F(posx,posy);
%     
%     subplot(232)
%     %plot3(x,y,val,'o'); 
%     title(['t = ' num2str(ii) ' m.y.'])
%     %subplot(223)
%     plot(x,val,'-');hold on
%     %subplot(224)
%     %plot(posy,gradient(profile),'o')
%     
%     %pause(0.5)
%     
% end
% 
% pause(0.5)
% 
% for ii = 1:100
%     clear dat val profile
%     eval(['! sh readscalars.sh ' loc2 num2str(ii-1) '.csv'])
%     dat = load('datarr');
%     x = dat(:,2);
%     y = dat(:,3);
%     val = dat(:,1);
%     %F = TriScatteredInterp(x,y,val)
%     %profile = F(posx,posy);
%     
%     subplot(235)
%     %plot3(x,y,val,'o'); 
%     title(['t = ' num2str(ii) ' m.y.'])
%     %subplot(223)
%     plot(x,val,'-');hold on
%     %subplot(224)
%     %plot(posy,gradient(profile),'o')
%     
%     %pause(0.5)
%     
% end
% 
% pause(0.5)
% 
% Tbstr= ['Tb\=1000'];
% loc  = ['BLNK89/b12.7clog128/' mstr Tbstr '/Pstar' depstr 'depth.']
% loc2 = ['HK03/' mstr Tbstr '/Pstar' depstr 'depth.']
% 
% for ii = 1:100
%     clear dat val profile
%     eval(['! sh readscalars.sh ' loc num2str(ii-1) '.csv'])
%     dat = load('datarr');
%     x = dat(:,2);
%     y = dat(:,3);
%     val = dat(:,1);
%     %F = TriScatteredInterp(x,y,val)
%     %profile = F(posx,posy);
%     
%     subplot(233)
%     %plot3(x,y,val,'o'); 
%     title(['t = ' num2str(ii) ' m.y.'])
%     %subplot(223)
%     plot(x,val,'-');hold on
%     %subplot(224)
%     %plot(posy,gradient(profile),'o')
%     
%     %pause(0.5)
%     
% end
% 
% pause(0.5)
% 
% for ii = 1:100
%     clear dat val profile
%     eval(['! sh readscalars.sh ' loc2 num2str(ii-1) '.csv'])
%     dat = load('datarr');
%     x = dat(:,2);
%     y = dat(:,3);
%     val = dat(:,1);
%     %F = TriScatteredInterp(x,y,val)
%     %profile = F(posx,posy);
%     
%     subplot(236)
%     %plot3(x,y,val,'o'); 
%     title(['t = ' num2str(ii) ' m.y.'])
%     %subplot(223)
%     plot(x,val,'-');hold on
%     %subplot(224)
%     %plot(posy,gradient(profile),'o')
%     
%     %pause(0.5)
%     
% end
% 
% 
% pause(0.5)