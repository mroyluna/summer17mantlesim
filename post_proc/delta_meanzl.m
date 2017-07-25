% script to analyze the variation of mean LAB depth in Matlab
% Mousumi Roy, July 17, 2014

clear all
flag = 1; % flag to save/print figures

figure (1); clf;
colormap(gray)
col = ['k','r','g','b','m'];
sym  = ['o','d','s'];

lims = [-20 100];
%
int = 5; %output interval
lw  = 1;  %default linewidth
%prange = [1.1, 1.5]*1.e9;
%gradprange = [-1000 1000];

%first, plot all results for meanzl in a manner where they can all be
%compared, namely, plot change in meanzl from starting value as a function 
% of time

% first choose folder where data are
locroot = ['test_steady_state/BLNK_LAB.85/']
Tbs = [800 1000 1300];
%Tbs = [1300];
k=1;

for Tbcount = 1:length(Tbs)
    Tb   = Tbs(Tbcount);
    Tbstr= ['Tb=' num2str(Tb)];
    
    subplot(151)
    mstr = ['mu=1e+23/'];
    name = [locroot mstr Tbstr '/meanz_1227.dat'];
    clear dat
    dat = load(name);
    
    t  = dat(:,1); % time in my
    zl = 1000-dat(:,2); % mean 1227 C isotherm depth in km
 %output
    %plot(t,zl,'marker',sym(Tbcount),'color',col(k)); hold on
    plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
    title(mstr); xlabel('Time, my');
    ylabel('Change in LAB Depth, km')
    set(gca,'ydir','reverse','ylim',lims,'fontname','Helvetica','fontsize',[14])
    k = k+1;
    
    subplot(152)
    mstr = ['mu=1e+22/'];
    name = [locroot mstr Tbstr '/meanz_1227.dat'];
    clear dat
    dat = load(name);
    
    t  = dat(:,1); % time in my
    zl = 1000-dat(:,2); % mean 1227 C isotherm depth in km
 %output
    %plot(t,zl,'marker',sym(Tbcount),'color',col(k)); hold on
    plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
    title(mstr); xlabel('Time, my');
    set(gca,'ydir','reverse','ylim',lims,'fontname','Helvetica','fontsize',[14])
    k = k+1;
    
    subplot(153)
    mstr = ['mu=1e+21/'];
    name = [locroot mstr Tbstr '/meanz_1227.dat'];
    clear dat
    dat = load(name);
    
    t  = dat(:,1); % time in my
    zl = 1000-dat(:,2); % mean 1227 C isotherm depth in km
 %output
    %plot(t,zl,'marker',sym(Tbcount),'color',col(k)); hold on
    plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
    title(mstr); xlabel('Time, my');
    set(gca,'ydir','reverse','ylim',lims,'fontname','Helvetica','fontsize',[14])
    k = k+1;
    
    subplot(154)
    mstr = ['mu=1e+20/'];
    name = [locroot mstr Tbstr '/meanz_1227.dat'];
    clear dat
    dat = load(name);
    
    t  = dat(:,1); % time in my
    zl = 1000-dat(:,2); % mean 1227 C isotherm depth in km
 %output
    %plot(t,zl,'marker',sym(Tbcount),'color',col(k)); hold on
    plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
    title(mstr); xlabel('Time, my');
    set(gca,'ydir','reverse','ylim',lims,'fontname','Helvetica','fontsize',[14])
    k = k+1;
    
    subplot(155)
    mstr = ['mu=1e+19/'];
    name = [locroot mstr Tbstr '/meanz_1227.dat'];
    clear dat
    dat = load(name);
    
    t  = dat(:,1); % time in my
    zl = 1000-dat(:,2); % mean 1227 C isotherm depth in km
 %output
    %plot(t,zl,'marker',sym(Tbcount),'color',col(k)); hold on
    plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
    title(mstr); xlabel('Time, my');
    set(gca,'ydir','reverse','ylim',lims,'fontname','Helvetica','fontsize',[14])
    
    subplot(151)
    %plot(t,zl,'marker',sym(Tbcount),'color',col(k)); hold on
    plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
    title(mstr); xlabel('Time, my');
    set(gca,'ydir','reverse','ylim',lims,'fontname','Helvetica','fontsize',[14])
    k = 1;
    
    %pause(5)
    
end
subplot(152)
legend('Tb=800','Tb=1000','Tb=1300')
%legend('Tb=1300')
    
if flag == 1
    WD = cd;
    cd(locroot);
    filename = ['DelMeanZl_t.pdf'];
    eval(['print -dpdf ' filename])
    cd(WD)
end
      

