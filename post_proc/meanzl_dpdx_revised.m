% script to analyze the variation of mean LAB depth in Matlab

% Mousumi Roy, July 17, 2014
% modified Oct 10, 14 -- to look at evolution of dpdx through time

clear all
flag = 0; % flag to save/print figures

col = ['k','r','g','b','m'];
sym  = ['o','d','s'];

lims = [0 500];
%
int = 5; %output interval
lw  = 3;  %default linewidth
lc  = 'g'; %line color
%prange = [1.1, 1.5]*1.e9;
%gradprange = [-1000 1000];


%first, plot all results for meanzl in a manner where they can all be
%compared, namely, plot change in meanzl from starting value as a function 
% of time


% first choose folder where data are
%locroot = ['test_steady_state/BLNK_LAB.5/longrun/']
%locroot = ['longruns_no_adiabat/test2/']
locroot = ['longruns_no_adiabat/test2/']

%Tbs = [800 1000 1300];
Tbs = [1300];
k=1;

mstr = ['mu=1e+18/'];

rho = 3300;
g   = 9.81;
rg  = rho*g;

Tbcount = 1;
%for Tbcount = 1:length(Tbs)
    Tb   = Tbs(Tbcount);
    Tbstr= ['Tb=' num2str(Tb)];
    
%%
ii=1;

for ii = 1:490
    %name = [locroot mstr Tbstr '/meanz_1227.dat'];
    name = [locroot mstr Tbstr '/contour_989.' num2str(ii) '.csv'];
    %name = [locroot mstr Tbstr '/contour_0.989.' num2str(ii) '.csv'];
    clear dat
    dat = csvread(name,1,0);

    contx = dat(:,1);
    conty = 400*dat(:,2)/0.4;
    meany = mean(conty);
    t(ii)  = 2*(ii-1);      % time in my
    zl(ii) = 400-meany; % mean 1227 C isotherm depth in km
    
%     Ra_int = dat(1:end-1,3); %interior Rayleigh number
%     
%     
%     k = 1;
%         name = [locroot mstr Tbstr '/dpdx_zll_t_' num2str(1)]; 
%         clear dat
%         dat  = load(name) ;
%         xpos = dat(:,1);
%         dpdx = dat(:,2);
%         maxdpdx(k,1) = max(dpdx); 
%         mdpdx(k,1)   = mean(abs(dpdx));
%         k=k+1;
%         %subplot(211)
%         %plot(xpos,dpdx,'--');hold on; pause(1)
%         name = [locroot mstr Tbstr '/dpdx_zll_t_' num2str(6)];
%         clear dat
%         dat  = load(name)
%         xpos = dat(:,1);
%         dpdx = dat(:,2);
%         maxdpdx(k,1) = max(dpdx);
%         mdpdx(k,1)   = mean(abs(dpdx));k=k+1;
        %subplot(211)
        %plot(xpos,dpdx,'--');hold on;pause(1)
    %for ii=10:5:799
%     for ii=2:1:799
%         name = [locroot mstr Tbstr '/dpdx_zll_t_' num2str(ii)];
%         clear dat
%         dat  = load(name);
%         xpos = dat(:,1);
%         dpdx = dat(:,2);
%         maxdpdx(k,1) = max(dpdx);
%         mdpdx(k,1)   = mean(abs(dpdx));k=k+1;
%         %plot(xpos,dpdx,'--');hold on;pause(1)
%     end
    
end

    % find moving averages
    %mmRa_int = movingmean(Ra_int,5);
    %mmdpdx   = movingmean(mdpdx,31); %or for long runs, I used window of 31
 %output
    figure(4);% clf
    %clf;
    %subplot(131); 
    plot(t,zl,'k-','linewidth',lw, 'color', lc); hold on
    %plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
    %legend('Tb=800','Tb=1000','Tb=1300');
    %legend('Tb=1300');
    %title(mstr); xlabel('Time, my');
    
    ylabel('LAB Depth, km')
    xlabel('Time, my')
    set(gca,'ydir','reverse','ylim',lims,'fontname','Helvetica','fontsize',[14])
    grid on; box on
    set(gca,'xlim',[0 max(t)],'ylim',[0 400],'ydir','reverse')
    plot([0 max(t)], [188 188],'r--','linewidth',[1])
    %plot([0 max(t)], [135 135],'r--','linewidth',[1])
    
    %figure(5); clf
%     subplot(132); hold on
%     plot(t,mmRa_int,'k-','linewidth',lw, 'color', lc); hold on
%     %plot(t,zl - zl(1),'marker',sym(Tbcount),'color',col(k)); hold on
%     %legend('Tb=800','Tb=1000','Tb=1300');
%     %legend('Tb=1300');
%     %title(mstr); xlabel('Time, my');
%     ylabel('Ra_i')
%     xlabel('Time, my')
%     set(gca,'ydir','reverse','fontname','Helvetica','fontsize',[14])
%     grid on; box on
%     plot([0 1800], [2.5 2.5]*1e7,'r--')
%     set(gca,'xlim',[0 1800])
    
    
    %figure(6); clf;
%     subplot(233); hold on
%     plot(t,mdpdx/(rg),'k-','linewidth',lw, 'color', lc); hold on
%     ylabel('Normalized dp/dx')
%     xlabel('Time, my');set(gca,'xlim',[0 1800])
%     set(gca,'fontname','Helvetica','fontsize',[14])
%     grid on; box on
%end

if flag == 1
  
    figure(4);
%     subplot(236);
%     text(19.5,555,'220','fontname','HelveticaBold','color','r','fontsize',[10],'fontweight','bold');hold on
%     plot(20,500,'.','markersize',[32],'color','r');
%     plot(20,300,'.','markersize',[32],'color','r');
%     plot(20,0,'.','markersize',[32],'color','r');
%     plot(19,500,'.','markersize',[32],'color','r');
%     plot(19,300,'.','markersize',[32],'color','r');
%     plot(19,0,'.','markersize',[32],'color','r');
%     plot(21,0,'.','markersize',[32],'color','k');
%     plot(21,300,'.','markersize',[32],'color','k');
%     plot(21,500,'.','markersize',[32],'color','k');
%     plot(22,0,'.','markersize',[32],'color','k');
%     plot(22,300,'.','markersize',[32],'color','k');
%     plot(22,500,'.','markersize',[32],'color','k');
%     plot(23,0,'.','markersize',[32],'color','k');
%     plot(23,300,'.','markersize',[32],'color','k');
%     plot(23,500,'.','markersize',[32],'color','k');
%     text(19.5,355,'224','fontname','HelveticaBold','color','r','fontsize',[10]);
%     text(19.5,55,'232','fontname','HelveticaBold','color','r','fontsize',[10]);
%     text(18.5,555,'18','fontname','HelveticaBold','color','r','fontsize',[10]);hold on
%     text(18.5,355,'23','fontname','HelveticaBold','color','r','fontsize',[10]);
%     text(18.5,55,'29','fontname','HelveticaBold','color','r','fontsize',[10]);
%     set(gca,'xlim',[18 23],'ylim',[-100 600],'ytick',[0 300 500],'xtick',[19 20 21 22]);
%     plot([20.5 20.5],ylim,'k--');
%     ylabel('\Delta T_{b}, C')
%     xlabel('Log10(\eta_{0} in Pa s)')
%     set(gca,'fontname','Helvetica','fontsize',[14]); box on    
%     figure(5); clf
%     subplot(131); hold on
%     plot(t,mdpdx,'k-','linewidth',lw, 'color', lc); hold on
%     
%     ylabel('Dpdx, Pa/m')
%     xlabel('Time, my');set(gca,'xlim',[0 100])
%     set(gca,'fontname','Helvetica','fontsize',[14])
%     grid on; box on
    
    

    %subplot(131)
    %line = [20 105; 20 230];
    %plot(line(:,1),line(:,2),'k--','linewidth',[0.5]); hold on
    line = [62 100; 62 250];
    plot(line(:,1),line(:,2),'k--','linewidth',[0.5])
    line = [398 100; 398 250];
    plot(line(:,1),line(:,2),'k--','linewidth',[0.5])
    line = [760 100; 760 250];
    plot(line(:,1),line(:,2),'k--','linewidth',[0.5])
% %     line = [1076 175; 1076 230];
%     plot(line(:,1),line(:,2),'k--','linewidth',[0.5])
%     line = [1164 175; 1164 230];
%     plot(line(:,1),line(:,2),'k--','linewidth',[0.5])
%    % gtext('3a');gtext('3b');gtext('3c');gtext('3d');gtext('3e');gtext('3f');
    
    WD = cd;
    cd(locroot);
    cd ..;
    filename = ['compare_longruns_LAB_t_rev.pdf'];
    eval(['print -dpdf ' filename])
    filename = ['compare_longruns_LAB_t_rev.eps'];
    eval(['print -depsc ' filename])
    figure(4);%
    savefig compare_longruns_LAB_t_rev.fig
%     figure(5);
%     filename = ['dpdx_through_t.pdf'];
%     eval(['print -dpdf ' filename])
    %figure(6);
    %filename = ['compare_longruns_Dpdx_t.pdf'];
    %eval(['print -dpdf ' filename])
    cd(WD)
end
   

