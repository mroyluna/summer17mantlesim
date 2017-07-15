% Code to compare and synthesize dynamic topography results
% that are outputted as .fig files by finduplift.m
% Mousumi Roy
% Feb 24, 2015

% root folder 
clear all
figure(1);clf
figure(2);clf

locroot = ['tanhstep_smallbox_400km_no_adiabat/']

%locroot = ['tanhstep/BLNK/w0.2h0.05sc0.05/'];
ii = 20 ; % time at which comparison is being made
mi = 1;
clf
for Tb = [1300]
    for muscale = [1.e+18 5.e+18 1.e+21 5.e+21]
        inp   = ['1e+18'; '5e+18'; '1e+21'; '5e+21'];
        ms   = [inp(mi,:)]; 
        mstr = ['mu=' ms ];
        Tbstr= ['Tb=' num2str(Tb)];
        f1  = [locroot mstr '/' Tbstr '/RockUplift_' num2str(ii) '.fig']
        f2  = [locroot mstr '/' Tbstr '/DynTopo_' num2str(ii) '.fig']
        openfig(f1);
        openfig(f2);
        %pause
        h1  = openfig(f1, 'invisible');
        ax1 = gca; % get handle to axes of figure
        h2  = openfig(f2, 'invisible');
        ax2 = gca;
       
        figure(2); 
            s1 = subplot(4,1,mi); box on%create and get handle to the subplot axes
            %s2 = subplot(2,1,2);
            fig1 = get(ax1,'children'); %get handle to all the children in the figure
            fig2 = get(ax2,'children');
            copyobj(fig2,s1); %copy children to new parent axes i.e. the subplot axe
            %copyobj(fig2,s2);
            tstr = ['Dynamic Topo at t = ' num2str(ii) ' my, Tb = ' num2str(Tb) ' C ' mstr ' Pa s'];
            %title(tstr);
            %pause
        mi = mi+1;
    end
end

%set(gca, 'fontname','Helvetica','fontsize',[14])
%legend('d\sigma_{xy}/dx','\sigma_{yy}', 'total')

savefig compdyntopo_revised.fig

% mi=1;
% for Tb = [800 1000 1300]
%     for muscale = [1.e+19]
%         ms   = ['1e+' num2str(log10(muscale))];
%         mstr = ['mu=' ms ];
%         Tbstr= ['Tb=' num2str(Tb)];
%         f1  = [locroot mstr '/' Tbstr '/RockUplift_' num2str(ii) '.fig']
%         f2  = [locroot mstr '/' Tbstr '/DynTopo_' num2str(ii) '.fig']
%         h1  = openfig(f1, 'invisible');
%         ax1 = gca; % get handle to axes of figure
%         h2  = openfig(f2, 'invisible');
%         ax2 = gca;
%         figure(4); 
%             s1 = subplot(3,1,mi); box on%create and get handle to the subplot axes
%             %s2 = subplot(2,1,2);
%             fig1 = get(ax1,'children'); %get handle to all the children in the figure
%             fig2 = get(ax2,'children');
%             copyobj(fig2,s1); %copy children to new parent axes i.e. the subplot axe
%             %copyobj(fig2,s2);
%             tstr = ['Dynamic Topo at t = ' num2str(ii) ' my, Tb = ' num2str(Tb) ' C ' mstr ' Pa s'];
%             %title(tstr);
%             pause(2)
%         mi = mi+1;
%     end
% end
% 
% %set(gca, 'fontname','Helvetica','fontsize',[14])
% legend('d\sigma_{xy}/dx','\sigma_{yy}', 'total')

%savefig compDynTopo_19.fig
