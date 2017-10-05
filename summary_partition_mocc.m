% Summary of variable mode occupancy on different sphere partition sizes

%% Results
n_calc=repmat([0.017,0.034,0.068,0.0085,0.0061,0.0044],[2,1]);

E_max=[
    0.89,0.79,0.64,0.96,0.97,0.98;...
    0.93,0.89,0.73,1,1,1;...
    ];

E_err=1e-2*[
    2,2,1,1,1,1;...
    3,3,3,0,0,0;...
    ];

config_str={'$\delta\theta = 2 w_{BB}$','$\delta\theta = w_{BB}$'};

%% Plot
% annotation configs
markersize=6;
linewidth=2;
% gray_col=0.5*ones(1,3);         % gray data points
namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
valarray={{linewidth,'w','r'},{linewidth,'w','k'}};                

figure();
for ii=1:2
    hold on;
    
    hh{ii}=ploterr(n_calc(ii,:),E_max(ii,:),0,E_err(ii,:),'o','hhxy',0);
    
    % annotate
    set(hh{ii}(1),namearray,valarray{ii},'MarkerSize',markersize,'DisplayName',config_str{ii});	% DATAPOINT
    set(hh{ii}(2),namearray,valarray{ii},'DisplayName','');                  % Y-err
    set(hh{ii}(3),namearray,valarray{ii},'DisplayName','');                  % X-err
    
    hold off;
end

legend([hh{1}(1),hh{2}(1)]);

xlabel('Mode occupancy');
ylabel('max $|E|$');

ax=gca;
ax.FontSize=13;

grid on;
box on;



% % poly fit
% pfit=polyfit(n_calc,E_max,2);
% nn_fit=linspace(min(n_calc),max(n_calc),100);
% Efit=polyval(pfit,nn_fit);
% 
% hold on;
% hfit=plot(nn_fit,Efit,'k--','LineWidth',linewidth);
% uistack(hfit,'bottom');
% 
% grid on;