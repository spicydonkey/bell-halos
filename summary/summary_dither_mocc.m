% Summary of dither on counts and calculated mode occupancy


%% Results
% wbb=0.03;
n_calc=[0,0.017,0.068,0.14,0.0042,0.27];

E_max=[0.95,0.90,0.75,0.5,0.95,0.45];
E_err=[0.05,0.05,0.05,0.1,0.05,0.05];


%% Plot
% annotation configs
markersize=6;
linewidth=2;
% gray_col=0.5*ones(1,3);         % gray data points
namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
valarray={linewidth,'w','r'};                

figure();
hh=ploterr(n_calc,E_max,0,E_err,'o','hhxy',0);

% annotate
set(hh(1),namearray,valarray,'MarkerSize',markersize,'DisplayName','Experiment');	% DATAPOINT
set(hh(2),namearray,valarray,'DisplayName','');                  % Y-err
set(hh(3),namearray,valarray,'DisplayName','');                  % X-err

xlabel('Mode occupancy');
ylabel('max $|E|$');

ax=gca;
ax.FontSize=13;

% poly fit
pfit=polyfit(n_calc,E_max,2);
nn_fit=linspace(min(n_calc),max(n_calc),100);
Efit=polyval(pfit,nn_fit);

hold on;
hfit=plot(nn_fit,Efit,'k--','LineWidth',linewidth);
uistack(hfit,'bottom');

grid on;