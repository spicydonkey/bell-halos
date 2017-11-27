% Summary of variable bin-size of spherical partition on correlation
% characterised bin-size by #MODES/#BINS


%% Results
n_modes_per_bin=1./[1.1,0.56,0.29,0.16,2.3,4.5,36,4.4,18];

E_max=[0.83,0.80,0.69,0.52,0.90,0.90,0.98,0.92,0.97];

E_err=1e-2*[3,2,2,2,4,5,10,5,9];

% config_str={'$\delta\theta = 2 w_{BB}$','$\delta\theta = w_{BB}$'};

%% Plot
% annotation configs
markersize=6;
linewidth=2;
gray_col=0.5*ones(1,3);         % gray data points
namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
valarray={linewidth,'w','r'};

figure();
hh=ploterr(n_modes_per_bin,E_max,0,E_err,'o','hhxy',0);

% annotate
set(hh(1),namearray,valarray,'MarkerSize',markersize,'DisplayName','Simulation');	% DATAPOINT
set(hh(2),namearray,valarray,'DisplayName','');                  % Y-err
set(hh(3),namearray,valarray,'DisplayName','');                  % X-err

xlabel('No. modes per bin');
ylabel('max $|E|$');

ax=gca;
ax.FontSize=13;

grid on;
box on;

ylim auto;
ylim_auto=ylim;
ylim([ylim_auto(1),1]);

% poly fit
pfit=polyfit(n_modes_per_bin,E_max,1);
nn_fit=linspace(min(n_modes_per_bin),max(n_modes_per_bin),100);
Efit=polyval(pfit,nn_fit);

hold on;
hfit=plot(nn_fit,Efit,'k--','LineWidth',linewidth);
uistack(hfit,'bottom');

grid on;