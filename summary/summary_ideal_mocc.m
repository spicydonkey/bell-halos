% Summary of ideal source - Bell correlations

%%% config used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% config_idealsource.m
% %% HALO
% configs.halo{1}.string='$m_F=0$';
% 
% configs.halo{2}.string='$m_F=1$';
% 
% %% Spherical zones
% configs.zone.nazim=100;
% configs.zone.nelev=50;
% 
% configs.zone.binmethod=1;
% configs.zone.binwidth=0.05;

%%%%%%%% gen_ideal_corr_source.m
% % configs
% n_shot=1e4;             % number of shots
% n_corr_pairs=60*ones(n_shot,1);        % number of correlated pairs generated from source
% det_qe=0.1*ones(1,2);       % detection efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% annotation configs
markersize=6;
linewidth=2;
% gray_col=0.5*ones(1,3);         % gray data points
namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
valarray={linewidth,'w','r'};                 % 90 deg (normal) data

%% main
nzones=50*100;      % number of "zones" == number of modes in halo

npair=[100,10,32,320,60,600];               % number of correlated pairs generated per shot
E=[-0.90,-0.99,-0.96,-0.73,-0.94,-0.59];    % evaluated correlation
Eerr=[0.01,0.01,0.01,0.01,0.01,0.01];       % SE error

nn=npair/nzones;        % mode occupancy

E_inv=abs(1./E);
E_inv_err=abs(Eerr./E).*E_inv;

figure();
% hh=ploterr(nn,1./abs(E),0,E_inv_err,'o','hhxy',0);
hh=ploterr(nn,abs(E),0,Eerr,'o','hhxy',0);
set(hh(1),namearray,valarray,'MarkerSize',markersize,'DisplayName','Experiment');	% DATAPOINT
set(hh(2),namearray,valarray,'DisplayName','');                  % Y-err
set(hh(3),namearray,valarray,'DisplayName','');                  % X-err

xlabel('Mode occupancy');
ylabel('$|E|$');

ax=gca;
ax.FontSize=13;

% linear fit
pfit=polyfit(nn,abs(E),1);
nn_fit=linspace(min(nn),max(nn),100);
Efit=polyval(pfit,nn_fit);

hold on;
hfit=plot(nn_fit,Efit,'k--','LineWidth',linewidth);
uistack(hfit,'bottom');

% inverse fit

