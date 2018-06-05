%% Violation of EPR-steering inequality with Bell+
%
% DKS
% 2018-05-29

% TODO
%   [ ] Check variability in mode occupancy
%   [ ] Update uncertainty in theta estimate
%   [x] Update with exp data for Ey  
%   [x] Update with all exp data for En; n in ZX-plane
%       1,2,3,VIOL,ideal_global_rotation_v4 (pi/4)

%% VIS config
font_siz_reg=12;
font_siz_sml=10;
font_siz_lrg=14;
mark_siz=7;
line_wid=1.1;

% [cc,clight,cdark]=palette(n_mf);
% mark_typ={'o','^','d'};



%% Experimental data (as code)
% linear approximation of Raman pulse duration and rotation angle
Tpi=10e-6;      % approx pi pulse duration

theta_frac_unc=3e-2;        % fractional uncertainty in rotation angle


%%% 1. <Jn(A).Jn(B)>/<N(A).N(B)>; n in ZX-plane ( n = z*cos(th) + x*sin(th) )
%%%% EXP5_BELL_YROT
% Run 1
idx=1;
T_raman_n_cell{idx}=1e-5*[0    0.0625    0.1250    0.1875    0.3125    0.3750    0.4375];
E_n_cell{idx}=[-0.8727   -0.8474   -0.5811   -0.1800    0.4924    0.6892    0.8560];
E_n_err_cell{idx}=[0.0898    0.0477    0.1095    0.1539    0.1650    0.1002    0.0776];
E0_n_cell{idx}=[-0.9792   -0.9189   -0.6240   -0.1938    0.5279    0.7312    0.9039];
E0_n_err_cell{idx}=[0.1064    0.0562    0.1164    0.1680    0.1795    0.1098    0.0799];

% Run 2
idx=2;
T_raman_n_cell{idx}=1e-5*[0.6250    0.7500    0.8750    1.0000];
E_n_cell{idx}=[0.6281    0.2059   -0.2843   -0.4375];
E_n_err_cell{idx}=[0.2837    0.2607    0.3886    0.2518];
E0_n_cell{idx}=[0.6659    0.2181   -0.3060   -0.4758];
E0_n_err_cell{idx}=[0.3248    0.2828    0.4360    0.2813];

% Run 3
idx=3;
T_raman_n_cell{idx}=1e-5*[0.5625    0.6875    0.8125    0.9375];
E_n_cell{idx}=[0.7830    0.2093   -0.2292   -0.5635];
E_n_err_cell{idx}=[0.1060    0.1439    0.1146    0.1505];
E0_n_cell{idx}=[0.8880    0.2380   -0.2607   -0.6488];
E0_n_err_cell{idx}=[0.1325    0.1642    0.1347    0.1788];


%%%% EXP1_BELL_MAX_VIOL
idx=4;
T_raman_n_cell{idx}=5e-6;
E_n_cell{idx}=0.85;
E_n_err_cell{idx}=0.05;
E0_n_cell{idx}=0.9395;
E0_n_err_cell{idx}=0.0725;


%%%% IDEAL_GLOBAL_ROTATION_V4 (slightly higher mode occupancy)
% picking out theta~=pi/4 (par_T=2.5us)
idx=5;
T_raman_n_cell{idx}=1e-5*[0.2500];
E_n_cell{idx}=[0.0933];
E_n_err_cell{idx}=[0.0815];
E0_n_cell{idx}=[0.1112];
E0_n_err_cell{idx}=[0.0978];


% collate cell-data as array
T_raman_n=horzcat(T_raman_n_cell{:});
E_n=horzcat(E_n_cell{:});
E_n_err=horzcat(E_n_err_cell{:});
E0_n=horzcat(E0_n_cell{:});
E0_n_err=horzcat(E0_n_err_cell{:});


%%% 2. <Jy(A).Jy(B)>/<N(A).N(B)>
E_y=0.835;
E_y_err=0.06;
E0_y=0.926;
E0_y_err=0.08;



%% Set up n-axis - theta
theta_n=T_raman_n/Tpi*pi;       % simple linear model
theta_n_err=theta_frac_unc*theta_n;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EPR-steering
%%% Experiment
% coefficient
S_epr=1/4*(E_y+E_n);          % ensure vectors are row OR col
S_epr_err=1/4*sqrt(E_y_err^2+E_n_err.^2);      % simple error estimate
S0_epr=1/4*(E0_y+E0_n);
S0_epr_err=1/4*sqrt(E0_y_err^2+E0_n_err.^2);


%%% Theory: Bell+
theta_th=linspace(-pi,2*pi,1e3);
S_epr_th=0.5*sin(theta_th).^2;


%%% Inequality
S_epr_max=sqrt(2)/4;


%% Data vis
%%% config
mark_typ='d';
line_col=[1 0 0];
face_col=colorspace('HSL->RGB',[1,1,1.7].*colorspace('RGB->HSL',line_col));
ptch_col=0.85;           % some gray to patch violation zone
ptch_alp=1;             % patch transparency
x_lim=pi*[-0.05,1.05];     % set to [] for auto
y_lim=[-0.05,0.5];

%%% graphics
h_epr=figure('Name','epr_steering');
ax=gca;
hold on;

if ~isempty(x_lim)
    xlim(x_lim)
end
if ~isempty(y_lim)
    ylim(y_lim)
end

%%% Experiment
h_epr_expdata=ploterr(theta_n,S_epr,theta_n_err,S_epr_err,mark_typ,'hhxy',0);

% annotation
% NOTE: handle indices (1): marker; (2): X-bar; (3): Y-bar;
set(h_epr_expdata(1),'Marker',mark_typ,'MarkerSize',mark_siz,...
    'Color',line_col,'LineWidth',line_wid,...
    'MarkerFaceColor',face_col,...
    'DisplayName','Experiment');
set(h_epr_expdata(2),'Color',line_col,'LineWidth',line_wid,...
    'DisplayName','');
set(h_epr_expdata(3),'Color',line_col,'LineWidth',line_wid,...
    'DisplayName','');


%%% Theory
p_epr_theory=plot(theta_th,S_epr_th,'k--','LineWidth',1.5,...
    'DisplayName','$\vert\Psi^+\rangle$');
uistack(p_epr_theory,'bottom');


%%% EPR-steering Inequality
p_epr_lim=line(ax.XLim,S_epr_max*ones(1,2),...
    'Color','k','LineStyle','-','LineWidth',1.5,...
    'DisplayName','EPR-steering limit');
uistack(p_epr_lim,'bottom');

% patch
p_epr_viol=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_epr_max,S_epr_max,ax.YLim(2),ax.YLim(2)],...
    ptch_col*ones(1,3),'FaceAlpha',ptch_alp,...
    'EdgeColor','none');
uistack(p_epr_viol,'bottom');


%%% annotation
set(gca,'Layer','Top');     % graphics axes should be always on top
% title('EPR-steering');
xlabel('$\theta$');
ylabel('$\mathcal{S}$');
box on;

% lgd=legend([h_epr_expdata(1),p_epr_theory,p_epr_lim]);

% fontsize
set(gca,'FontSize',font_siz_reg);
% set(lgd,'FontSize',font_siz_reg);

% X-ticks
xticks(0:pi/4:pi);
xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bell correlations

%%% Experiment
% Bell correlation
%   same as E_n* values


%%% Theory: Bell+
theta_bell_th=linspace(-pi,2*pi,1e3);
E_bell_th=-cos(2*theta_bell_th);


%%% Bell inequality
E_bell_max=1/sqrt(2);


%% data vis
%%% config
mark_typ='o';
line_col=[0 0 1];
face_col=colorspace('HSL->RGB',[1,1,1.7].*colorspace('RGB->HSL',line_col));
ptch_col=0.85;           % some gray to patch violation zone
ptch_alp=1;             % patch transparency
x_lim=pi*[-0.05,1.05];     % set to [] for auto
y_lim=[-1,1];


%%% graphics
h_bell=figure('Name','bell_test');
ax=gca;
hold on;

if ~isempty(x_lim)
    xlim(x_lim)
end
if ~isempty(y_lim)
    ylim(y_lim)
end

%%% Experiment
h_bell_expdata=ploterr(theta_n,E_n,theta_n_err,E_n_err,mark_typ,'hhxy',0);


% annotation
% NOTE: handle indices (1): marker; (2): X-bar; (3): Y-bar;
set(h_bell_expdata(1),'Marker',mark_typ,'MarkerSize',mark_siz,...
    'Color',line_col,'LineWidth',line_wid,...
    'MarkerFaceColor',face_col,...
    'DisplayName','Experiment');
set(h_bell_expdata(2),'Color',line_col,'LineWidth',line_wid,...
    'DisplayName','');
set(h_bell_expdata(3),'Color',line_col,'LineWidth',line_wid,...
    'DisplayName','');

%%% Theory
p_bell_theory=plot(theta_bell_th,E_bell_th,'k--','LineWidth',1.5,...
    'DisplayName','$\vert\Psi^+\rangle$');
uistack(p_bell_theory,'bottom');

%%% Bell Inequality
p_bell_lim=line(ax.XLim,E_bell_max*ones(1,2),...
    'Color','k','LineStyle','-','LineWidth',1.5,...
    'DisplayName','Bell inequality limit');
uistack(p_bell_lim,'bottom');

% patch
p_bell_viol=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [E_bell_max,E_bell_max,ax.YLim(2),ax.YLim(2)],...
    ptch_col*ones(1,3),'FaceAlpha',ptch_alp,...
    'EdgeColor','none');
uistack(p_bell_viol,'bottom');


%%% annotation
set(gca,'Layer','Top');     % graphics axes should be always on top
% title('Bell correlation');
xlabel('$\theta$');
ylabel('$\mathcal{E}$');
box on;

% lgd=legend([h_bell_expdata(1),p_bell_theory,p_bell_lim]);

% fontsize
set(gca,'FontSize',font_siz_reg);
% set(lgd,'FontSize',font_siz_reg);

% X-ticks
xticks(0:pi/4:pi);
xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});