%% Violation of EPR-steering inequality with Bell+
%
% DKS
% 2018-05-29

% TODO
%   [ ] Check variability in mode occupancy
%   [ ] Update uncertainty in theta estimate
%   [x] Update with exp data for Ey  
%   [x] Update with all exp data for En; n in ZX-plane
%       1,2,3,VIOL


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
C_epr=1/4*(E_y+E_n);          % ensure vectors are row OR col
C_epr_err=1/4*sqrt(E_y_err^2+E_n_err.^2);      % simple error estimate
C0_epr=1/4*(E0_y+E0_n);
C0_epr_err=1/4*sqrt(E0_y_err^2+E0_n_err.^2);


%%% Theory: Bell+
theta_th=linspace(-pi,2*pi,1e3);
C_epr_th=0.5*sin(theta_th).^2;


%%% Inequality
C_epr_max=sqrt(2)/4;


%% Data vis
h_epr=figure('Name','epr_steering');
ax=gca;
hold on;

xlim([-0.05,1.05]);
ylim([-0.05,0.5]);

%%% Experiment
h_epr_expdata=ploterr(theta_n/pi,C_epr,theta_n_err/pi,C_epr_err,'o','hhxy',0);

% annotation
% NOTE: handle indices (1): marker; (2): X-bar; (3): Y-bar;
set(h_epr_expdata(1),'Marker','o','MarkerSize',7,...
    'Color','r','LineWidth',1.2,...
    'MarkerFaceColor','w',...
    'DisplayName','Experiment');
set(h_epr_expdata(2),'Color','r','LineWidth',1.2,...
    'DisplayName','');
set(h_epr_expdata(3),'Color','r','LineWidth',1.2,...
    'DisplayName','');


%%% Theory
p_epr_theory=plot(theta_th/pi,C_epr_th,'k--','LineWidth',1.5,...
    'DisplayName','$\vert\Psi^+\rangle$ Bell triplet');


% EPR-steering Inequality
p_epr_lim=line(ax.XLim,C_epr_max*ones(1,2),...
    'Color','k','LineStyle','-','LineWidth',1.5,...
    'DisplayName','EPR-steering');


%%% annotation
title('EPR-steering');
xlabel('$\theta/\pi$');
ylabel('$\mathcal{E}$');
box on;

lgd=legend([h_epr_expdata(1),p_epr_theory,p_epr_lim]);


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
mark_siz=7;
line_col='b';
face_col='w';
line_wid=1.2;


%%% graphics
h_bell=figure('Name','bell_test');
ax=gca;
hold on;

xlim([-0.05,1.05]);
ylim([-1,1]);

%%% Experiment
h_bell_expdata=ploterr(theta_n/pi,E_n,theta_n_err/pi,E_n_err,'o','hhxy',0);

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
p_bell_theory=plot(theta_bell_th/pi,E_bell_th,'k--','LineWidth',1.5,...
    'DisplayName','$\vert\Psi^+\rangle$ Bell triplet');


% EPR-steering Inequality
p_bell_lim=line(ax.XLim,E_bell_max*ones(1,2),...
    'Color','k','LineStyle','-','LineWidth',1.5,...
    'DisplayName','Bell inequality');


%%% annotation
title('Bell correlation');
xlabel('$\theta/\pi$');
ylabel('$\mathcal{B}$');
box on;

lgd=legend([h_bell_expdata(1),p_bell_theory,p_bell_lim]);