% Nonlocal spin correlations
%   fitted with realistic state + rotation model
% 
% DKS 2019 - 2020


%% CONFIG: data
data_path = 'C:\Users\David\Dropbox\PhD\projects\bell_witness\paper\fig4\data_20181023.mat';
% data_path = 'C:\Users\hE BEC\Dropbox\PhD\projects\bell_witness\paper\fig4\data_20181023.mat';

vars_to_load = {'B','B_se','Th','tau','Tpi','S_epr','S_epr_se','Th_epr'};

S = load(data_path,vars_to_load{:});


%% CONFIG: VIS
% general
f_units='centimeters';
f_pos=[0 0 12.5 8];
f_ren='painters';

% [c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
% str_ss={'$(1,1)$','$(0,0)$','$(1,0)$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
ax_lwid=1.2;
fontsize=12;

% custom
% col_Bcorr=[0.5 0.3 0.7];               % main color for correlation
col_Bcorr=[0 0 1];               % main color for correlation
[col_Bcorr_l,col_Bcorr_d]=colshades(col_Bcorr);      

col_S = [1,0,0];
[col_S_l,col_S_d]=colshades(col_S);     

% col_gray=0.9*ones(1,3);
col_psip = 'k';
% col_bogo=0.85*ones(1,3);
col_qent=0.92*ones(1,3);
col_epr=0.8*ones(1,3);
col_alpha = 0.33;     

% c_patch=0.85;         	% some gray to patch violation zone
% a_patch=1;       		% patch transparency

lim_th=[-pi/20,pi+pi/20];   % theta-axis limits


%%% nonlocality
S_ent_min = 1;
S_epr_min=sqrt(2);


%% Model and fit
B_mdl.fun = @mdl_rhop_offres;       % @Bcorr_mdl
B_mdl.cname = {'p','delta'};
B_mdl.par0 = [0.9,0.1];               %[0.9,0.9*pi/2];
B_mdl.fopts=statset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e2);


% weighted fit: weight is 1/var of y
B_fitmdl = fitnlm(S.Th,S.B,B_mdl.fun,B_mdl.par0,'Weights',(1./S.B_se).^2,...
    'CoefficientNames',B_mdl.cname,'Options',B_mdl.fopts);


% get fitted params
fit_params_mean = B_fitmdl.Coefficients.Estimate;
fit_params_se = B_fitmdl.Coefficients.SE;


%% Model prediction
% config
pred_nsig = 1;                  % n-sdev range: mu ± n*sigma
pred_conflvl = erf(pred_nsig/sqrt(2));  % confidence level
pred_alpha = 1-pred_conflvl;    % alpha: 100(1 – alpha)%

th_fit=linspace(-2*pi,2*pi,1e4)';       % rotation angles for model prediction


%%% B correlator
% [B_pred,B_pred_ci]=predict(Bcorr_fit,th,'Alpha',pred_alpha);
[B_pred,B_pred_ci]=predict(B_fitmdl,th_fit,'Alpha',pred_alpha,'Simultaneous',true);


%%% S (EPR-steering) parameter
% get indices to where theta = 0 and pi/2
[~,I_0]=min(abs(th_fit));
[~,I_pi2]=min(abs(th_fit-pi/2));
Didx_pi2 = I_pi2-I_0;     % find diff in positions for a pi/2 gap

dth = th_fit - circshift(th_fit,-Didx_pi2);
b_ok = dth+pi/2 < 2*mean(diff(th_fit));     % indices where wrap errors do not occur


% evaluate pi/2-shifted difference
S_pred =  abs(B_pred - circshift(B_pred,-Didx_pi2));     % S(x,x+pi/2) = |B(x) - B(x+pi/2)|
S_pred_ci = abs(B_pred_ci - circshift(B_pred_ci,-Didx_pi2));
S_pred(~b_ok) = NaN;            % set errorneous shift to NaN
S_pred_ci(~b_ok,:) = NaN;


%% Psip
par_psip = [1,pi/2];        % ideal Psi+ state parameter

B_psip = Bcorr_mdl(par_psip,th_fit);

Scorr_mdl = @(beta,th) abs(Bcorr_mdl(beta,th) - Bcorr_mdl(beta,th+pi/2));

S_psip = Scorr_mdl(par_psip,th_fit);



%% plot: B vs theta
h_Bcorr = figure('Name','Bcorr','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

% data
pB_data = ploterr(S.Th/pi,S.B,[],S.B_se,'o','hhxy',0);
set(pB_data(1),'MarkerEdgeColor',col_Bcorr,'MarkerFaceColor',col_Bcorr_l,...
    'LineWidth',line_wid,'DisplayName','data');
set(pB_data(2),'Color',col_Bcorr,'LineWidth',line_wid);

% fit
pB_fit = plot(th_fit/pi,B_pred,'-','Color',col_Bcorr,'LineWidth',line_wid);
tstr = sprintf('$\\hat{\\rho}_+(p = %0.2g)$,\n$\\delta = %0.2g$',fit_params_mean(1),fit_params_mean(2));
set(pB_fit,'DisplayName',tstr);
pB_fit.Color(4) = col_alpha;
% pB_fit_ci = plot(th_fit/pi,B_pred_ci,'--','Color',col_Bcorr);       % confidence interval

% Psip
pB_psip = plot(th_fit/pi,B_psip,'Color',col_psip,'LineStyle','--',...
    'DisplayName','$\vert\Psi^+\rangle$');


%%% annotation
ax=gca;
set(ax,'Layer','Top');     % graphics axes should be always on top
uistack(pB_fit,'bottom');
% uistack(pB_fit_ci,'bottom');
uistack(pB_psip,'bottom');

box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Correlator $\mathcal{B}(\theta)$');

xlim(1.1*[-1,1]);

% title('Fit of correlator');
lgd = legend([pB_fit,pB_psip],'Location','best','AutoUpdate','off');
lgd.Box = 'off';
lgd.Position=[0.2148    0.1528    0.2547    0.1542];

% %%% violation regions
% % Bogoliubov nonlocality signal
% B_bogo_min = 1/sqrt(2);
% p_bogo=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
%     [B_bogo_min,B_bogo_min,ax.YLim(2),ax.YLim(2)],...
%     col_bogo,...
%     'EdgeColor','none');
% uistack(p_bogo,'bottom');      % this should REALLY be bottom - to not cover any other graphics

% text(0,0.5*(B_bogo_max+1),sprintf('Nonlocal (QM)'),'FontSize',fontsize-1,'VerticalAlignment','middle');


%% plot: S vs theta
% %%% config
% col_S = 'r';
% col_S_theory = 'k';
% col_epr = col_gray;
% col_ent = col_gray;

%%% plot
h_S = figure('Name','Sepr','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

% data
pS_data = ploterr(S.Th_epr/pi,S.S_epr,[],S.S_epr_se,'o','hhxy',0);
set(pS_data(1),'MarkerEdgeColor',col_S,'MarkerFaceColor',col_S_l,...
    'LineWidth',line_wid,'DisplayName','data');
set(pS_data(2),'Color',col_S,'LineWidth',line_wid);

% fit
pS_fit = plot(th_fit/pi,S_pred,'-','Color',col_S,'LineWidth',line_wid);
tstr = sprintf('$\\hat{\\rho}_+(p = %0.2g)$,\n$\\delta = %0.2g$',fit_params_mean(1),fit_params_mean(2));
set(pS_fit,'DisplayName',tstr);
pS_fit.Color(4) = col_alpha;
% pS_fit_ci = plot(th_fit/pi,S_pred_ci,'--','Color',col_S);

% Psip
pS_psip = plot(th_fit/pi,S_psip,'Color',col_psip,'LineStyle','--',...
    'DisplayName','$\vert\Psi^+\rangle$');


%%% annotation
ax=gca;
uistack(pS_fit,'bottom');
% uistack(pS_fit_ci,'bottom');
uistack(pS_psip,'bottom');

box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Steering parameter $\mathcal{S}(\theta,\theta + \pi/2)$');

xlim([-1,1]);
ylim([0,2]);        % algebraic range


%%% violation region
% EPR
p_epr1=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_epr_min,S_epr_min,2,2],col_epr,'EdgeColor','none',...
    'DisplayName','EPR-steering');
uistack(p_epr1,'bottom');
% hatch to distinguish
% p_epr2=hatchfill2(p_epr1,'single','HatchAngle',45,'HatchDensity',30,...
%     'HatchColor','k');      %,'HatchLineWidth',1
set(gca,'Layer','Top');     % graphics axes should be always on top
% text(pi/4,0.5*(S_epr_min+2),sprintf('Bell correlation witness'),'FontSize',fontsize-1,...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'BackgroundColor',col_epr);

% entanglement
p_ent=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_ent_min,S_ent_min,2,2],col_qent,'EdgeColor','none',...
    'DisplayName','Entangled');
uistack(p_ent,'bottom');   
% text(pi/4,0.5*(S_ent_min+S_epr_min),sprintf('Entanglement'),'FontSize',fontsize-1,...
%     'HorizontalAlignment','center','VerticalAlignment','middle');

set(ax,'Layer','Top');     % graphics axes should be always on top

% lgd = legend([pS_fit,pS_psip,p_ent,p_epr1],'Location','best','AutoUpdate','off');
lgd = legend([p_ent,p_epr1],'Location','best','AutoUpdate','off');
lgd.Box = 'off';
lgd.Location="southwest";