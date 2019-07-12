% demo: fitting experimental correlator data
% 
% DKS
% 2019-07-11

%% data
data_path = 'C:\Users\David\Dropbox\PhD\projects\bell_witness\paper\fig4\data_20181023.mat';
% data_path = 'C:\Users\hE BEC\Dropbox\PhD\projects\bell_witness\paper\fig4\data_20181023.mat';

vars_to_load = {'B','B_se','Th','tau','Tpi','S_epr','S_epr_se','Th_epr'};

S = load(data_path,vars_to_load{:});


%% Model and fit
B_mdl.fun = @Bcorr_mdl;
B_mdl.cname = {'p','alpha'};
B_mdl.par0 = [0.9,0.9*pi/2];
B_mdl.fopts=statset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e2);


% weighted fit: weight is 1/var of y
B_fit = fitnlm(S.Th,S.B,B_mdl.fun,B_mdl.par0,'Weights',(1./S.B_se).^2,...
    'CoefficientNames',B_mdl.cname,'Options',B_mdl.fopts);


% get fitted params
fit_params_mean = B_fit.Coefficients.Estimate;
fit_params_se = B_fit.Coefficients.SE;


%% Model prediction
% config
pred_nsig = 1;                  % n-sdev range: mu ± n*sigma
pred_conflvl = erf(pred_nsig/sqrt(2));  % confidence level
pred_alpha = 1-pred_conflvl;    % alpha: 100(1 – alpha)%

th=linspace(-2*pi,2*pi,1e3)';       % rotation angles for model prediction


%%% B correlator
% [B_pred,B_pred_ci]=predict(Bcorr_fit,th,'Alpha',pred_alpha);
[B_pred,B_pred_ci]=predict(B_fit,th,'Alpha',pred_alpha,'Simultaneous',true);


%%% S (EPR-steering) parameter
% get indices to where theta = 0 and pi/2
[~,I_0]=min(abs(th));
[~,I_pi2]=min(abs(th-pi/2));
Didx_pi2 = I_pi2-I_0;     % find diff in positions for a pi/2 gap

dth = th - circshift(th,-Didx_pi2);
b_ok = dth+pi/2 < 2*mean(diff(th));     % indices where wrap errors do not occur

% evaluate pi/2-shifted difference
S_pred =  abs(B_pred - circshift(B_pred,-Didx_pi2));     % S(x,x+pi/2) = |B(x) - B(x+pi/2)|
S_pred_ci = abs(B_pred_ci - circshift(B_pred_ci,-Didx_pi2));
S_pred(~b_ok) = NaN;            % set errorneous shift to NaN
S_pred_ci(~b_ok,:) = NaN;


%% Graphics config common
fig_lwidth=1;
col_gray=0.9*ones(1,3);


%% plot: B vs theta
%%% config
col_B = 'b';
col_B_theory = 'k';
col_bogo=col_gray;


%%% plot
H=figure('Name','B_vs_theta');
hold on;

% data
pB_data = ploterr(S.Th/pi,S.B,[],S.B_se,'o','hhxy',0);
set(pB_data(1),'MarkerEdgeColor',col_B,'MarkerFaceColor','w',...
    'LineWidth',fig_lwidth,'DisplayName','data');
set(pB_data(2),'Color',col_B,'LineWidth',fig_lwidth);

% fitted theory curve
pB_fit = plot(th/pi,B_pred,'-','Color',col_B_theory);
tstr = sprintf('$p = %0.2g$; $\\alpha = %0.2g$',fit_params_mean(1),fit_params_mean(2));
set(pB_fit,'DisplayName',tstr);

pB_fit_ci = plot(th/pi,B_pred_ci,'--','Color',col_B_theory);       % confidence interval

% annotation
ax=gca;
uistack(pB_fit,'bottom');
uistack(pB_fit_ci,'bottom');

box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Correlator $\mathcal{B}(\theta)$');
% title('Fit of correlator');
lgd = legend([pB_data(1),pB_fit],'Location','best','AutoUpdate','off');


%%% violation regions
% Bogoliubov nonlocality signal
B_bogo_min = 1/sqrt(2);
p_bogo=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [B_bogo_min,B_bogo_min,ax.YLim(2),ax.YLim(2)],...
    col_bogo,...
    'EdgeColor','none');
uistack(p_bogo,'bottom');      % this should REALLY be bottom - to not cover any other graphics

% text(0,0.5*(B_bogo_max+1),sprintf('Nonlocal (QM)'),'FontSize',fontsize-1,'VerticalAlignment','middle');

set(ax,'Layer','Top');     % graphics axes should be always on top

%% plot: S vs theta
%%% config
col_S = 'r';
col_S_theory = 'k';
col_epr = col_gray;
col_ent = col_gray;

%%% plot
H=figure('Name','S_vs_theta');
hold on;

% data
pS_data = ploterr(S.Th_epr/pi,S.S_epr,[],S.S_epr_se,'o','hhxy',0);
set(pS_data(1),'MarkerEdgeColor',col_S,'MarkerFaceColor','w',...
    'LineWidth',fig_lwidth,'DisplayName','data');
set(pS_data(2),'Color',col_S,'LineWidth',fig_lwidth);

% fitted theory curve
pS_fit = plot(th/pi,S_pred,'-','Color',col_S_theory);
tstr = sprintf('$p = %0.2g$; $\\alpha = %0.2g$',fit_params_mean(1),fit_params_mean(2));
set(pS_fit,'DisplayName',tstr);

pS_fit_ci = plot(th/pi,S_pred_ci,'--','Color',col_S_theory);       % confidence interval

% annotation
ax=gca;
uistack(pS_fit,'bottom');
uistack(pS_fit_ci,'bottom');

box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Steering parameter $\mathcal{S}(\theta,\theta + \pi/2)$');
lgd = legend([pS_data(1),pS_fit],'Location','best','AutoUpdate','off');

xlim(1.5*[-1,1]);
ylim([0,2]);        % algebraic range


%%% violation region
% EPR
S_epr_min=sqrt(2);
p_epr=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_epr_min,S_epr_min,2,2],...
    col_epr,'EdgeColor','none');
uistack(p_epr,'bottom');
% hatch to distinguish
hatchfill2(p_epr,'single','HatchAngle',45,'HatchDensity',30,...
    'HatchColor','k');      %,'HatchLineWidth',1
set(gca,'Layer','Top');     % graphics axes should be always on top
% text(pi/4,0.5*(S_epr_min+2),sprintf('Bell correlation witness'),'FontSize',fontsize-1,...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'BackgroundColor',col_epr);

% entanglement
S_ent_min = 1;
p_ent=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_ent_min,S_ent_min,2,2],...
    col_ent,'EdgeColor','none');
uistack(p_ent,'bottom');   
% text(pi/4,0.5*(S_ent_min+S_epr_min),sprintf('Entanglement'),'FontSize',fontsize-1,...
%     'HorizontalAlignment','center','VerticalAlignment','middle');

set(ax,'Layer','Top');     % graphics axes should be always on top