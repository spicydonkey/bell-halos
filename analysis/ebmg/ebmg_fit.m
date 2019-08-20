%% Fitting reconstructed EBMG dBdr data to model
% DKS 2019

%% config
%TODO: need to get data with correction made

config.path_dir = 'C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\exp4_bell_t_evo\out\postproc';
config.data_name = 'out_20190301_134703.mat';

config.path_fig = fullfile(config.path_dir,config.data_name);


%% load data
S = load(config.path_fig);


%% Cartesian MODEL
mdl_ebmg_cart.fun = @(gradB,x) abs(sum(gradB.*sph2xyz(x(:,1),x(:,2),1),2));
mdl_ebmg_cart.cname = {'dBx'  'dBy'  'dBz'};
mdl_ebmg_cart.par0 = [5,0,0];
mdl_ebmg_cart.opts = statset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e2);

mdl_ebmg_cart.x = [S.gaz(:),S.gel(:)];
mdl_ebmg_cart.y = S.dBdr(:);

mdl_ebmg_cart.fit = fitnlm(mdl_ebmg_cart.x,mdl_ebmg_cart.y,mdl_ebmg_cart.fun,...
    mdl_ebmg_cart.par0,'CoefficientNames',mdl_ebmg_cart.cname,'Options',mdl_ebmg_cart.opts);


%% Cartesian MODEL: prediction
% config
pred_nsig = 1;                          % n-sdev range: mu ± n*sigma
pred_conflvl = erf(pred_nsig/sqrt(2));  % confidence level
pred_alpha = 1-pred_conflvl;            % alpha: 100(1 – alpha)%


[mdl_ebmg_cart.gaz_pred,mdl_ebmg_cart.gel_pred] = grid_azel(1e3,5e2);

mdl_ebmg_cart.xpred = [mdl_ebmg_cart.gaz_pred(:),mdl_ebmg_cart.gel_pred(:)];
[mdl_ebmg_cart.ypred,mdl_ebmg_cart.ypred_ci] = predict(mdl_ebmg_cart.fit,mdl_ebmg_cart.xpred,'Alpha',pred_alpha);
mdl_ebmg_cart.dBdr_pred = reshape(mdl_ebmg_cart.ypred,size(mdl_ebmg_cart.gaz_pred));
mdl_ebmg_cart.dBdr_pred_sdev = reshape(diff(mdl_ebmg_cart.ypred_ci,[],2)/2,size(mdl_ebmg_cart.gaz_pred));


%% VIS
%%% DATA
H_data = figure('Name','ebmg_data');

pp=plotFlatMapWrappedRad(S.gaz,S.gel,S.dBdr,'rect');


ax=gca;
set(ax,'Layer','top');
axis tight;
ylim([-90,90]);
xlabel('$\theta$');
ylabel('$\phi$');
box on;

cbar = colorbar();
cbar.Label.Interpreter = 'latex';
cbar.TickLabelInterpreter='latex';
cbar.Label.String = '$d\textrm{B}/dr$ (G/m)';


%%% MODEL
% fitted prediction
H_model = figure('Name','ebmg_model_fit');

pp=plotFlatMapWrappedRad(mdl_ebmg_cart.gaz_pred,mdl_ebmg_cart.gel_pred,mdl_ebmg_cart.dBdr_pred,'rect');


ax=gca;
set(ax,'Layer','top');
axis tight;
ylim([-90,90]);
xlabel('$\theta$');
ylabel('$\phi$');
box on;

cbar = colorbar();
cbar.Label.Interpreter = 'latex';
cbar.TickLabelInterpreter='latex';
cbar.Label.String = '$d\textrm{B}/dr$ (G/m)';


% prediction uncertainty
H_model_unc = figure('Name','ebmg_model_unc');

pp=plotFlatMapWrappedRad(mdl_ebmg_cart.gaz_pred,mdl_ebmg_cart.gel_pred,mdl_ebmg_cart.dBdr_pred_sdev,'rect');


ax=gca;
set(ax,'Layer','top');
axis tight;
ylim([-90,90]);
xlabel('$\theta$');
ylabel('$\phi$');
box on;

cbar = colorbar();
cbar.Label.Interpreter = 'latex';
cbar.TickLabelInterpreter='latex';
cbar.Label.String = 'unc. $d\textrm{B}/dr$ (G/m)';


%% equatorial prediction
mdl_ebmg_cart.az_eq = linspace(0,pi,1e3)';
[mdl_ebmg_cart.dBdr_eq_pred,mdl_ebmg_cart.dBdr_eq_pred_ci] = predict(mdl_ebmg_cart.fit,[mdl_ebmg_cart.az_eq,0*mdl_ebmg_cart.az_eq],'Alpha',pred_alpha);
mdl_ebmg_cart.dBdr_eq_pred_sdev = diff(mdl_ebmg_cart.dBdr_eq_pred_ci,[],2)/2;


%%% VIS
H_eq = figure('Name','ebmg_equator');
hold on;

% DATA
b_plot = (S.az>=0)&(S.az<=pi);

p_data = ploterr(S.az(b_plot),S.dBdr_eq(b_plot),S.configs.bins.alpha,S.dBdrerr_eq(b_plot),'hhxy',0);
set(p_data(1),'Marker','s','LineStyle','none',...
    'MarkerFaceColor','w','MarkerEdgeColor','k',...
    'DisplayName','data');
set(p_data(2),'Color','k');
set(p_data(3),'Color','k');


% % PRED
% p_model_orig = plot(mdl_ebmg_cart.az_eq,mdl_ebmg_cart.dBdr_eq_pred,'r--');

% smoothed prediction
%   moving average
l_smooth = 2*round(S.configs.bins.alpha/mean(diff(mdl_ebmg_cart.az_eq)));   % smoothing full width = bin size
yy_sm = smooth(mdl_ebmg_cart.dBdr_eq_pred,l_smooth,'moving');
p_model_sm = plot(mdl_ebmg_cart.az_eq,yy_sm,'r-','DisplayName','smoothed');


% annotation
ax = gca;

xlabel('$\theta$');
ylabel('$d\textrm{B}/dr$ (G/m)');

xlim([0,pi]);

set(ax,'Layer','top');

lgd = legend([p_data(1),p_model_sm]);