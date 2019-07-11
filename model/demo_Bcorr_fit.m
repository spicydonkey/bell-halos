% demo: fitting experimental correlator data
% 
% DKS
% 2019-07-11

%% data
data_path = 'C:\Users\David\Dropbox\PhD\projects\bell_witness\paper\fig4\data_20181023.mat';

vars_to_load = {'B','B_se','Th','tau','Tpi','S_epr','S_epr_se','Th_epr'};

S = load(data_path,vars_to_load{:});

%% model
% b = [p, alpha]
Bcorr_mdl.fun = @(b,x) -b(1) * (cos(b(2))^2 + sin(b(2))^2*cos(x)).^2 ...
    + b(1) * ( (sin(b(2))*cos(b(2))*(1-cos(x))).^2 + (sin(b(2))*sin(x)).^2);
Bcorr_mdl.cname = {'p','alpha'};
Bcorr_mdl.par0 = [0.9,0.8562*pi/2];
Bcorr_mdl.fopts=statset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e2);

% do fit!
Bcorr_fit = fitnlm(S.Th,S.B,Bcorr_mdl.fun,Bcorr_mdl.par0,...
    'CoefficientNames',Bcorr_mdl.cname,'Options',Bcorr_mdl.fopts);

fit_params_mean = Bcorr_fit.Coefficients.Estimate;
fit_params_se = Bcorr_fit.Coefficients.SE;

% eval fitted oscillation
% th_fit=linspace(0,2*pi,1e3);
th_fit=linspace(0,2*pi,1e3);
B_fit=feval(Bcorr_fit,th_fit);


%% plot: B vs theta
H=figure('Name','B_vs_theta');
hold on;

% data
col_B = 'b';
p_data = ploterr(S.Th/pi,S.B,[],S.B_se,'o','hhxy',0);
set(p_data(1),'MarkerEdgeColor',col_B,'MarkerFaceColor','w','DisplayName','data');
set(p_data(2),'Color',col_B);

% fitted theory curve
col_theory = 'k';
p_fit = plot(th_fit/pi,B_fit,'-','Color',col_theory);
tstr = sprintf('$p = %0.2g$; $\\alpha = %0.2g$',fit_params_mean(1),fit_params_mean(2));
set(p_fit,'DisplayName',tstr);

% annotation
box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('Correlator $\mathcal{B}(\theta)$');
% title('Fit of correlator');
lgd = legend([p_data(1),p_fit]);


%% plot: S vs theta
% get indices to where theta = 0 and pi/2
[~,I_0]=min(abs(th_fit));
[~,I_pi2]=min(abs(th_fit-pi/2));
Didx_pi2 = I_pi2-I_0;     % find diff in positions for a pi/2 gap

S_fit =  B_fit - circshift(B_fit,-Didx_pi2);     % S(x,x+pi/2) = B(x) - B(x+pi/2)

% plot
H=figure('Name','S_vs_theta');
hold on;

% data
col_S = 'r';
p_data = ploterr(S.Th_epr/pi,S.S_epr,[],S.S_epr_se,'o','hhxy',0);
set(p_data(1),'MarkerEdgeColor',col_S,'MarkerFaceColor','w','DisplayName','data');
set(p_data(2),'Color',col_S);

% fitted theory curve
p_fit = plot(th_fit/pi,abs(S_fit),'k-');
tstr = sprintf('$p = %0.2g$; $\\alpha = %0.2g$',fit_params_mean(1),fit_params_mean(2));
set(p_fit,'DisplayName',tstr);

% annotation
box on;
xlabel('Rotation angle $\theta/\pi$');
ylabel('$\mathcal{S}(\theta,\theta + \pi/2)$');
lgd = legend([p_data(1),p_fit]);

