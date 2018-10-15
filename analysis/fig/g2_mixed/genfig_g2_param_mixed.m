%% g2 parameters under unitary rotation (mixing spins)
% DKS
% 20181009

% config
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot\bell_signal';
data_fname='bell_20181009_bs_fix.mat';      % correct g2 stderr from bootstrap
path_data=fullfile(path_dir,data_fname);

% load data
load(path_data);    % load into workspace

%% collate data for vis
tau=cat(1,S.par_T);     % rotation pulse duration

%%% g2 PARAMS
%   AMPLITUDE
%   eval'd value at dk=0 ==> used in analysis
%       g2_anti=cat(1,S.g2anti_par);
%       g2_corr=cat(1,S.g2corr_par);
%
%   3D GAUSSIAN fit:
%     y ~ 1 + amp*exp( - 0.5*(((x1 - mu1)/sig1)^2 + ((x2 - mu2)/sig2)^2 + ((x3 - mu3)/sig3)^2))
%   CoefficientNames: ('amp','mu1','mu2','mu3','sig1','sig2','sig3')
%

n_dataset=numel(S);
g2_amp_cell=cell(n_dataset,1);      % amplitude
g2_se_cell=cell(n_dataset,1);       % stderr of amp
g2_fit_param_cell=cell(n_dataset,1);     % fit params
for ii=1:n_dataset
    tn_par=numel(S(ii).par_T);
    g2_amp_cell{ii}=NaN(tn_par,3);
    g2_se_cell{ii}=NaN(tn_par,3);
    g2_fit_param_cell{ii}=cell(tn_par,3);
    tidx0=S(ii).idx_dk0;
    for jj=1:tn_par
        % get amp/SE from different S-S corr types
        g2_amp_cell{ii}(jj,:)=cellfun(@(x) x(tidx0,tidx0,tidx0),S(ii).g2{jj});
        g2_se_cell{ii}(jj,:)=cellfun(@(x) x(tidx0,tidx0,tidx0),S(ii).g2_sdev{jj});
        
        % get model params (param-"vectors" from corr-type idxd by cell)
        g2_fit_param_cell{ii}(jj,:)=cellfun(@(f) f.Coefficients.Estimate',S(ii).g2mdl{jj},'UniformOutput',false);
    end
end
% form into array
g2_amp=cat(1,g2_amp_cell{:});       
g2_se=cat(1,g2_se_cell{:});       

g2_fit_param=cell(1,3);     % a little more tedious for model params
temp_g2_fit_param=cell(n_dataset,3);     
for ii=1:n_dataset
    for jj=1:3
        temp_g2_fit_param{ii,jj}=cat(1,g2_fit_param_cell{ii}{:,jj});
    end
end
for ii=1:3
    g2_fit_param{ii}=cat(1,temp_g2_fit_param{:,ii});
end
sig_bb=cellfun(@(p) abs(p(:,5:7)),g2_fit_param,'UniformOutput',false);	% fitted g2BB sigma


%% evaluate normalised g2
g2_amp_tot=0.5*(sum(g2_amp,2)+g2_amp(:,3));
g2_amp_norm=g2_amp./g2_amp_tot;
g2_se_norm=g2_se./g2_amp_tot;       % scale err


%%% fit
ninterp=1e4;
tau_fit=linspace(min(tau),max(tau),ninterp);

% %%% polynomial fit
% polyfit_gnorm=cell(1,3);
% gnorm_fit_poly=NaN(3,ninterp);
% for ii=1:3
%     polyfit_gnorm{ii}=polyfit(tau,g2_amp_norm(:,ii),6);
%     gnorm_fit_poly(ii,:)=polyval(polyfit_gnorm{ii},tau_fit);
% end

%%% damped sine fit
fo = statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);
% % Model 1: amp,offset,phase,decayrate
% modelfun='y~c+a*exp(-x1/xc)*sin(2*pi*0.1*x1+phi)';
% modelcoef={'a','c','phi','xc'};
% par0={[0.5 0.5 -pi/2 20], [0.5 0.5 -pi/2 20], [0.5 0.5 pi/2 20]};
% 
% sinefit_gnorm=cell(1,3);
% gnorm_fit_sine=NaN(3,ninterp);
% for ii=1:3
%     sinefit_gnorm{ii}=fitnlm(1e6*tau,g2_amp_norm(:,ii),modelfun,par0{ii},'CoefficientNames',modelcoef,'Options',fo);
%     gnorm_fit_sine(ii,:)=feval(sinefit_gnorm{ii},1e6*tau_fit);
% end

% Model 2: amp,offset,decayrate (i.e. freq, phase are fixed for simplicity)
modelfun=cell(1,3);
modelfun{1}='y~c+a*exp(-x1/xc)*sin(2*pi*0.1*x1-pi/2)';
modelfun{2}=modelfun{1};
modelfun{3}='y~c+a*exp(-x1/xc)*sin(2*pi*0.1*x1+pi/2)';
modelcoef={'a','c','xc'};
par0={[0.5 0.5 20], [0.5 0.5 20], [0.5 0.5 20]};

sinefit_gnorm=cell(1,3);
gnorm_fit_sine=NaN(3,ninterp);
for ii=1:3
    sinefit_gnorm{ii}=fitnlm(1e6*tau,g2_amp_norm(:,ii),modelfun{ii},par0{ii},'CoefficientNames',modelcoef,'Options',fo);
    gnorm_fit_sine(ii,:)=feval(sinefit_gnorm{ii},1e6*tau_fit);
end

%% evaluate correlation volume
% TODO:
%   need to check fits

% mean BB corr length
sig_bb_mean=cellfun(@(s) geomean(s,2),sig_bb,'UniformOutput',false);
sig_bb_mean=cat(2,sig_bb_mean{:});

% approx BB corr volume
beta_mode=(2*pi)^(3/2);     % geo-factor by Gauss approx of BEC
V_bb=beta_mode*sig_bb_mean.^3;


%% vis
% grahics configs
[c0,clight,cdark]=palette(3);
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
% str_ss={'$(1,1)$','$(0,0)$','$(1,0)$'};
mark_ss={'o','s','^'};
mark_siz=8;
line_wid=2;

%%% g2 amplitude
hfig=figure('Name','g2_vs_theta',...
    'Units','normalized','Position',[0.2,0.2,0.5,0.25],'Renderer','painters');
hold on;
pleg=NaN(1,3);      % data for legend
for ii=1:3
    tp=myploterr(1e6*tau,g2_amp(:,ii),[],g2_se(:,ii),mark_ss{ii},c0(ii,:));
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
    set(tp(2),'LineWidth',line_wid,'DisplayName','');
end
ax=gca;
lgd=legend(pleg);
lgd.FontSize=12;
box on;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('$g^{(2)}_\textrm{BB}$ amplitude');
ax.FontSize=12;
ax.LineWidth=1.2;



%%% NORMalised g2 amp G
hfig=figure('Name','G_vs_theta',...
    'Units','normalized','Position',[0.2,0.2,0.5,0.25],'Renderer','painters');
hold on;
pleg=NaN(1,3);      % data for legend
for ii=1:3
    tp=myploterr(1e6*tau,g2_amp_norm(:,ii),[],g2_se_norm(:,ii),mark_ss{ii},c0(ii,:));
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
    set(tp(2),'LineWidth',line_wid,'DisplayName','');
end

%%% fitted line
%%% poly-fit: WARNING: simple to implement but often displays misleading
%%% features!
% p_gnorm_poly=NaN(1,3);
% for ii=1:3
%     p_gnorm_poly(ii)=plot(1e6*tau_fit,gnorm_fit_poly(ii,:),'LineWidth',2.2,'Color',0.5*(c0(ii,:)+clight(ii,:)));
%     uistack(p_gnorm_poly(ii),'bottom');
% end

% damped sine fit: closer to underlying physics
p_gnorm_sine=NaN(1,3);
for ii=1:3
    p_gnorm_sine(ii)=plot(1e6*tau_fit,gnorm_fit_sine(ii,:),'LineWidth',2.2,'Color',0.5*(c0(ii,:)+clight(ii,:)));
    uistack(p_gnorm_sine(ii),'bottom');
end


% annotation
ax=gca;
lgd=legend(pleg);
lgd.FontSize=12;
box on;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Normalised correlation $\mathcal{G}$');
ax.FontSize=12;
ax.LineWidth=1.2;
