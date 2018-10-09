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


%% evaluate normalised g2
g2_amp_tot=0.5*(sum(g2_amp,2)+g2_amp(:,3));
g2_amp_norm=g2_amp./g2_amp_tot;
g2_se_norm=g2_se./g2_amp_tot;       % scale err

% smoothed fit to data
polyfit_gnorm=cell(1,3);
ninterp=1e4;
tau_sm_fit=linspace(min(tau),max(tau),ninterp);
gnorm_sm_fit=NaN(3,ninterp);
for ii=1:3
    polyfit_gnorm{ii}=polyfit(tau,g2_amp_norm(:,ii),6);
    gnorm_sm_fit(ii,:)=polyval(polyfit_gnorm{ii},tau_sm_fit);
end

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

% fitted line
p_gnorm_sm=NaN(1,3);
for ii=1:3
    p_gnorm_sm(ii)=plot(1e6*tau_sm_fit,gnorm_sm_fit(ii,:),'LineWidth',2.2,'Color',0.5*(c0(ii,:)+clight(ii,:)));
    uistack(p_gnorm_sm(ii),'bottom');
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
