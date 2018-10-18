%% Uncertainty estimation of g2 correlation (g2BB) by bootstrapping
% DKS
% 20181018

%% CONFIG
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';
D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_ismat=cellfun(@(s) strcmp(fileext(s),'.mat'),filename);
dataname=filename(b_ismat);
n_data=numel(dataname);

%% main
% preallocate var size (as cell to be arrayed)
tau=cell(n_data,1);
g2BB=cell(n_data,1);
g2BB_samp=cell(n_data,1);
g2BB_se=cell(n_data,1);
g2BB_mean=cell(n_data,1);
param_bs=cell(n_data,1);    % [n_frac_samp, n_subset]

for ii=1:n_data     % SKIP exp1_collated
    %% load data
    tdataname=dataname{ii};
    tS=load(fullfile(path_dir,tdataname));
    tidx0=tS.idx_dk0;
    
    %% routine
    tau{ii}=tS.par_T;
    %%% get g2(0)
    tg2BB=cellfun(@(c) cellfun(@(k) k(tidx0,tidx0,tidx0),c),tS.g2,'UniformOutput',false);
    tg2BB=cat(1,tg2BB{:});
    g2BB{ii}=tg2BB;
    %%% get g2(0) from bootstrapping
    tg2BB_samp=cellfun(@(c) cellfun(@(k) squeeze(k(tidx0,tidx0,tidx0,:)),c,'UniformOutput',false),...
        tS.g2_bootstrap,'UniformOutput',false);
    tg2BB_samp=cellfun(@(c) cat(2,c{:}),tg2BB_samp,'UniformOutput',false);
    g2BB_samp{ii}=tg2BB_samp;
    %%% get error estimate (bootstrapping)
    tg2BB_se=cellfun(@(x) sqrt(tS.n_frac_samp)*std(x,0,1),tg2BB_samp,'UniformOutput',false);
    tg2BB_se=cat(1,tg2BB_se{:});
    g2BB_se{ii}=tg2BB_se;
    % mean
    tg2BB_mean=cellfun(@(x) mean(x,1),tg2BB_samp,'UniformOutput',false);
    tg2BB_mean=cat(1,tg2BB_mean{:});
    g2BB_mean{ii}=tg2BB_mean;
    % bootstrap parameter
    param_bs{ii}=[tS.n_frac_samp, tS.n_subset];
end
clearvars tS;

% tidy vars structure
tau=cat(1,tau{:});
g2BB=cat(1,g2BB{:});
g2BB_samp=cat(1,g2BB_samp{:});
g2BB_se=cat(1,g2BB_se{:});
g2BB_mean=cat(1,g2BB_mean{:});
param_bs=cat(1,param_bs{:});

% sort
[tau,Isort]=sort(tau);
g2BB=g2BB(Isort,:);
g2BB_samp=g2BB_samp(Isort,:);
g2BB_se=g2BB_se(Isort,:);
g2BB_mean=g2BB_mean(Isort,:);

%% vis
% configs
[c,cl,cd]=palette(3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;

%% g2BB
h=figure('Name','g2BB',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
hold on;
pleg=NaN(1,3);
for ii=1:3
    tp=myploterr(1e6*tau,g2BB(:,ii),[],g2BB_se(:,ii),mark_typ{ii},c(ii,:));
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
    set(tp(2),'LineWidth',line_wid,'DisplayName','');
end
ax=gca;
lgd=legend(pleg);
lgd.FontSize=12;
box on;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('$g^{(2)}_\textrm{BB}(0)$');
ax.FontSize=12;
ax.LineWidth=1.2;

%% UNCERTAINTY g2BB
h=figure('Name','g2BB_uncertainty',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
hold on;
pleg=NaN(1,3);
for ii=1:3
    tp=plot(1e6*tau,g2BB_se(:,ii)./g2BB(:,ii),'Color','none',...
        'Marker',mark_typ{ii},'MarkerEdgeColor',c(ii,:),...
        'MarkerFaceColor',cl(ii,:),'LineWidth',line_wid,...
        'MarkerSize',mark_siz,'DisplayName',str_ss{ii});
    pleg(ii)=tp;     % line data to show in legend
end
ax=gca;
lgd=legend(pleg);
lgd.FontSize=12;
box on;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Frac. unc. of $g^{(2)}_\textrm{BB}(0)$');
ax.FontSize=12;
ax.LineWidth=1.2;