%% Integrated g2 correlation
%   to be free of correlation length/volume
% DKS
% 20181019

%% CONFIG
dk_intg=0.1;     % half-width around dk=0 for g2 integration (ensure it is in range of g2)

flag_save_data=false;

% path_dir='C:\Users\David\Dropbox\PhD\data\bell_epr_2018\proc\exp5_bell_yrot';
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';

% get data (.mat) files to load
D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_ismat=cellfun(@(s) strcmp(fileext(s),'.mat'),filename);
dataname=filename(b_ismat);
n_data=numel(dataname);

%% main
tau=cell(n_data,1);     % pulse dur (s)
G=cell(n_data,1);       % volume integrated correlation
G_se=cell(n_data,1);    % STDERR of G
E=cell(n_data,1);       % correlator eval'd from G
E_se=cell(n_data,1);    % err in E

for ii=1:n_data
    %% load data
    tdataname=dataname{ii};
    tS=load(fullfile(path_dir,tdataname));
    tidx0=tS.idx_dk0;
    tau{ii}=tS.par_T;
    tg2=tS.g2;
    
    %% routine
    tn_dk=numel(tS.dk_cent_vec);
    tdk_bin=tS.dk_cent_vec(2)-tS.dk_cent_vec(1);
    tIlim=round(dk_intg/tdk_bin);
    tidx_intg=tidx0+[-tIlim:tIlim];
    tV_intg=(tdk_bin*(2*tIlim+1))^3;    % integration vol (norm units)
    n_pix_int=numel(tidx_intg)^3;       % number of integrated pixels
    
    % integrated g2
    tG=cellfun(@(s) cellfun(@(g) meanall(g(tidx_intg,tidx_intg,tidx_intg)),s),tg2,'UniformOutput',false);
    tG=cat(1,tG{:});
    G{ii}=tG;
    
    % error (bootstrapping)
    tfrac_bs=tS.n_frac_samp;
    tg2_bs=tS.g2_bootstrap;
    tnsamp_bs=tS.n_subset;      % num samps repeated for bootstrapping
    tG_bs=cellfun(@(s) ...
        cellfun(@(c) ...
            cellfun(@(g) ...
                meanall(g(tidx_intg,tidx_intg,tidx_intg),'omitnan'),...     % integrated g2
            squeeze(mat2cell(c,tn_dk,tn_dk,tn_dk,ones(1,tnsamp_bs)))),...   % format BS # into cell array
        s,'UniformOutput',false),...
    tS.g2_bootstrap,'UniformOutput',false);
    tG_bs=cellfun(@(g) cat(2,g{:}),tG_bs,'UniformOutput',false);    % tidy format: array-cat corrtype
    tG_se=cellfun(@(x) sqrt(tfrac_bs)*std(x,[],1),tG_bs,'UniformOutput',false);    % std err of mean from bootstrapping
    tG_se=cat(1,tG_se{:});
    G_se{ii}=tG_se;     % store var
    
    % integrated correlator
    tE=Ecorr(tG(:,1),tG(:,2),tG(:,3),tG(:,3));
    E{ii}=tE;
    
    % intg corr from bs
    tE_bs=cellfun(@(g) Ecorr(g(:,1),g(:,2),g(:,3),g(:,3)), tG_bs,'UniformOutput',false);
    tE_se=cellfun(@(e) sqrt(tfrac_bs)*std(e),tE_bs);
    E_se{ii}=tE_se;
end
clearvars tS;

tau=cat(1,tau{:});
G=cat(1,G{:});
G_se=cat(1,G_se{:});
E=cat(1,E{:});
E_se=cat(1,E_se{:});

%% vis
% configs
[c,cl,cd]=palette(3);
c_gray=0.6*ones(1,3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;

%% Integrated g2
h=figure('Name','integrated_g2',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');

hold on;

pleg=NaN(3,1);
for ii=1:3
    tp=myploterr(1e6*tau,G(:,ii),[],G_se(:,ii),mark_typ{ii},c(ii,:));
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
    set(tp(2),'LineWidth',line_wid,'DisplayName','');
end
box on;
ax=gca;
% lgd=legend(pleg,'Location','southeast');
% lgd.FontSize=fontsize;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Integrated correlation $\bar{g}^{(2)}_{\mathrm{BB}}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
tylim=ax.YLim;
ylim([0,tylim(2)]);

%% Sum g2
h=figure('Name','sum_integrated_g2',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
G_sum=sum(G,2)+G(:,3);
G_sum_se=vnorm(cat(2,G_se,G_se(:,3)),2);
tp=myploterr(1e6*tau,G_sum,[],G_sum_se,'o',[0,0,0]);
set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'MarkerFaceColor',c_gray);
set(tp(2),'LineWidth',line_wid,'DisplayName','');

box on;
ax=gca;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('$\sum \bar{g}^{(2)}_{\mathrm{BB}}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;

%% correlator
h=figure('Name','correlator',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');

hold on;

tp=myploterr(1e6*tau,E,[],E_se,mark_typ{1},c(1,:));
set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName','Data');
pleg=tp(1);     % line data to show in legend
set(tp(2),'LineWidth',line_wid,'DisplayName','');

box on;
ax=gca;
% lgd=legend(pleg);
% lgd.FontSize=fontsize;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Integrated Correlator $\bar{\mathcal{B}}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
axis tight;
ylim([-1,1]);
