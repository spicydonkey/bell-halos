%% Scattered number in halos
% DKS
% 20181018

%% CONFIG
path_dir='C:\Users\David\Dropbox\PhD\data\bell_epr_2018\proc\exp5_bell_yrot';

% get data (.mat) files to load
D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_ismat=cellfun(@(s) strcmp(fileext(s),'.mat'),filename);
dataname=filename(b_ismat);
n_data=numel(dataname);

%%% analysis
det_qe=0.1;

%% main
% preallocate
tau=cell(n_data,1);
Fcorr=cell(n_data,1);
Ndet_Vsc=cell(n_data,1);
Ndet_Vsc_mean=cell(n_data,1);
Ndet_Vsc_std=cell(n_data,1);
Ndet_Vsc_se=cell(n_data,1);
Nsc_mean=cell(n_data,1);
Nsc_std=cell(n_data,1);
Nsc_se=cell(n_data,1);

for ii=1:n_data
    %% load data
    tdataname=dataname{ii};
    tS=load(fullfile(path_dir,tdataname));
    n_shot_par=shotSize(tS.k_par);
    
    %% routine
    % general
    tau{ii}=tS.par_T;
    % get number in captured halo volume
    tNdet_Vsc=cellfun(@(c) cellfun(@(k) size(k,1),c),tS.k_par,'UniformOutput',false)';
    % stats
    tNdet_Vsc_mean=cellfun(@(n) mean(n,1),tNdet_Vsc,'UniformOutput',false);
    tNdet_Vsc_mean=cat(1,tNdet_Vsc_mean{:});
    tNdet_Vsc_std=cellfun(@(n) std(n,0,1),tNdet_Vsc,'UniformOutput',false);
    tNdet_Vsc_std=cat(1,tNdet_Vsc_std{:});
    tNdet_Vsc_se=tNdet_Vsc_std./sqrt(n_shot_par)';
    % total scattered number
    tVsc_geo_factor=1;   % TODO geometric factor to correct for complete halo
    tFcorr=(1/det_qe)*tVsc_geo_factor;    % complete corection factor
    
    tNsc_mean=tFcorr*tNdet_Vsc_mean;
    tNsc_std=tFcorr*tNdet_Vsc_std;
    tNsc_se=tFcorr*tNdet_Vsc_se;
    
    % store vars
    Fcorr{ii}=tFcorr;
    Ndet_Vsc{ii}=tNdet_Vsc;
    Ndet_Vsc_mean{ii}=tNdet_Vsc_mean;
    Ndet_Vsc_std{ii}=tNdet_Vsc_std;
    Ndet_Vsc_se{ii}=tNdet_Vsc_se;
    Nsc_mean{ii}=tNsc_mean;
    Nsc_std{ii}=tNsc_std;
    Nsc_se{ii}=tNsc_se;
    
end
clearvars tS;

% tidy var structure
tau=cat(1,tau{:});
Fcorr=cat(1,Fcorr{:});
Ndet_Vsc=cat(1,Ndet_Vsc{:});
Ndet_Vsc_mean=cat(1,Ndet_Vsc_mean{:});
Ndet_Vsc_std=cat(1,Ndet_Vsc_std{:});
Ndet_Vsc_se=cat(1,Ndet_Vsc_se{:});
Nsc_mean=cat(1,Nsc_mean{:});
Nsc_std=cat(1,Nsc_std{:});
Nsc_se=cat(1,Nsc_se{:});

% sort
[tau,Isort]=sort(tau);
Ndet_Vsc=Ndet_Vsc(Isort,:);
Ndet_Vsc_mean=Ndet_Vsc_mean(Isort,:);
Ndet_Vsc_std=Ndet_Vsc_std(Isort,:);
Ndet_Vsc_se=Ndet_Vsc_se(Isort,:);
Nsc_mean=Nsc_mean(Isort,:);
Nsc_std=Nsc_std(Isort,:);
Nsc_se=Nsc_se(Isort,:);

n_tau=numel(tau);

%% vis
% configs
[c,cl,cd]=palette(3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;

%% total scattered number 
h=figure('Name','scattered_number',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
hold on;
pleg=NaN(1,2);
for ii=1:2
    tp=myploterr(1e6*tau,Nsc_mean(:,ii),[],Nsc_se(:,ii),mark_typ{ii},c(ii,:));
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
    set(tp(2),'LineWidth',line_wid,'DisplayName','');
end
ax=gca;
lgd=legend(pleg);
lgd.FontSize=fontsize;
box on;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('$N_{\textrm{sc}}$');
ax.FontSize=fontsize;
ax.LineWidth=1.2;
tylim=ax.YLim;
ylim([0,tylim(2)]);

%% detected number distribution
H=NaN(n_tau,1);
for ii=1:n_tau
    %%
    tfigname=sprintf('Ndet_dist_%0.2g',1e6*tau(ii));
    H(ii)=figure('Name',tfigname,...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
    hold on;
    pleg=NaN(2,1);
    for jj=1:2
        pleg(jj)=histogram(Ndet_Vsc{ii}(:,jj),'DisplayName',str_ss{jj});
    end
    ax=gca;
    box on;
    xlabel('$N_{\textrm{det}}$');
    ylabel('$\#$ observations');
    ax.FontSize=fontsize;
    ax.LineWidth=1.2;
    lgd=legend(pleg);
    lgd.FontSize=fontsize;
end