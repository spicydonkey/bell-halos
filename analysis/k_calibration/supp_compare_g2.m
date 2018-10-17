%% SUPP: Compare g2 correlation values
% (1) g2 obtained across BB regions in calibrated "k"-space
% (2) amplitude of fitted g2 function
%
% DKS
% 20181017

%% configs
dk_bb=0.04;     % BB corr length [norm unit]

path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';


%% load data
D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_ismat=cellfun(@(s) strcmp(fileext(s),'.mat'),filename);
dataname=filename(b_ismat);
n_data=numel(dataname);

S=cell(n_data,1);
for ii=1:n_data
    % get data
    S{ii}=load(fullfile(path_dir,dataname{ii}),'par_T','g2','g2mdl','nparam','idx_dk0');
end

%% get g2
% CELL as temporary var holder
tau=cell(n_data,1);
g2_BB_0=cell(n_data,1);
g2_amp=cell(n_data,1);
g2_param=cell(n_data,1);
g_ratio=cell(n_data,1);

for ii=1:n_data
    % load vars from struct
    tS=S{ii};
    par_T=tS.par_T;
    g2=tS.g2;
    g2mdl=tS.g2mdl;
    nparam=tS.nparam;
    idx_dk0=tS.idx_dk0;
    
    %%% evaluate
    tau{ii}=par_T;
    
    % (1)
    tg2_BB_0=NaN(nparam,3);
    for jj=1:nparam
        tg2_BB_0(jj,:)=cellfun(@(g) g(idx_dk0,idx_dk0,idx_dk0),g2{jj});
    end
    
    % (2)
    tg2_amp=NaN(nparam,3);
    tg2_param=cell(nparam,3);
    for jj=1:nparam
        % fitted model params
        tg2param=cellfun(@(m) m.Coefficients.Estimate, g2mdl{jj},'UniformOutput',false);
        % check fit ok
        for kk=1:3      % corr type
            ttg2param=tg2param{kk};
            tg2_param{jj,kk}=ttg2param;
            tamp=ttg2param(1);
            tmu=ttg2param(2:4);
            tsig=ttg2param(5:7);
            if ~(vnorm(tmu')>2*dk_bb || geomean(abs(tsig))>2*dk_bb)
                % passes a simple test
                tg2_amp(jj,kk)=tamp;
            end
        end
    end
    % Comparison
    tg_ratio=tg2_amp./tg2_BB_0;
    
    %%% store results
    g2_BB_0{ii}=tg2_BB_0;
    g2_amp{ii}=tg2_amp;
    g_ratio{ii}=tg_ratio;
    g2_param{ii}=tg2_param;
end
% TIDY var structure: to array!
tau=cat(1,tau{:});
g2_BB_0=cat(1,g2_BB_0{:});
g2_amp=cat(1,g2_amp{:});
g2_param=cat(1,g2_param{:});
g_ratio=cat(1,g_ratio{:});

% sort vars in TAU
[tau,Isort]=sort(tau);
g2_BB_0=g2_BB_0(Isort,:);
g2_amp=g2_amp(Isort,:);
g2_param=g2_param(Isort,:);
g_ratio=g_ratio(Isort,:);

%%% processed - correlator
% sum correlations
G_BB=sum(g2_BB_0,2)+g2_BB_0(:,3);       % anti-corr twice!
G_amp=sum(g2_amp,2)+g2_amp(:,3);

%% vis
%%% configs
[c,cl,cd]=palette(3);
c_gray=0.6*ones(1,3);
% line_sty={'-','--',':'};
line_sty={'-','-','-'};
mark_typ={'o','s','^'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};

%%% g2 ratio
figname='g2_ratio';

h=figure('Name',figname,'Units','normalized',...
    'Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
hold on;
p=NaN(3,1);
for ii=1:3
    p(ii)=plot(1e6*tau,g_ratio(:,ii),...
        'Color',c(ii,:),'LineStyle',line_sty{ii},'LineWidth',1.5,...
        'Marker',mark_typ{ii},'MarkerFaceColor',cl(ii,:),'MarkerSize',7,...
        'DisplayName',str_ss{ii});
end
ax=gca;
box on;
xlabel('$\tau$ [$\mu$s]');
ylabel('Ratio of $g^{(2)}$ correlation');
ax.FontSize=12;
ax.LineWidth=1.2;
legend(p,'Location','southeast');

%% sum g2
figname='g2_sum';
h=figure('Name',figname,'Units','normalized',...
    'Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
hold on;
p=[];

% BB
p(1)=plot(1e6*tau,G_BB,...
    'Color','k','LineStyle','-','LineWidth',1.5,...
    'Marker','o','MarkerFaceColor',c_gray,'MarkerSize',7,...
    'DisplayName','BB');

% fitted amplitude
p(2)=plot(1e6*tau,G_amp,...
    'Color','k','LineStyle','--','LineWidth',1.5,...
    'Marker','d','MarkerFaceColor',c_gray,'MarkerSize',7,...
    'DisplayName','fit amp');

ax=gca;
box on;
xlabel('$\tau$ [$\mu$s]');
ylabel('$\sum g^{(2)}_{S^{A}S^{B}}$');
ax.FontSize=12;
ax.LineWidth=1.2;
legend(p,'Location','northeast');
