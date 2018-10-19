%% asymmetric detector saturation effect from MCP charge relaxation
% DKS
% 20181019

% study by looking for any sign of relaxation in detection efficiency after
% BEC transient

%% CONFIG
flag_save_data=false;

path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';

%%% vis
lim_phi_plot=asin(0.9*[-1,1]);      % z-limits to phi-lims

% get data (.mat) files to load
D=dir(path_dir);
tfnames={D.name};
b_isfile=cellfun(@(s) is_file(fullfile(path_dir,s)), tfnames);
filename=tfnames(b_isfile);
b_ismat=cellfun(@(s) strcmp(fileext(s),'.mat'),filename);
dataname=filename(b_ismat);
n_data=numel(dataname);

%% main
tau=cell(n_data,1);
ks=cell(n_data,1);

n_phi=cell(n_data,1);
ed_phi=cell(n_data,1);
phi=cell(n_data,1);

for ii=1:n_data
    %% load data
    tdataname=dataname{ii};
    tS=load(fullfile(path_dir,tdataname));
    n_shot_par=shotSize(tS.k_par);
    tnparam=tS.nparam;
    
    %% routine
    % general
    ttau=tS.par_T;
    tau{ii}=tS.par_T;
    
    tk_par=tS.k_par;
    % collate
    tk=cellfun(@(c) collate_shots(c),tk_par,'UniformOutput',false);
    tk=cat(1,tk{:});
    % to sph-pol coord
    tks=cellfun(@(z) zxy2sphpol(z),tk,'UniformOutput',false); 
    ks{ii}=tks;
    
    % histogram
    [tN_phi,ted_phi]=cellfun(@(s) histcounts(s(:,2)),tks,'UniformOutput',false);    % second col of spol is elev
    tphi=cellfun(@(e) e(1:end-1)+0.5*diff(e),ted_phi,'UniformOutput',false);
    tdphi=cellfun(@(p) p(2)-p(1),tphi,'UniformOutput',false);
    tdV=cellfun(@(p,dp) 2*pi*dp*cos(p),tphi,tdphi,'UniformOutput',false);    % diff area (vol) of dphi section
    tn_phi=cellfun(@(n,v) n./v,tN_phi,tdV,'UniformOutput',false);
    % normalise
    tnorm_phi=cellfun(@(n,dx) sum(n)*dx,tn_phi,tdphi,'UniformOutput',false);      % normalisation for PDF
    tn_phi=cellfun(@(x,c) x/c,tn_phi,tnorm_phi,'UniformOutput',false);
    % store
    n_phi{ii}=tn_phi;
    ed_phi{ii}=ted_phi;
    phi{ii}=tphi;
end
clearvars tS;

% tidy
tau=cat(1,tau{:});
ks=cat(1,ks{:});
n_tau=numel(tau);

n_phi=cat(1,n_phi{:});
ed_phi=cat(1,ed_phi{:});
phi=cat(1,phi{:});

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
facealpha=0.6;

%% elev angle histogram
H=cell(n_tau,1);
for ii=1:n_tau
    %%
    tfigname=sprintf('phi_%0.2g',1e6*tau(ii));
    H{ii}=figure('Name',tfigname,...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');
    hold on;
    thist=cell(2,1);
    for jj=1:2
%         thist{jj}=histogram(ks{ii,jj}(:,2),...
%             'Normalization','pdf',...
%             'DisplayName',str_ss{jj});
        thist{jj}=bar(phi{ii,jj},n_phi{ii,jj},'hist');
        thist{jj}.FaceColor=c(jj,:);
        thist{jj}.FaceAlpha=facealpha;
        
        thist{jj}.DisplayName=str_ss{jj};
    end
    ax=gca;
    box on;
    xlabel('Polar angle $\phi$ [rad]');
    ylabel('Norm. polar dist. $\rho(\phi)$');
    ax.FontSize=fontsize;
    ax.LineWidth=1.2;
    lgd=legend([thist{:}]);
    lgd.FontSize=fontsize;
    xlim(lim_phi_plot);
    
    str_tau=sprintf('%s %0.2g %s','$\tau=$',1e6*tau(ii),'$\mu$s');
    text('Units','normalized','Position',[0.5 1.03],'String',str_tau,...
        'FontSize',fontsize+2,'VerticalAlignment','baseline','HorizontalAlignment','center');
end


%% save data
if exist('flag_save_data','var') && flag_save_data
    %% manually run this section
    path_save_dir='C:\Users\HE BEC\Dropbox\phd\thesis\projects\bell_epr_steering\figS\det_qe';
    
    % workspace
    savefname=sprintf('supp_asymqe_%s',getdatetimestr);
    save(fullfile(path_save_dir,savefname));
    
    % phi dist
    for ii=1:n_tau
        savefigname=sprintf('fig_nphi_%0.2g_%s',1e6*tau(ii),getdatetimestr);
        fpath=fullfile(path_save_dir,savefigname);
        
        saveas(H{ii},strcat(fpath,'.fig'),'fig');
        print(H{ii},strcat(fpath,'.svg'),'-dsvg');
    end
end