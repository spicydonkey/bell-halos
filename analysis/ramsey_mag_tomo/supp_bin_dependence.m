% Supplemetary: halo metrology
% investigate spatial integration bin dependence on performance
%
% DKS
% 2019-02-19

%% config
% path to .mat file to load preprocessed halo data
path_mat_file='C:\Users\HE BEC\Dropbox\phd\projects\halo_metrology\analysis\ramsey\out\out_preproc_20190218_120000';

% spatial integration bin -------------------------------------------------
% bin size
c_alpha=logspace(log10(0.01),log10(1),200);       % factors of alpha
alpha=pi*c_alpha;

beta_ramsey=pi*0.0619;   % PSF rms width in angle (rad)


% sampling locations
az=pi*[0,1];
el=pi/2*[0];

% Phys-constants
Cphys=physConsts;
C_gymag=Cphys.He_gymag;         % He* gyromagnetic ratio (gamma) [Hz/G]

% figure
config_fig=loadFigureConfig;

%% load data
Sdata=load(path_mat_file);

% use halo data incl BECs
% k = Sdata.k_all_ramsey_phi0;
k = Sdata.k_halo_ramsey_phi0;

n_shot=cellfun(@(x) size(x,1),k);       % number of shots per config

tau = Sdata.tau;
tt=arrayfun(@(x,n) x*ones(n,1),tau,n_shot,'UniformOutput',false);
tt=vertcat(tt{:});

%% Ramsey model
Pramsey_mdl='y~amp*cos(om*x+phi)';
Pramsey_cname={'amp','om','phi'};
Pramsey_par0=[0.5,2*pi*1.5,0];

idx_om_Pramsey=2;

%% Analysis
% az-el location to sample
[gaz,gel]=ndgrid(az,el);
n_az=size(gaz,1);
n_el=size(gaz,2);
n_zone=numel(gaz);

% init
n_alpha=length(alpha);
B_alpha=NaN(n_az,n_el,n_alpha);
Berr_alpha=NaN(n_az,n_el,n_alpha);


% loop analysis for different bin size ------------------------------------
for kk=1:n_alpha
    talpha=alpha(kk);
    
    % get atom number ------------------------------------
    Nm_k=cell(length(k),1);     % #atoms in k,mj-zone categorise by exp param
    for ii=1:numel(k)
        % bin shots atom numbers in bin
        tk=k{ii};        % exp data for this param
        
        % get num in zone/shot
        tN=arrayfun(@(th,ph) cellfun(@(x) size(inCone(x,th,ph,talpha),1),tk),...
            gaz,gel,'UniformOutput',false);
        
        % format into multi-dim array: AZ X EL X SHOT X MJ
        ttN=NaN(cat(2,size(gaz),size(tN{1})));
        for jj=1:n_zone
            [iaz,iel]=ind2sub([n_az,n_el],jj);
            
            ttN(iaz,iel,:,:)=tN{iaz,iel};
        end
        Nm_k{ii}=ttN;
    end
    Ninv_k=cellfun(@(x) -diff(x,1,4),Nm_k,'UniformOutput',false);   % num-diff
    N_k=cellfun(@(x) sum(x,4),Nm_k,'UniformOutput',false);          % total number

    % polarisation ----------------------------------------
    P_k=cellfun(@(n,N) n./N,Ninv_k,N_k,'UniformOutput',false);
    P_k_avg=cellfun(@(x) squeeze(mean(x,3,'omitnan')),P_k,'UniformOutput',false);
    P_k_std=cellfun(@(x) squeeze(std(x,0,3,'omitnan')),P_k,'UniformOutput',false);
    % multidim array: T_DELAY X AZ X EL
    P_k_avg=cell2mat(cellfun(@(a) reshape(a,[1,size(a)]),P_k_avg,'UniformOutput',false));
    P_k_std=cell2mat(cellfun(@(a) reshape(a,[1,size(a)]),P_k_std,'UniformOutput',false));
    P_k_se=P_k_std./sqrt(n_shot);
    
    % fit Ramsey model -------------------------------------
    Pramseyk_fit=cell(n_az,n_el);   % dim: AZ X EL
    for ii=1:n_zone
        [iaz,iel]=ind2sub([n_az,n_el],ii);
        
        % collate data to fit
        yy=cellfun(@(x) squeeze(x(iaz,iel,:)),P_k,'UniformOutput',false);
        yy=vertcat(yy{:});
        
        % fit
        Pramseyk_fit{iaz,iel}=fitnlm(1e6*tt,yy,Pramsey_mdl,...
            Pramsey_par0,'CoefficientNames',Pramsey_cname);
    end
    
    % get fit params
    Pramseyk_fpar=arrayfun(@(I) cellfun(@(f) f.Coefficients.Estimate(I),Pramseyk_fit),...
        1:numel(Pramsey_cname),'UniformOutput',false);
    Pramseyk_fpar=cat(ndims(Pramseyk_fpar{1})+1,Pramseyk_fpar{:});
    
    Pramseyk_fparerr=arrayfun(@(I) cellfun(@(f) f.Coefficients.SE(I),Pramseyk_fit),...
        1:numel(Pramsey_cname),'UniformOutput',false);
    Pramseyk_fparerr=cat(ndims(Pramseyk_fparerr{1})+1,Pramseyk_fparerr{:});
    % dim: AZ X EL X PAR#
    
    % get freq
    omk_Pramsey=Pramseyk_fpar(:,:,idx_om_Pramsey);
    omerrk_Pramsey=Pramseyk_fparerr(:,:,idx_om_Pramsey);
    % dim: AZ X EL

    % magnetic field -----------------------------------------
    Bk_Pramsey=1e6*omk_Pramsey/(2*pi*C_gymag);
    Berrk_Pramsey=1e6*omerrk_Pramsey/(2*pi*C_gymag);
    
    % store 
    B_alpha(:,:,kk)=Bk_Pramsey;
    Berr_alpha(:,:,kk)=Berrk_Pramsey;
end

%% VIS: binsize vs mean
% config
col={'b','k'};      % manual colors

% plot
h=figure('Name','binsize_vs_B','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
ax=gca;
hold on;

pleg=[];
counter=1;
for ii=1:n_az
    for jj=1:n_el
        tI=sub2ind([n_az,n_el],ii,jj);

        tp=shadedErrorBar(alpha/pi,squeeze(B_alpha(ii,jj,:)),squeeze(Berr_alpha(ii,jj,:)));
        tp.mainLine.LineStyle=config_fig.line_sty{tI};
        tp.mainLine.DisplayName=sprintf('%0.3g, %0.3g',az(ii),el(jj));
        tp.patch.FaceAlpha=0.2;
        
        % color
        tp.patch.FaceColor=col{counter};
        tp.mainLine.Color=col{counter};
        set(tp.edge,'Color',col{counter});
        
        pleg(end+1)=tp.mainLine;
        
        counter=counter+1;
    end
end

set(ax,'Layer','top');
set(ax,'XScale','log');
xlabel('bin size $\alpha/\pi$');
ylabel('$B$ (G)');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

% % annotation ----------------------
% line_resol=line([beta_ramsey,beta_ramsey]/pi,ax.YLim,'LineStyle','--','Color','k');
% uistack(line_resol,'bottom');


% legend -------------------------------
lgd=legend(pleg);
lgd.Title.String='$\theta, \phi$';
lgd.Location='best';
lgd.Box='off';


%% VIS: binsize vs error
h=figure('Name','binsize_vs_B_err','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
ax=gca;
hold on;
counter=1;
for ii=1:n_az
    for jj=1:n_el
        tI=sub2ind([n_az,n_el],ii,jj);
%         tp=plot(alpha/pi,squeeze(Berr_alpha(ii,jj,:)));
        tp=plot(alpha/pi,squeeze(Berr_alpha(ii,jj,:)));
        tp.LineStyle=config_fig.line_sty{tI};
        tp.Marker='none';
%         tp.Color='k';
        tp.Color=col{counter};
    
        counter=counter+1;
    end
    
end
set(ax,'Layer','top');
set(ax,'XScale','log');
set(ax,'YScale','log');
xlabel('bin size $\alpha/\pi$');
ylabel('$\Delta B$ (G)');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

% % annotation ----------------------
% line_resol=line([beta_ramsey,beta_ramsey]/pi,ax.YLim,'LineStyle','--','Color','k');
% uistack(line_resol,'bottom');


%% VIS: binsize vs error (scaling by half-cone angle)
% %   NOTE: Vol \propto alpha^2 relation will break down for large alpha
% h=figure('Name','binsize_vs_B_err','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
% ax=gca;
% hold on;
% counter=1;
% for ii=1:n_az
%     for jj=1:n_el
%         tI=sub2ind([n_az,n_el],ii,jj);
% 
%         y_scaled = (alpha/pi)'.*squeeze(Berr_alpha(ii,jj,:));
%         y_scaled = y_scaled/y_scaled(1);
%         
%         tp=plot(alpha/pi,y_scaled);
%         tp.LineStyle=config_fig.line_sty{tI};
%         tp.Marker='none';
%         tp.Color=col{counter};
%         
%         counter=counter+1;
%     end
% end
% set(ax,'Layer','top');
% set(ax,'XScale','log');
% set(ax,'YScale','log');
% xlabel('bin size $\alpha/\pi$');
% ylabel('$\alpha\cdot\Delta B$ (arb. u.)');
% box on;
% ax.FontSize=config_fig.ax_fontsize;
% ax.LineWidth=config_fig.ax_lwid;
% 
% % % annotation ----------------------
% % line_resol=line([beta_ramsey,beta_ramsey]/pi,ax.YLim,'LineStyle','--','Color','k');
% % uistack(line_resol,'bottom');


%% VIS: binsize vs error (scaling by sqrt(cone solid-angle))
%   NOTE:   Vol \propto solidangle(cone)
%           N \propto Vol
       
h=figure('Name','binsize_vs_B_err','Units',config_fig.units,'Position',[0,0,8.6,4.5],'Renderer',config_fig.rend);
ax=gca;
hold on;
counter=1;
for ii=1:n_az
    for jj=1:n_el
        tI=sub2ind([n_az,n_el],ii,jj);

        y_scaled = sqrt(cone_solang(alpha))'.*squeeze(Berr_alpha(ii,jj,:));
        y_scaled = y_scaled/y_scaled(1);
        
        tp=plot(alpha/pi,y_scaled);
        tp.LineStyle=config_fig.line_sty{tI};
        tp.Marker='none';
        tp.Color=col{counter};
        
        counter=counter+1;
    end
end
set(ax,'Layer','top');
set(ax,'XScale','log');
set(ax,'YScale','log');
xlabel('bin size $\alpha/\pi$');
ylabel('$\sqrt{\delta V} \cdot\Delta B$ (arb. u.)');
box on;
ax.FontSize=config_fig.ax_fontsize;
ax.LineWidth=config_fig.ax_lwid;

% % annotation ----------------------
% line_resol=line([beta_ramsey,beta_ramsey]/pi,ax.YLim,'LineStyle','--','Color','k');
% uistack(line_resol,'bottom');


%% save ===============================================
vars_to_save={'path_mat_file','config_fig',...
    'alpha','az','el','n_az','n_el',...
    'B_alpha','Berr_alpha',...
    };

% check exists
for ii=1:numel(vars_to_save)
    tvarname=vars_to_save{ii};
    if ~exist(tvarname,'var')
        warning('variable %s does not exist.',tvarname);
    end
end

% save
save(['supp_binsize_',getdatetimestr,'.mat'],vars_to_save{:});

