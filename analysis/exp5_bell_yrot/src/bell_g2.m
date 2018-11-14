%% summary of Bell signal and g2 under global rotation
%
% DKS
% 20180814

%% configs
%%% Rabi
Tpi=10e-6;        % pi-pulse duration [s]
om_rabi=2*pi*50.3e3;    % Rabi frequency [rad/s]

%%% quantum correlations
% BOGO: Bogoluibov theory --> violation of Bell inequality
B_bogo_max=1/sqrt(2);

% QENT: quantum entanglement
B_qent_max_pkpk=1;      % max pk-pk range of the diag correlator B

% steering param
S_epr_min=sqrt(2);
S_ent_min=1;

%%% VIS
% general
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_ren='painters';

[c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
% str_ss={'$(1,1)$','$(0,0)$','$(1,0)$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
ax_lwid=1.2;
fontsize=12;

% custom
cB=[0.5 0.3 0.7];               % main color for correlation
[cB_l,cB_d]=colshades(cB);      % get shades

c_epr=0.85*ones(1,3);
c_bogo=0.85*ones(1,3);
c_qent=[0.92,0.92,1.0];

% c_patch=0.85;         	% some gray to patch violation zone
% a_patch=1;       		% patch transparency

lim_th=[-pi/20,pi+pi/20];   % theta-axis limits

%% load collated data
load('C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot\bell_signal\bell_20181009_bs_fix.mat');
% load('C:\Users\David\Dropbox\PhD\data\bell_epr_2018\proc\exp5_bell_yrot\bell_signal\bell_20181009_bs_fix.mat');

% OR have dataset loaded as "S"


%% post processing 
ndatasets=numel(S);

%%% g2 profiles
g2_0=cell(ndatasets,1);     % g2 at dk=0
vk=cell(ndatasets,1);       % ~g2 dk-volume
G2=cell(ndatasets,1);       % g2 amplitude: mode-volume corrected 
for ii=1:ndatasets
    % preallocate 
    g2_0{ii}=NaN(length(S(ii).par_T),3);
    vk{ii}=NaN(length(S(ii).par_T),3);
    G2{ii}=NaN(length(S(ii).par_T),3);
    
    tidx_0=S(ii).idx_dk0;
    for jj=1:length(S(ii).par_T)
        for kk=1:3
            % datapoint at dk=0
            g2_0{ii}(jj,kk)=S(ii).g2{jj}{kk}(tidx_0,tidx_0,tidx_0);

%             % fitted amplitude
%             g2_0{ii}(jj,kk)=S(ii).g2mdl{jj}{kk}.Coefficients.Estimate(1)+1;
            
            % fitted corr-volume
            % NOTE: product of 3D gaussian sigmas
            vk{ii}(jj,kk)=prod(abs(S(ii).g2mdl{jj}{kk}.Coefficients.Estimate(4:7)));
            
        end
    end
    % mode-volume corrected
    G2{ii}=g2_0{ii}.*vk{ii};
end

%%% normalising g2
%   NOTE: total g2 needs to take into account that there is 01/10 pairing
%           normalisation is scaled s.t. (g2_11+g2_00)/2 + g2_10 = 1
g2_tot=cellfun(@(g) sum(g,2)+g(:,3),g2_0,'UniformOutput',false);
g2_norm=cellfun(@(g0,gtot) 2*g0./gtot,g2_0,g2_tot,'UniformOutput',false);

%% tidy vars
%%% concat as array (lose run# from cell struct)
tau=cat(1,S.par_T);     % pulse duration [s]
Th=tau*om_rabi;         % rotation angle

B=cat(1,S.E_par);       % B correlator
B_se=cat(1,S.E_bootstrap_sdev);     % stderr of B

g2_tot=cat(1,g2_tot{:});
g2_norm=cat(1,g2_norm{:});

%%% sort
[tau,Isort]=sort(tau);
Th=Th(Isort);
B=B(Isort);
B_se=B_se(Isort);
g2_tot=g2_tot(Isort,:);
g2_norm=g2_norm(Isort,:);

%%% clear raw data
clearvars S;    

%% VIS: Norm g2 vs theta
h=figure('Name','g2_vs_theta','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

p=[];
for jj=1:3
    p(jj)=plot(Th,g2_norm(:,jj),mark_typ{jj},...
        'MarkerEdgeColor',c(jj,:),'MarkerFaceColor',cl(jj,:),...
        'DisplayName',str_ss{jj});
end

% annotation
box on;

ax=gca;
ax.FontSize=fontsize;
ax.LineWidth=1.2;

lgd=legend(p);
lgd.FontSize=fontsize;

xlabel('Rotation angle $\theta$');
ylabel('Normalised correlation $\mathcal{G}$');

xlim(lim_th);
xticks(0:pi/4:pi);
xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});

%% VIS: Bell corr vs theta
h=figure('Name','B_vs_theta','Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

p=myploterr(Th,B,[],B_se,'o',cB);
set(p(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName','Data');
set(p(2),'LineWidth',line_wid,'DisplayName','');

% smoothed fit to data
hold on;
polyfit_B=polyfit(Th,B,6);
Th_fit=linspace(lim_th(1),lim_th(2),1e4);
B_fit=polyval(polyfit_B,Th_fit);
p_B_sm=plot(Th_fit,B_fit,'LineWidth',line_wid,'Color',0.5*(cB+cB_l(1,:)));
uistack(p_B_sm,'bottom');

%%% Theory
% Bell+: Ideal rotation
th=linspace(lim_th(1),lim_th(2),1e3);
B_th_psi=-cos(2*th);
p_B_psi=plot(th,B_th_psi,'Color',0.2*ones(1,3),'LineStyle','- -','LineWidth',line_wid,...
    'DisplayName','$\vert\Psi^+\rangle$');
uistack(p_B_psi,'bottom');

% annotation
box on;

ax=gca;
ax.FontSize=fontsize;
ax.LineWidth=ax_lwid;

xlabel('Rotation angle $\theta$');
ylabel('Correlator $\mathcal{B}(\theta)$');

xlim(lim_th);
xticks(0:pi/4:pi);
xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});
ylim([-1,1]);
% yticks(-1:0.5:1);

%%% inequality region
%%%% BOGO
% ABOVE
p_bogo=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [B_bogo_max,B_bogo_max,ax.YLim(2),ax.YLim(2)],...
    c_bogo,...
    'EdgeColor','none');
uistack(p_bogo,'bottom');      % this should REALLY be bottom - to not cover any other graphics

text(0,0.5*(B_bogo_max+1),sprintf('Nonlocal (QM)'),'FontSize',fontsize-1,'VerticalAlignment','middle');

% BELOW
p_bogo_2=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [ax.YLim(1),ax.YLim(1),-B_bogo_max,-B_bogo_max],...
    c_bogo,...
    'EdgeColor','none');
uistack(p_bogo_2,'bottom');      
set(gca,'Layer','Top');
% text(0,-0.5*(B_bogo_max+1),sprintf('Nonlocal (QM)'),'FontSize',fontsize-1,'VerticalAlignment','middle');

% %%%% QENT
% B_cent=0.5*(max(B)+min(B));   % centre value
% p_qent=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
%     B_cent*ones(1,4)+(0.5*B_qent_max_pkpk)*[-1,-1,+1,+1],...
%     ptch_col_qent,'FaceAlpha',ptch_alp,...
%     'EdgeColor','none');
% uistack(p_qent,'bottom');      % this should REALLY be bottom - to not cover any other graphics
% set(gca,'Layer','Top');     % graphics axes should be always on top
% text(0,B_cent+0.5*B_qent_max_pkpk,sprintf('Separable'),'FontSize',fontsize,'VerticalAlignment','top');

lgd=legend([p(1), p_B_psi]);
set(gca,'Layer','Top');     % graphics axes should be always on top

%% EPR-steering parameter
% 20181022

%%% EPR-steering
Th_epr=Th(1:9);
S_epr=abs(B(1:9)-B(9:end));
S_epr_se=vnorm(cat(2,B_se(1:9),B_se(9:end)),2);

% from correlator fit
dth=Th_fit(2)-Th_fit(1);
dI=round((pi/2)/dth);

Th_epr_fit=Th_fit(1:end-dI+1);
S_epr_fit=abs(B_fit(1:end-dI+1) - B_fit(dI:end));


%% vis: EPR-steering
h_epr=figure('Name','epr_steering',...
    'Units',f_units,'Position',f_pos,'Renderer',f_ren);
hold on;

% data
tp=myploterr(Th_epr,S_epr,[],S_epr_se,'o',c(1,:));
set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName','');
pleg=tp(1);     % line data to show in legend
set(tp(2),'LineWidth',line_wid,'DisplayName','');

% fit
tp=plot(Th_epr_fit,S_epr_fit,'Color',c(1,:),'LineWidth',line_wid,...
    'LineStyle','-');
uistack(tp,'bottom');

box on;
ax=gca;
xlabel('Rotation angle $\theta$');
ylabel('$\mathcal{S}\left(\theta,\theta+\pi/2\right)$')
ax.FontSize=fontsize;
ax.LineWidth=ax_lwid;
xlim([-0.1,pi/2+0.1]);
ylim([0,2]);
xticks(0:pi/8:pi/2);
xticklabels({'$0$','$\pi/8$','$\pi/4$','$3\pi/8$','$\pi/2$'});

%%% violation region
% EPR
p_epr=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_epr_min,S_epr_min,2,2],...
    c_epr,'EdgeColor','none');
uistack(p_epr,'bottom');
% hatch to distinguish
hatchfill2(p_epr,'single','HatchAngle',45,'HatchDensity',30,...
    'HatchColor','k','HatchLineWidth',1);
set(gca,'Layer','Top');     % graphics axes should be always on top
text(pi/4,0.5*(S_epr_min+2),sprintf('Bell correlation witness'),'FontSize',fontsize-1,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'BackgroundColor',c_epr);

% entanglement
p_ent=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_ent_min,S_ent_min,2,2],...
    c_epr,'EdgeColor','none');
uistack(p_ent,'bottom');   
text(pi/4,0.5*(S_ent_min+S_epr_min),sprintf('Entanglement'),'FontSize',fontsize-1,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

