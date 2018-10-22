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
ptch_col_bogo=0.85*ones(1,3);

% QENT: quantum entanglement
B_qent_max_pkpk=1;      % max pk-pk range of the diag correlator B
ptch_col_qent=[0.92,0.92,1.0];

%%% vis
font_siz_reg=13;
font_siz_sml=11;
font_siz_lrg=15;
mark_siz=9;
line_wid=2;

% ptch_col=0.85;         	% some gray to patch violation zone
ptch_alp=1;       		% patch transparency

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


%% VIS: g2 vs theta
[c0,clight,cdark]=palette(3);
str_ss={'$(1,1)$','$(0,0)$','$(1,0)$'};
mark_ss={'o','s','^'};

h=figure('Name','g2_vs_theta');

p=[];
for ii=1:ndatasets
    hold on;
    %     p(1)=plot(S(ii).par_T*om_rabi,S(ii).g2anti_par,'o',...
    %         'MarkerEdgeColor',c0(1,:),'MarkerFaceColor',clight(1,:),...
    %         'DisplayName','anti');
    %     p(2)=plot(S(ii).par_T*om_rabi,S(ii).g2corr_par,'^',...
    %         'MarkerEdgeColor',c0(2,:),'MarkerFaceColor',clight(2,:),...
    %         'DisplayName','corr');
    
    % spin-spin pairing
%     for jj=1:3
%         p(jj)=plot(S(ii).par_T*om_rabi,g2_0{ii}(:,jj),mark_ss{jj},...
%             'MarkerEdgeColor',c0(jj,:),'MarkerFaceColor',clight(jj,:),...
%             'DisplayName',str_ss{jj});
%     end
    
%     % mode-volume corrected
%     for jj=1:3
%         p(jj)=plot(S(ii).par_T*om_rabi,G2{ii}(:,jj),mark_ss{jj},...
%             'MarkerEdgeColor',c0(jj,:),'MarkerFaceColor',clight(jj,:),...
%             'DisplayName',str_ss{jj});
%     end

    % normalised
    for jj=1:3
        p(jj)=plot(S(ii).par_T*om_rabi,g2_norm{ii}(:,jj),mark_ss{jj},...
            'MarkerEdgeColor',c0(jj,:),'MarkerFaceColor',clight(jj,:),...
            'DisplayName',str_ss{jj});
    end
end

% annotation
lgd=legend(p);

box on;
xlabel('$\theta$');
ylabel('$\mathcal{G}$');

xlim(lim_th);
xticks(0:pi/4:pi);
xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});

set(gca,'FontSize',font_siz_reg);


%% VIS: Bell corr vs theta
cB=[0.5 0.3 0.7];       % main color for correlation
[cBlight,cBdark]=colshades(cB);      % get shades

h=figure('Name','B_vs_theta');

% TODO: collate and plot as whole
%   * use B_col rather than parsing through individual exp (S)
%   * integrate with the smoothed plot below
theta_col=om_rabi*cat(1,S.par_T);
B_col=cat(1,S.E_par);       % collate correlator from all experiments

%%% data
for ii=1:ndatasets
    hold on;
    p=ploterr(S(ii).par_T*om_rabi,S(ii).E_par,...
        [],S(ii).E_bootstrap_sdev,'o','hhxy',0);
        %         'MarkerEdgeColor',c0(1,:),'MarkerFaceColor',clight(1,:));
    
    set(p(1),'Marker','o','MarkerSize',mark_siz,...
        'Color',cB,'LineWidth',line_wid,...
        'MarkerFaceColor',cBlight(1,:),...
        'DisplayName','Data');
    set(p(2),'Color',cB,'LineWidth',line_wid,...
        'DisplayName','');
%     set(p(3),'Color',c0(1,:),'LineWidth',line_wid,...
%         'DisplayName','');
end

% smoothed fit to data
hold on;
polyfit_B=polyfit(theta_col,B_col,6);
theta_sm_fit=linspace(lim_th(1),lim_th(2),1e4);
B_sm_fit=polyval(polyfit_B,theta_sm_fit);
p_B_sm=plot(theta_sm_fit,B_sm_fit,'LineWidth',2.2,'Color',0.5*(cB+cBlight(1,:)));
uistack(p_B_sm,'bottom');


%%% Theory
% Bell+: Ideal rotation
th=linspace(lim_th(1),lim_th(2),1e3);
B_th_ideal=-cos(2*th);
p_B_th_ideal=plot(th,B_th_ideal,'Color',0.2*ones(1,3),'LineStyle','- -','LineWidth',2.2,...
    'DisplayName','$\vert\Psi^+\rangle$');
%     'DisplayName','Ideal $\vert\Psi^+\rangle$');
uistack(p_B_th_ideal,'bottom');
% text(3/4*pi,-0.75,'$\vert\Psi^+\rangle$','HorizontalAlignment','left',...
%     'FontSize',font_siz_lrg);

% annotation
ax=gca;
box on;
% axis square;
xlabel('Rotation angle $\theta$');
ylabel('Spin correlator $\mathcal{B}(\theta)$');

xlim(lim_th);
xticks(0:pi/4:pi);
xticklabels({'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'});
ylim([-1,1]);
% yticks(-1:0.5:1);

set(gca,'FontSize',font_siz_reg);
set(gca,'LineWidth',1.1);

%%% inequality region
%%%% BOGO
p_bogo=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [B_bogo_max,B_bogo_max,ax.YLim(2),ax.YLim(2)],...
    ptch_col_bogo,'FaceAlpha',ptch_alp,...
    'EdgeColor','none');
uistack(p_bogo,'bottom');      % this should REALLY be bottom - to not cover any other graphics
set(gca,'Layer','Top');     % graphics axes should be always on top
text(0,B_bogo_max,sprintf('Nonlocal\n(based on QM)'),'FontSize',font_siz_reg,'VerticalAlignment','bottom');

% %%%% QENT
% B_cent=0.5*(max(B_col)+min(B_col));   % centre value
% p_qent=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
%     B_cent*ones(1,4)+(0.5*B_qent_max_pkpk)*[-1,-1,+1,+1],...
%     ptch_col_qent,'FaceAlpha',ptch_alp,...
%     'EdgeColor','none');
% uistack(p_qent,'bottom');      % this should REALLY be bottom - to not cover any other graphics
% set(gca,'Layer','Top');     % graphics axes should be always on top
% text(0,B_cent+0.5*B_qent_max_pkpk,sprintf('Separable'),'FontSize',font_siz_reg,'VerticalAlignment','top');

% legend
lgd=legend([p(1), p_B_th_ideal]);

%%% misc
h.Renderer='painters';

%% EPR-steering parameter
% 20181022

tau=cat(1,S.par_T);
Th=tau*om_rabi;

B=cat(1,S.E_par);
B_se=cat(1,S.E_bootstrap_sdev);

% sort
[tau,Isort]=sort(tau);
Th=Th(Isort);
B=B(Isort);
B_se=B_se(Isort);

%%% EPR-steering
Th_epr=Th(1:9);
S_epr=abs(B(1:9)-B(9:end));
S_epr_se=vnorm(cat(2,B_se(1:9),B_se(9:end)),2);

% get fitted correlator
Th_fit=theta_sm_fit;
B_fit=B_sm_fit;

dth=Th_fit(2)-Th_fit(1);
dI=round((pi/2)/dth);

Th_epr_fit=Th_fit(1:end-dI+1);
S_epr_fit=abs(B_fit(1:end-dI+1) - B_fit(dI:end));



%% vis
% configs
[c,cl,cd]=palette(3);
% c_gray=0.6*ones(1,3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ptch_col_epr=0.85*ones(1,3);

% bounds
S_epr_min=sqrt(2);
S_ent_min=1;

% figure
h_epr=figure('Name','epr_steering',...
    'Units','normalized','Position',[0.2,0.2,0.2,0.3],'Renderer','painters');

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
ylabel('$\mathcal{S}\left(\theta\right)$')
ax.FontSize=fontsize;
ax.LineWidth=1.2;
xlim([-0.1,pi/2+0.1]);
ylim([0,2]);
xticks(0:pi/8:pi/2);
xticklabels({'$0$','$\pi/8$','$\pi/4$','$3\pi/8$','$\pi/2$'});

%%% violation region
% EPR
p_epr=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_epr_min,S_epr_min,2,2],...
    ptch_col_epr,'EdgeColor','none');
uistack(p_epr,'bottom');
% hatch fill to distinguish
hatchfill2(p_epr,'single','HatchAngle',45,'HatchDensity',30,...
    'HatchColor',0.25*ones(1,3),'HatchLineWidth',1);
set(gca,'Layer','Top');     % graphics axes should be always on top
text(pi/2,0.5*(S_epr_min+2),sprintf('EPR-steering'),'FontSize',font_siz_reg-2,...
    'HorizontalAlignment','right','VerticalAlignment','middle',...
    'BackgroundColor',ptch_col_epr);

% entanglement
p_ent=patch([ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1)],...
    [S_ent_min,S_ent_min,2,2],...
    ptch_col_epr,'EdgeColor','none');
uistack(p_ent,'bottom');   
text(pi/2,0.5*(S_ent_min+S_epr_min),sprintf('Entanglement'),'FontSize',font_siz_reg-2,...
    'HorizontalAlignment','right','VerticalAlignment','middle');

