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
% load('C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot\bell_signal\bell_20180815_1.mat');

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