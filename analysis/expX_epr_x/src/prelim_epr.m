%% Prelim multiplicative variance EPR-steering
%
% DKS
% 2018-11-19
%

%% CONFIGS
% fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\expX_epr_x\prelim_20181119\k_zz_20181119.mat';
% fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\expX_epr_x\prelim_20181119\k_xx_20181119.mat';
% fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\expX_epr_x\prelim_20181119\k_yy_20181119.mat';
fdata='C:\Users\HE BEC\Dropbox\phd\data\bell_epr_2018\proc\expX_epr_x\prelim_20181119\k_mm_20181119.mat';

% k-mode (A,B - spherically opposite)
alpha=0.08;

n_az=40;                	% equispaced bins
n_el=20;

az_disp=deg2rad([0,45,90]);     % azim sections (great circles) to display
el_disp=deg2rad([-30,0,30]);    % elev/lat zones to display

% EPR-steering
%   Note: symmetry from the post-selected spin-1/2 condition (RHS of
%   criteria)
g_inf = {-1, -1, 1};    % for Jx, Jy, Jz measurements

% vis
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';

cviridis=viridis(5);
[c,cl,cd]=palette(3);
c_gray=0.8*ones(1,3);
line_sty={'-','--',':','-.'};
mark_typ={'o','s','^','d'};
str_mm={'$x$','$y$','$z$'};
% str_ss={'$\vert\!\uparrow\rangle$','$\vert\!\downarrow\rangle$'};
% str_ss={'$m_J = 1$','$m_J = 0$'};
str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
mark_siz=7;
line_wid=1.5;
fontsize=12;
ax_lwidth=1.2;

%% load data
load(fdata);
% k_mm = {k_xx,k_yy,k_zz};

%% construct k-modes
lim_az=[0,pi];              % limit of azim angle (exc +pi)
phi_max=pi/4;
lim_el=[-phi_max,phi_max];

Vaz=linspace(lim_az(1),lim_az(2),n_az+1);
Vaz=Vaz(1:end-1);       % exclude last
Vel=linspace(lim_el(1),lim_el(2),n_el);

[vaz,vel]=ndgrid(Vaz,Vel);    % AZ-EL grid
n_zone=numel(vaz);

% special angles
[~,iaz_disp]=arrayfun(@(x) min(abs(Vaz-x)),az_disp);     % idx to displayable azim
[~,iel_disp]=arrayfun(@(x) min(abs(Vel-x)),el_disp);     % idx to displayable elev
naz_disp=length(iaz_disp);
nel_disp=length(iel_disp);
[~,iel_0]=min(abs(Vel));        % idx to ~zero elev angle (equator)


%% main
% preallocate
Jm_AB=cell(numel(g_inf),n_az,n_el);
std_Jm_AB=NaN(numel(g_inf),n_az,n_el);
n_det_events=NaN(numel(g_inf),n_az,n_el);

for mm=1:numel(g_inf)      % for measurement configuration
    tk=k_mm{mm};
    
    progressbar();
    for kk=1:n_zone
        [iaz,iel]=ind2sub([n_az,n_el],kk);
        
        taz=vaz(kk);
        tel=vel(kk);
        
        zone_k{1}=[taz,tel];       % az-el vector to A
        zone_k{2}=[taz+pi,-tel];   % B: sph-opposite to A
        
        % get num counts in A/B
        N=cell(1,2);        % spin-resolved numbers in A/B
        Ntot=cell(1,2);
        for ii=1:2
            [~,~,N{ii}]=cellfun(@(x) inCone(x,zone_k{ii}(1),zone_k{ii}(2),alpha),tk,'UniformOutput',false);
            N{ii}=cellfun(@(x) x,N{ii});
            
            Ntot{ii}=sum(N{ii},2);      % get total counts
        end
        
        % post-selection: single-atom detection at A/B
        b_1atom_k=cellfun(@(n) n==1, Ntot, 'UniformOutput',false);
        b_1atom=and(b_1atom_k{:});
        
        N_ps=cellfun(@(x) x(b_1atom,:),N,'UniformOutput',false);
        
        % angular momentum
        %   0.5*(N_up - N_down)
        Jm=cellfun(@(n) 0.5*(n(:,1)-n(:,2)),N_ps,'UniformOutput',false);
        
        % collective measurement
        tJm_AB=g_inf{mm}*Jm{1}+Jm{2};       % linear estimation
        tn_det_events=numel(tJm_AB);        % successfully post-selected events
        tstd_Jm_AB=std(tJm_AB);             % variation in coll meas of spin
        
        % store
        Jm_AB{mm,iaz,iel}=tJm_AB;
        n_det_events(mm,iaz,iel)=tn_det_events;
        std_Jm_AB(mm,iaz,iel)=tstd_Jm_AB;
        
        progressbar(kk/n_zone);
    end
end


%% EPR-steering
S_epr_lim=1/4;      % EPR-steering bound (spin-1/2)
S_epr=squeeze(std_Jm_AB(1,:,:).*std_Jm_AB(3,:,:));       % [Jx^B,Jz^B] EPR-S parameter

%% vis: EPR-steering
figname=sprintf('EPR');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

% inversion symmetry
[vaz_sym,vel_sym,S_epr_sym]=autofill_cent_symm(vaz,vel,S_epr);

p=plotFlatMap(rad2deg(vel_sym),rad2deg(vaz_sym),S_epr_sym,'eckert4','texturemap');
colormap('viridis');

% annotation
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
ax.LineWidth=ax_lwidth;

% titlestr=sprintf('%s%0.2f %s; %s%0.1g(%0.1g)','$\tau=$',tau(ii),'ms',...
%     '$\overline{\mathrm{SE}}=$',mu_seG,sig_seG);  % label expparam
% title(titlestr);

cbar=colorbar('SouthOutside');
cbar.TickLabelInterpreter='latex';
cbar.Label.Interpreter='latex';
cbar.Label.String='$\Delta_{\mathrm{inf}} S_{z}^{(B)} \Delta_{\mathrm{inf}} S_{x}^{(B)}$';
cbar.Label.FontSize=fontsize;


%% vis: histogram EPR-steering parameter
figname=sprintf('EPR_histogram');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

H=histogram(S_epr(:),'Normalization','pdf');
H.FaceColor=cviridis(4,:);
H.FaceAlpha=1;      % opaque
hold on;

% EPR-steering limit
ax=gca;
ax_ylim=ax.YLim;
ax_xlim=ax.XLim;
p_epr=patch([ax_xlim(1),S_epr_lim,S_epr_lim,ax_xlim(1)],...
    [ax_ylim(1),ax_ylim(1),ax_ylim(2),ax_ylim(2)],...
    c_gray,'EdgeColor','none');
uistack(p_epr,'bottom');
text(S_epr_lim,mean(ax_ylim),sprintf('EPR-steering'),...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'Rotation',90,'FontSize',fontsize-1);


% annotation
set(ax,'Layer','Top');
ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
xlabel('Inferred uncertainty $\Delta_{\mathrm{inf}} S_{z}^{(B)} \Delta_{\mathrm{inf}} S_{x}^{(B)}$');
ylabel('PDF');
% lgd=legend(H);    


%% vis: collective momenta
figname=sprintf('collective_momenta');
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

hold on;
H=[];
for ii=1:numel(g_inf)
    H(ii)=histogram(cat(1,Jm_AB{ii,:,:}),'Normalization','pdf',...
        'DisplayName',str_mm{ii});
end

% annotation
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;
xlabel('Collective spin measurement');
ylabel('PDF');
lgd=legend(H);    




%%
% %% vis
% %% observer
% h=figure;
%
% subplot(1,2,1);
% hold on;
% H(1)=histogram(n_A(:,1),'DisplayName','$\uparrow$');
% H(2)=histogram(n_A(:,2),'DisplayName','$\downarrow$');
%
% xlabel('$N_{m}$');
% ylabel('\# observations');
% lgd=legend(H);
% title('Alice');
%
% subplot(1,2,2);
% hold on;
% H(1)=histogram(n_B(:,1),'DisplayName','$\uparrow$');
% H(2)=histogram(n_B(:,2),'DisplayName','$\downarrow$');
%
% xlabel('$N_{m}$');
% ylabel('\# observations');
% lgd=legend(H);
% title('Bob');
%
% %% angular momentum
% h=figure;
% hold on;
%
% H(1)=histogram(Jn_A,'DisplayName','Alice');
% H(2)=histogram(Jn_B,'DisplayName','Bob');
%
% xlabel('Angular momentum $J_{i}$');
% ylabel('\# observations');
% lgd=legend(H);
%
% %% tot atoms
% h=figure;
% hold on;
%
% H(1)=histogram(ntot_A,'DisplayName','Alice');
% H(2)=histogram(ntot_B,'DisplayName','Bob');
%
% xlabel('$N_{\uparrow} + N_{\downarrow}$');
% ylabel('\# observations');
% lgd=legend(H);
