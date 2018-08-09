% analyse automated two-pulse Raman process
%   Ramsey interferometry for long delay
%   
%   approx optimal polarization to eliminate sigma+ component
%   global rotation - beam very large compared to spatial distribution of atoms
%   rotation on (1,1) halo source
%
%   
%
%   NOTE:
%       * significant instability in PUSH/NULLER-SO sequence: mf=1 position variability
%
% 20180409
% DKS
%


%% load raw data
run('config_v1');

%%%% load param log
param_log=load_logfile(configs.path.paramlog);
param_array = paramlog2array(param_log);

% get unique param-vecs and tag each shot with param-ID
[params,~,Ipar] = unique(param_array(:,2:end),'rows');

param_id=param_array(:,1);
nparam=size(params,1);      % number of unique param-set

% group shot-ids by exp-param
id_in_param=cell(1,nparam);
for ii=1:nparam
    id_in_param{ii}=param_id(Ipar==ii);
end

%% START debug - param search manual
%%% get searched param
% unique params
upar_dt=unique(params(:,1));

%% END debug - param search manual

%%%% load txy
[txy,fout]=load_txy(configs.load.path,configs.load.id,configs.load.window,configs.load.mincount,configs.load.maxcount,[],1,0,0);

zxy=txy2zxy(txy,vz);
shot_id=fout.id_ok;

% build boolean array selector for scanned param-set
b_paramset=cellfun(@(idx) ismember(shot_id,idx),...
    id_in_param,'UniformOutput',false);
b_paramset=horzcat(b_paramset{:});


% DEBUG
h_txy_raw=figure;
plot_zxy(txy,1e5);
axis equal;
xlabel('x');
ylabel('y');
zlabel('t');
view([0,0]);

h_zxy_raw=figure;
plot_zxy(zxy,1e5);
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);



%% distinguish mF and capture halo
n_shot=size(zxy,1);
n_mf=numel(configs.mf);

% preallocate
p_bec=cell(1,n_mf);     % marker BEC center positions
p_bec0=cell(1,n_mf);    % marker BEC in halo centered ZXY coord sys
p_halo=cell(1,n_mf);   % approx halo center (mid-point of marker BECs)
for ii=1:n_mf
    p_bec{ii}=NaN(n_shot,3,2);
    p_halo{ii}=NaN(n_shot,3);
end
zxy0=cell(n_shot,n_mf); % halo centerised zxy


for ii=1:n_shot
    tzxy=zxy{ii};
    for jj=1:n_mf
        tzxy_mf=boxcull(tzxy,configs.mf(jj).window);
        tp_bec0=configs.mf(jj).p_bec0;
        tr_bec0=configs.mf(jj).r_bec0;
        
        for kk=1:size(tp_bec0,1)
            tp_bec=capture_bec(tzxy_mf,tp_bec0(kk,:),tr_bec0,0);
            p_bec{jj}(ii,:,kk)=tp_bec;
        end
        tp_halo=mean(p_bec{jj}(ii,:,:),3);
        p_halo{jj}(ii,:)=tp_halo;
        zxy0{ii,jj}=tzxy_mf-tp_halo;
        p_bec0{jj}(ii,:,:)=p_bec{jj}(ii,:,:)-tp_halo;
    end
end

% DEBUG
figure(h_zxy_raw);
hold on;
plot_zxy(p_halo,[],50,'mkg');

h_zxy0=figure();
plot_zxy(zxy0,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);




%% filter data
zxy0_filt=zxy0;     % initialise filtered data

%%% BEC + thermal
r_thermal=configs.filt.r_ball;
for ii=1:n_shot
    for jj=1:n_mf
        for kk=1:2      % pair of marker BECs
            tp_bec0=p_bec0{jj}(ii,:,kk);
            % get atoms outside ball centered at BEC
            zxy0_filt{ii,jj}=cropBall(zxy0_filt{ii,jj},r_thermal,tp_bec0,false);
        end
    end
end

% DEBUG
h_zxy0_filt_1=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%%% radial 
% here we rely on "average radius" of halo determined by the mean marker
% BEC locations

% estimate halo radius from marker BECs
r_crop=configs.filt.r_crop;
r_halo=NaN(n_shot,n_mf);
for ii=1:n_mf
    r_halo(:,ii)=0.5*vnorm(p_bec0{ii}(:,:,2)-p_bec0{ii}(:,:,1));
end
r_halo_avg=mean(r_halo,1);

% filter
for ii=1:n_mf
    tr_lim=r_crop*r_halo_avg(ii);       % absolute norm limits
	zxy0_filt(:,ii)=cfilter_norm(zxy0_filt(:,ii),tr_lim(1),tr_lim(2));
end

% DEBUG
h_zxy0_filt_2=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%%% polar caps
% build box window for halo between caps
z_cap=configs.filt.z_cap;
window_z_filt=cell(1,n_mf);
for ii=1:n_mf
    window_z_filt{ii}={z_cap*r_halo_avg(ii)*[-1,1],[],[]};
end

% filter
for ii=1:n_mf
    tbox_lim=window_z_filt{ii};
    zxy0_filt(:,ii)=cellfun(@(x) boxcull(x,tbox_lim),zxy0_filt(:,ii),'UniformOutput',false);
end

% DEBUG
h_zxy0_filt_3=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%% k-space and distortion cancellation
%%% unit spherise
k_halo=zxy0_filt;       % initialise atoms in k-space (i.e. atoms lie on unit-sphere)

for ii=1:n_mf
    k_halo(:,ii)=cellfun(@(x) x/r_halo_avg(ii),k_halo(:,ii),'UniformOutput',false);
end

%%% ellipsoid fit to sphere
for ii=1:n_mf
    k_halo(:,ii)=map2usph(k_halo(:,ii));
end
    
% DEBUG
scatter_halo(k_halo);


%% filter post-processed data
%%% radial
r_crop2=configs.filt2.r_crop;
k_halo_filt=cfilter_norm(k_halo,r_crop2(1),r_crop2(2));

%%% z-cap
% build box window for halo between caps
z_cap2=configs.filt2.z_cap;
window_z_filt2=cell(1,n_mf);
for ii=1:n_mf
    window_z_filt2{ii}={z_cap2*[-1,1],[],[]};
end

% filter
for ii=1:n_mf
    tbox_lim=window_z_filt2{ii};
    k_halo_filt(:,ii)=cellfun(@(x) boxcull(x,tbox_lim),k_halo_filt(:,ii),'UniformOutput',false);
end


% DEBUG
scatter_halo(k_halo_filt);


%% categorise data by exp-params
k_par=cell(1,nparam);
for ii=1:nparam
    k_par{ii}=k_halo_filt(b_paramset(:,ii),:);      % get all halo data
    %from this param-set and store
end

% % DEBUG
% figure;
% for ii=1:nparam
%     subplot(1,nparam,ii);
%     plot_zxy(k_par{ii});
%     
%     axis equal;
%     xlabel('kx');
%     ylabel('ky');
%     zlabel('kz');
%     view([0,0]);
% end


%%% END OF PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ANALYSIS
dt=upar_dt;


%% Rabi oscillation
Nmf=cellfun(@(x) shotSize(x),k_par,'UniformOutput',false);

% statistics
Nmf_avg=cellfun(@(n) mean(n,1),Nmf,'UniformOutput',false);
Nmf_avg=vertcat(Nmf_avg{:});

Nmf_std=cellfun(@(n) std(n,0,1),Nmf,'UniformOutput',false);
Nmf_std=vertcat(Nmf_std{:});

Nmf_se=Nmf_std./sqrt(cellfun(@(x)size(x,1),Nmf))';


%% data vis
figure('Name','2-pulse_Rabi_osc');
hold on;

p1=ploterr(1e6*dt,Nmf_avg(:,1),[],Nmf_se(:,1),'ro');
p2=ploterr(1e6*dt,Nmf_avg(:,2),[],Nmf_se(:,2),'bo');

p1(1).DisplayName='1';
p2(1).DisplayName='0';

titlestr=sprintf('2-pulse Raman: pi/2 with delay');
title(titlestr);
xlabel('Pulse delay (us)');
ylabel('Number in ROI');
box on;
ax=gca;
xlim([0,ax.XLim(2)]);
ylim([0,ax.YLim(2)]);

lgd=legend([p1(1),p2(1)]);
lgd.Title.String='$m_F$';


%% non-dimensionalised time: Larmor period
% f_larmor=
T_larmor=0.65e-6;       % larmor precession period (s)

dtau=dt/T_larmor;

figure('Name','ramsey_decay');
hold on;

p1=ploterr(dtau,Nmf_avg(:,1),[],Nmf_se(:,1),'ro');
p2=ploterr(dtau,Nmf_avg(:,2),[],Nmf_se(:,2),'bo');

p1(1).DisplayName='1';
p2(1).DisplayName='0';

titlestr=sprintf('2-pulse Raman: pi/2 with delay');
title(titlestr);
xlabel('$\tau$');
ylabel('Number in ROI');
box on;
ax=gca;
xlim([0,ax.XLim(2)]);
ylim([0,ax.YLim(2)]);

lgd=legend([p1(1),p2(1)]);
lgd.Title.String='$m_F$';


%% Analysis of decoherence
Nmf_max=max(Nmf_avg,[],1);      % max number for normalisation

% normalise number s.t. max = 1
nmf_avg=Nmf_avg./Nmf_max;
nmf_se=Nmf_se./Nmf_max;

nmf_tot=sum(nmf_avg,2);


% vis
figure('Name','normalised_number');
hold on;

p1=ploterr(dtau,nmf_avg(:,1),[],nmf_se(:,1),'ro');
p2=ploterr(dtau,nmf_avg(:,2),[],nmf_se(:,2),'bo');

p1(1).DisplayName='1';
p2(1).DisplayName='0';

titlestr=sprintf('number normalised');
title(titlestr);
xlabel('$\tau$');
ylabel('n');
box on;
ax=gca;
xlim([0,ax.XLim(2)]);
ylim([0,ax.YLim(2)]);

lgd=legend([p1(1),p2(1)]);
lgd.Title.String='$m_F$';


%% population fraction
P=nmf_avg./nmf_tot;
P_err=(vnorm(nmf_se,2)./nmf_tot);      % simple estimate of error

% vis
figure('Name','population_fraction');
hold on;

p1=ploterr(dtau,P(:,1),[],P_err,'ro');
p2=ploterr(dtau,P(:,2),[],P_err,'bo');

p1(1).DisplayName='1';
p2(1).DisplayName='0';

titlestr=sprintf('Rabi oscillation decay');
title(titlestr);
xlabel('$\tau$');
ylabel('P');
box on;
ax=gca;
xlim([0,ax.XLim(2)]);
ylim([0,ax.YLim(2)]);

lgd=legend([p1(1),p2(1)]);
lgd.Title.String='$m_F$';


%% decay model
pp=P(:,1);


% v1
% modelfun = @(b,x) b(1)*exp(-b(2)*x).*cos(b(3)*(2*pi)*x) + b(4);
% beta0 = [1,1e-2,1,0.5];

% v2
modelfun = @(b,x) 0.5*exp(-b(1)*x).*cos(b(2)*(2*pi)*x) + 0.5;
beta0=[1e-2,1];

mdl = fitnlm(dtau,pp,modelfun,beta0);

%%% eval fit
beta_fit=mdl.Coefficients.Estimate;

tfit=linspace(0,1200,1e5);
pfit=modelfun(beta_fit,tfit);


% vis
figure('Name','population_fraction');
hold on;

p1=ploterr(dtau,P(:,1),[],P_err(:,1),'ro');
p1(1).DisplayName='data';

titlestr=sprintf('Rabi oscillation decay');
title(titlestr);
xlabel('$\tau$');
ylabel('$P_1$');
box on;
ax=gca;
xlim([0,ax.XLim(2)]);
ylim([0,ax.YLim(2)]);

p_mdl=plot(tfit,pfit,'-','Color',0.85*[1,1,1],'DisplayName','fit');
uistack(p_mdl,'bottom')

lgd=legend([p1(1),p_mdl]);
% lgd.Title.String='$m_F$';