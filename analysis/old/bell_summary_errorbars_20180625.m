%% Summary of spin-corrleations in the scattered pairs of BEC collision:
%   * global rotation
%   * variable delay
%
%   with error bars!
%
% 20180507
% DKS


%% CONFIGS
idx_to_plot=[];        % set to [] to plot all


% linear approximation of Raman pulse duration and rotation angle
Tpi=10e-6;      % approx pi pulse duration


% vis
markershape={'o','^','s','*','x','d','o','^','s','*','x','d'};
markercolor=distinguishable_colors(15);


%% theory
% Bell triplet
theta0_th=linspace(0,2*pi,1e3);
E_th=-cos(2*theta0_th);

%%% Bell inequality limit
Emax_bell=1/sqrt(2);



%% Experimental Data - as code
% NOTE:
%   * error is evaluated from bootstrapping - see each analysis config
%   FORMAT
%
%%% EXP #. DETAILS
% idx=[];
% 
% T_raman{idx}=1e-5*[];      
% E_raw{idx}=[];
% E_raw_err{idx}=[];
% E0{idx}=[];
% E0_err{idx}=[];
% t_evo(idx)=[];
% dataname{idx}='';


%%% 1. low mode occupancy (exp5_low_n_theta_v1); 0.8 ms
idx=1;

T_raman{idx}=1e-05*[0    0.0625    0.1250    0.1875    0.3125    0.3750    0.4375];
E_raw{idx}=[-0.8727   -0.8474   -0.5811   -0.1800    0.4924    0.6892    0.8560];
E_raw_err{idx}=[0.0852    0.0614    0.1408    0.1369    0.1197    0.1274    0.0792];
E0{idx}=[-0.9792   -0.9189   -0.6240   -0.1938    0.5279    0.7312    0.9039];
E0_err{idx}=[0.0817    0.0956    0.1247    0.1519    0.1519    0.1266    0.0782];
t_evo(idx)=0.8e-3;       % delay SRC-ROT in seconds
dataname{idx}='0.8 ms (exp5.1)';


%%% 2. violation
idx=2;

T_raman{idx}=5e-6;      % pi/2
E_raw{idx}=[0.85];
E_raw_err{idx}=[0.05];
E0{idx}=[0.9395];
E0_err{idx}=[0.0725];
t_evo(idx)=0.8e-3;
dataname{idx}='0.8 ms (viol)';


%%% 3. ideal_global_rotation_v4: 0.8 ms 
idx=3;

T_raman{idx}=1e-5*[0.2500    0.5000];      
E_raw{idx}=[0.0254    0.7975];
E_raw_err{idx}=[0.0619    0.0412];
E0{idx}=[0.0310    0.9512];
E0_err{idx}=[0.0833    0.0420];
t_evo(idx)=0.8e-3;
dataname{idx}='0.8 ms (id-grot-v4)';


%%% 4. low mode occupancy (exp5_low_n_theta_v2); 0.8ms
idx=4;

T_raman{idx}=1e-5*[0.6250    0.7500    0.8750    1.0000];      
E_raw{idx}=[0.3412   -0.3363   -0.7481   -0.8261];
E_raw_err{idx}=[0.2053    0.1424    0.1579    0.0543];
E0{idx}=[];
E0_err{idx}=[];
t_evo(idx)=0.8e-3;
dataname{idx}='0.8 ms (exp5.2)';


%%% 5. Y-axis med-mocc; 0.8 ms
idx=5;

T_raman{idx}=1e-5*[0.0833    0.1667    0.2500    0.3333    0.4167    0.5000];      
E_raw{idx}=[-0.3672   -0.1348    0.1198    0.4447    0.5844    0.5854];
E_raw_err{idx}=[0.0843    0.0840    0.0609    0.0932    0.0759    0.0715];
E0{idx}=[];
E0_err{idx}=[];
t_evo(idx)=0.8e-3;
dataname{idx}='0.8 ms (y-v1-2)';


%%% 6. low mode occupancy (exp5_low_n_theta_v3); 0.8ms
idx=6;

T_raman{idx}=1e-5*[0.5625    0.6875    0.8125    0.9375];      
E_raw{idx}=[0.6829   -0.0093   -0.3331   -0.6548];
E_raw_err{idx}=[0.0667    0.1248    0.1055    0.0614];
E0{idx}=[];
E0_err{idx}=[];
t_evo(idx)=0.8e-3;
dataname{idx}='0.8 ms (exp5.3)';


%%% 7. t-evolution: 0.95ms
idx=7;

T_raman{idx}=1e-5*[0.5];      
E_raw{idx}=[0.6131];
E_raw_err{idx}=[0.0516];
E0{idx}=[1.0097];
E0_err{idx}=[0.1025];
t_evo(idx)=0.95e-3;
dataname{idx}='0.95 ms';


%%% 8. t-evolution: 1.1ms
idx=8;

T_raman{idx}=1e-5*[0.5];      
E_raw{idx}=[0.4756];
E_raw_err{idx}=[0.0562];
E0{idx}=[0.7614];
E0_err{idx}=[0.1008];
t_evo(idx)=1.1e-3;
dataname{idx}='1.1 ms';


%%% 9. t-evolution: 1.25ms
idx=9;

T_raman{idx}=1e-5*[0.5];      
E_raw{idx}=[0.3770];
E_raw_err{idx}=[0.0418];
E0{idx}=[0.6523];
E0_err{idx}=[0.0884];
t_evo(idx)=1.25e-3;
dataname{idx}='1.25 ms';


%%% 10. t-evolution: 1.7ms
idx=10;

T_raman{idx}=1e-5*[0.5];      
E_raw{idx}=[0.1735];
E_raw_err{idx}=[0.0600];
E0{idx}=[0.2648];
E0_err{idx}=[0.0987];
t_evo(idx)=1.7e-3;
dataname{idx}='1.7 ms';


%%% 11. t-evolution: 1.4ms: ideal_global_rotation_v3 (pi/2)
idx=11;

T_raman{idx}=1e-5*[0.5];      
E_raw{idx}=[0.5237];
E_raw_err{idx}=[0.0582];
E0{idx}=[0.5731];
E0_err{idx}=[0.0616];
t_evo(idx)=1.4e-3;
dataname{idx}='1.4 ms';



%% Calculation
% calculate rotation angle 
theta0=cellfun(@(t) (t/Tpi)*pi,T_raman,'UniformOutput',false);



%% Data visualisation
if isempty(idx_to_plot)
    idx_to_plot=1:length(T_raman);
end


%%% 1. RAW DATA
h_raw=figure('Name','corr_raw');
ax=gca;
hold on;

hplot=[];   % array to store plot handles

% Experimental data
for ii=1:numel(idx_to_plot)
    tii=idx_to_plot(ii);

    hh=ploterr(theta0{tii}/pi,E_raw{tii},...
        [],E_raw_err{tii},...
        markershape{ii},'hhxy',0);
    
    %%% annotation
    % NOTE: hh(1): marker; hh(2): Y-bar; hh(3): X-bar;
    set(hh(1),'Marker',markershape{ii},'MarkerSize',7,...
        'Color',markercolor(ii,:),'LineWidth',1.2,...
        'MarkerFaceColor','w',...               % FILLS with white
        'DisplayName',dataname{tii});
    set(hh(2),'Color',markercolor(ii,:),'LineWidth',1.2,...
        'DisplayName','');  % Y-err
%     set(hh(3));       % X-err

    hplot=[hplot,hh(1)];
    
end

% theory: Bell triplet 
p_theory=plot(theta0_th/pi,E_th,'k--','LineWidth',1.5,...
    'DisplayName','$\vert\Psi^+\rangle$ Bell triplet');

% Bell inequality limit
p_bell_lim=line(ax.XLim,Emax_bell*ones(1,2),...
    'Color','k','LineStyle','-','LineWidth',1.5,...
    'DisplayName','Bell');


%%% annotate
% ax=gca;
ax.FontSize=12;

max_theta0_plot=ceil(max(cellfun(@(th) max(th), theta0(idx_to_plot)))/pi);
xlim([0,max_theta0_plot]);
% xlim([0,2]);
ylim([-1,1]);

h_raw.Name='raw_corr';

title('Raw correlation');

xlabel('$\theta_0/\pi$');
ylabel('$E$');

% lgd=legend('$\vert\Psi^+\rangle$ Bell triplet','Bell',...
%     'Location','NorthEastOutside');
lgd=legend([hplot, p_theory, p_bell_lim],...
    'Location','NorthEast');
lgd.FontSize=10;

box on;


%%% 2. Normalised
h_norm=figure;
ax=gca;
hold on;

hplot=[];   % array to store plot handles

% Experimental data
for ii=1:numel(idx_to_plot)
    tii=idx_to_plot(ii);

    % TODO: fill out E0 data
    if isempty(E0{tii})
        continue;
    end
    
    hh=ploterr(theta0{tii}/pi,E0{tii},...
        [],E0_err{tii},...
        markershape{ii},'hhxy',0);
    
    %%% annotation
    % NOTE: hh(1): marker; hh(2): Y-bar; hh(3): X-bar;
    set(hh(1),'Marker',markershape{ii},'MarkerSize',7,...
        'Color',markercolor(ii,:),'LineWidth',1.2,...
        'MarkerFaceColor','w',...               % FILLS with white
        'DisplayName',dataname{tii});
    set(hh(2),'Color',markercolor(ii,:),'LineWidth',1.2,...
        'DisplayName','');  % Y-err
%     set(hh(3));       % X-err

    hplot=[hplot,hh(1)];
    
end

% theory: Bell triplet 
p_theory=plot(theta0_th/pi,E_th,'k--','LineWidth',1.5,...
    'DisplayName','$\vert\Psi^+\rangle$ Bell triplet');

% Bell inequality limit
p_bell_lim=line(ax.XLim,Emax_bell*ones(1,2),...
    'Color','k','LineStyle','-','LineWidth',1.5,...
    'DisplayName','Bell');


%%% annotate
% ax=gca;
ax.FontSize=12;

max_theta0_plot=ceil(max(cellfun(@(th) max(th), theta0(idx_to_plot)))/pi);
xlim([0,max_theta0_plot]);
% xlim([0,2]);
ylim([-1,1]);

h_raw.Name='corr_normalised';

title('Normalised correlation');

xlabel('$\theta_0/\pi$');
ylabel('$\bar{E}$');

% lgd=legend('$\vert\Psi^+\rangle$ Bell triplet','Bell',...
%     'Location','NorthEastOutside');
lgd=legend([hplot, p_theory, p_bell_lim],...
    'Location','NorthEast');
lgd.FontSize=10;

box on;




%% Time evolution of correlation @ pi/2 rotation
% get all data for pi/2 rotation
Tpi=10e-6;
Tpi2=Tpi/2;

is_pi2=cell(numel(idx_to_plot),1);

E0_pi2=NaN(numel(idx_to_plot),1);
E0_pi2_err=NaN(numel(idx_to_plot),1);

for ii=1:numel(idx_to_plot)
    tii=idx_to_plot(ii);
    
    is_pi2{tii}=(T_raman{tii}==Tpi2);
    
    if sum(is_pi2{tii})>0 && ~isempty(E0{tii})
        E0_pi2(tii)=E0{tii}(is_pi2{tii});
        E0_pi2_err(tii)=E0_err{tii}(is_pi2{tii});
    else
        E0_pi2(tii)=NaN;
        E0_pi2_err(tii)=NaN;
    end
    
end


%%% data vis
h = figure('Name','time_evolution_E_pi/2');

hh=ploterr(1e3*t_evo(idx_to_plot)',E0_pi2(idx_to_plot),...
    [],E0_pi2_err(idx_to_plot),...
    'o','hhxy',0);

% annotation
% NOTE: hh(1): marker; hh(2): Y-bar; hh(3): X-bar;
set(hh(1),'Marker','o','MarkerSize',7,...
    'Color','k','LineWidth',1.2,...
    'MarkerFaceColor','w');
set(hh(2),'Color','k','LineWidth',1.2,...
    'DisplayName','');  % Y-err
%     set(hh(3));       % X-err

title('Time-evolution of normalised spin-correlation under $\pi/2$ rotation');

xlabel('$t$ [ms]');
ylabel('$\bar{E}(\pi/2,\pi/2)$');