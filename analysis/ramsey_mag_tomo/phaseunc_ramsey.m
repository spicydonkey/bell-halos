% Ramsey interferometry
% phase uncertainty analysis
%
% 2019-07-04
% DKS
%

%% load data
load('C:\Users\HE BEC\Documents\MATLAB\bell-halos\analysis\ramsey_mag_tomo\out\dataforpaper\out_20190221_094756.mat');


%% analysis
% estimate phase from polarisation measurements
phi_est = cellfun(@(p) acos(p), P_k, 'uni', 0);

% statistics
phi_est_avg = cellfun(@(ph) mean(ph,3,'omitnan'), phi_est, 'uni', 0);
phi_est_avg = permute(cat(3,phi_est_avg{:}),[3,1,2]);   % form array as [tau X az X el]

phi_est_std = cellfun(@(ph) std(ph,0,3,'omitnan'), phi_est, 'uni', 0);
phi_est_std = permute(cat(3,phi_est_std{:}),[3,1,2]);


%%% theory: phase uncertainty 
eta_qe = 0.1;       % detector efficiency
N = 68;            % num atoms per phase estimation (estimate of actually scattered)
Dphi_theory = @(phi,N,eta) sqrt((1/N) * (1 + (1./sin(phi)).^2 * (1/eta - 1)));


%% plot
% config
iel_disp = iel_0;
iaz_disp = round(linspace(1,length(az),5));
iaz_disp = iaz_disp(1:end-1);

col_type = {'b','r','k'};
col_scatter = 0.3*ones(1,3);        
mark_type = {'o','^','s'};
lwidth = 1.5;
mark_size_small = 3;

% plot
n_disp = length(iaz_disp);

H = figure;

for ii=1:n_disp
    subplot(n_disp,1,ii)
    hold on;

    %%% theory
    tp_theory = ploterr(1e6*tau, phi_est_avg(:,iaz_disp(ii),iel_disp), [], Dphi_theory(phi_est_avg(:,iaz_disp(ii),iel_disp),N,eta_qe),'o','hhxy',1);
    set(tp_theory(1),'Visible','off');
    set(tp_theory(2),'Color',col_type{2},'DisplayName','theory');
    
    %%% exp
    % stat: avg and sdev
    tp_exp = ploterr(1e6*tau, phi_est_avg(:,iaz_disp(ii),iel_disp), [], phi_est_std(:,iaz_disp(ii),iel_disp),'o','hhxy',0.5);
    set(tp_exp(1),'Color',col_type{1},'LineWidth',lwidth,'Marker',mark_type{1},'MarkerFaceColor','w','DisplayName','mean, sdev');
    set(tp_exp(2),'Color',col_type{1},'LineWidth',lwidth);
    
    % scatter
    for jj=1:length(tau)
        ttau = tau(jj);
        tphi = squeeze(phi_est{jj}(iaz_disp(ii),iel_disp,:));
        
        tp_scatter = plot(1e6*ttau,tphi,'o','MarkerSize',mark_size_small,'Color',col_scatter,'MarkerFaceColor',col_scatter,'DisplayName','data');
    end
    
    
    % annotation
    xlim([1.8,3.6]);
    ylim([0,pi]);
    box on;
    
    xlabel('$\tau$ ($\mu$s)');
    ylabel('$\phi_{\textrm{est}}$ (rad)');
    
    titlestr = sprintf('$(\\theta, \\phi)$ = (%0.2f, %0.2f)',az(iaz_disp(ii)), el(iel_disp));
    title(titlestr);
    
    % legend
    legend([tp_scatter(1), tp_exp(1),tp_theory(2)],'Location','southwest');
end

set(H,'Renderer','painters');