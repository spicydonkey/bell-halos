%% Pair-source g2 revisited


%% analysis - revisited
% to evaluate like Bell:
%   g2 for all data and configs: g2, dk, g2mdl
%   bootstrapping statistics:    g2_bootstrap, g2_mean, g2_sdev

% load dataset from pairsource

clearvars nparam k_par;

nparam=1;
k_par{1}=halo_k0;
par_T=NaN;


%% RUN analysis of exp5
% manually run the analysis part in main.m of exp5

% summary
g2_amp=g2mdl{1}{3}.Coefficients.Estimate(1);
g2_amp_sdev=g2_sdev{1}{3}(idx_dk0,idx_dk0,idx_dk0);


%% save revised results
% there's a lot of unnecessary outputs but we'll keep them

