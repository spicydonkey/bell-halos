% Generate ideal correlated source
% * mf=0,1 halos
%
% generates halo_k0 ideal source of momentum-spin correlated atoms
%

clear all; clc; close all;

% configs
n_shot=1e4;             % number of shots
n_corr_pairs=600*ones(n_shot,1);        % number of correlated pairs generated from source
det_qe=0.1*ones(1,2);       % detection efficiency

% construct ideal distinguishable halo experiment
halo_k0=cell(n_shot,2);
for ii=1:n_shot
    % build correlated source for this shot
    halo_temp=get_rand_usph(n_corr_pairs(ii));      % create user defined number of atoms on u-sphere
    
    halo_k0{ii,1}=halo_temp;        % mF=0 atoms
    halo_k0{ii,2}=-halo_temp;       % mF=1 atoms are just BB-pairs
    
    % random sampling for non-ideal detector efficiency
    for jj=1:2
        is_detected=(rand(n_corr_pairs(ii),1)<det_qe(jj));
        halo_k0{ii,jj}=halo_k0{ii,jj}(is_detected,:);
    end
end

%% cull shots with empty halo
% TODO - not doing it for now - will reduce correlation

%% load configs
run('config_idealsource.m');

%% evaluate correlations
run('eval_bell_corr.m');

