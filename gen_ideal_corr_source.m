% Generate ideal correlated source
% * mf=0,1 halos
%
% generates halo_k0 ideal source of momentum-spin correlated atoms
%

clear all; clc; close all;

% configs
n_shot=1e4;             % number of shots
n_corr_pairs=100*ones(n_shot,1);        % number of correlated pairs generated from source
det_qe=0.1*ones(1,2);       % detection efficiency

% define local operation
fun_localoper=@(th,phi) sin(phi);
% fun_localoper=@(th,phi) 0;              % no local operation

% construct ideal distinguishable halo experiment
halo_k0=cell(n_shot,2);
for ii=1:n_shot
    %% build correlated source for this shot
    halo_temp=get_rand_usph(n_corr_pairs(ii));      % create user defined number of atoms on u-sphere
    
    halo_k0{ii,1}=halo_temp;        % mF=0 atoms
    halo_k0{ii,2}=-halo_temp;       % mF=1 atoms are just BB-pairs (INDEX BB MATCHED!)
    
    %% local operation
    theta_this=cell(2,1);
    for jj=1:2
        th_phi_this=zxy2pol_mat(halo_k0{ii,jj});
        th_phi_this=th_phi_this(:,[2,3]);   % grab azim,elev angles
        
        % evaluate rotation angle at each atom
        theta_this{jj}=fun_localoper(th_phi_this(:,1),th_phi_this(:,2));
    end
    dtheta_this=theta_this{2}-theta_this{1};    % evaluate relative rot angle for this shot
    Pr_flip=sin(dtheta_this).^2;        % probability to spin-flip
    
    % state rotation - quantum entanglement-like
    %   do a 50/50 coin toss on which state to apply DTHETA rotation
    collapse_sel=(rand(n_corr_pairs(ii),1)<0.5);	% collapse one of the pairs - mF
    halo_k0_temp=cell(1,2);     % temporary variable to hold post-mixing counts
    for jj=0:1
        % jj is mF; (jj+1) is index; (~jj)+1 gives index to other!
        bool_rotate=(collapse_sel==~jj);      % atoms to rotate are the ones that aren't selected as reference
        
        % append reference collapsed atoms
        halo_k0_temp{jj+1}=vertcat(halo_k0_temp{jj+1},halo_k0{ii,jj+1}(~bool_rotate,:));
        
        % rotation - monte-carlo
        bool_flip=(rand(n_corr_pairs(ii),1)<Pr_flip);       % atoms to flip (need to filter with ones to rotate)
        % rotated + unflipped
        halo_k0_temp{jj+1}=vertcat(halo_k0_temp{jj+1},halo_k0{ii,jj+1}(bool_rotate&(~bool_flip),:));
        % rotated + flipped
        halo_k0_temp{(~jj)+1}=vertcat(halo_k0_temp{(~jj)+1},halo_k0{ii,jj+1}(bool_rotate&bool_flip,:));
    end
    
    % update halo_k0 to spin-rotated: now of unequal length from spin rotation
    for jj=1:2
        halo_k0{ii,jj}=halo_k0_temp{jj};
    end
    
    %% random sampling for non-ideal detector efficiency
    for jj=1:2
%         is_detected=(rand(n_corr_pairs(ii),1)<det_qe(jj));
        is_detected=(rand(size(halo_k0{ii,jj},1),1)<det_qe(jj));
        halo_k0{ii,jj}=halo_k0{ii,jj}(is_detected,:);
    end
end
% leave empty halos as is

%% load configs
run('config_idealsource.m');

%% evaluate correlations
run('eval_bell_corr.m');

