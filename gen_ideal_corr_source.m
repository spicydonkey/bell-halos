% Generate ideal correlated source
% * mf=0,1 halos
%
% generates halo_k0 ideal source of momentum-spin correlated atoms
%

% configs
n_shot=1e3;             % number of shots
n_corr_pairs=4*ones(n_shot,1);        % number of correlated pairs generated from source
det_qe=[0.1,0.1];       % detection efficiency

% construct ideal distinguishable halo experiment
halo_k0=cell(n_shot,2);
for ii=1:n_shot
    % build correlated source for this shot
    halo_temp=get_rand_usph(n_corr_pairs(ii));      % create user defined number of atoms on u-sphere
    
    for jj=1:2
        halo_k0{ii,1}=halo_temp;        % mF=0 atoms
        halo_k0{ii,2}=-halo_temp;       % mF=1 atoms are just BB-pairs
    end
end