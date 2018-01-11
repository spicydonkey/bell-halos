function nn_halo=haloZoneCount(halo_k,az,el,sig,lim,histtype)
% Counts number of vectors around sph-polar zones
%   a wrapper for testing different zone-types
%
% NN_HALO = HALOZONECOUNT(HALO_K, AZ, EL, SIG, LIM, HISTTYPE)
%
% HALO_K: Nx3 array (shot) (TODO - not for latlon)
% 

% Count each zone
%
% TODO
% [ ] compare against the lat-lon counting
% [x] how to handle BAD regions: bad regions are post-processed by NaN-padding
% [x] ONLY accept SHOT data (array)
%   [x] gauss
%       [x] TEST
%   [x] simple
%       [ ] TEST
% [x] supply latlon zones
% [ ] document code


switch histtype
    case 'simple'
        % simple atom counting in bin: convolution with top-hat/rectangle filter
        % i.e. counts atoms "in" a defined solid-angle 
        % TODO - radial
        
		% get bin widths
        dpsi_max=sig(1);    % max relative angle from bin-vector
        
        % counting atoms in bins
        k_sph=zxy2sphpol(halo_k);   % cart to sph-polar
        nn_halo=zeros(size(az));  % preallocate bin counts
        for ii=1:numel(az)
            dpsi=diffAngleSph(k_sph(:,1),k_sph(:,2),az(ii),el(ii));
            nn_halo(ii)=sum(dpsi<dpsi_max);     % number of atoms in this zone
        end
        
    case 'gauss'
        warning('gauss mode: grid must be linearly spaced and complete.');
        [nAz,nEl]=size(az);
        nn_halo=wHaloDensity(halo_k,nAz,nEl,sig,lim);        % gaussian weighted counting
            
    otherwise
        error('unknown histtype. histtype should be either simple or gauss.');
end

end
