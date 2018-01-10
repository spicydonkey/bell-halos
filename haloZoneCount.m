function [nn_halo,z_az,z_el]=haloZoneCount(halo_k,nAz,nEl,sig,lim,histtype)
% Counts number of vectors around sph-polar zones
%   a wrapper for testing different zone-types
%
% [NN_HALO, Z_AZ, Z_EL] = HALOZONECOUNT(HALO_K, NAZ, NEL, SIG, LIM, HISTTYPE)
%
% HALO_K: Nx3 array (shot) (TODO - not for latlon)
% 

% Count each zone
% TODO
% [] compare against the lat-lon counting
% [x] how to handle BAD regions: bad regions are post-processed by NaN-padding
% [] accept SHOT data (array)
%   [x] gauss
%       [x] TEST
%   [x] simple
%       [ ] TEST
%   [] latlon
% [] document code


switch histtype
    case 'simple'
        % simple atom counting in bin: convolution with top-hat/rectangle filter
        % i.e. counts atoms "in" a defined solid-angle 
        % TODO - radial
        
        dpsi_max=sig(1);
%         dr_max=sig(2);    % TODO - something like this
        
        % define spherical momentum zones
        az=linspace(-pi,pi,nAz);
        az=az(1:end-1);             % unique angles only
        el=linspace(-pi/2,pi/2,nEl);
        [z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing
        
        % counting atoms in bins
        k_sph=zxy2sphpol(halo_k);   % cart to sph-polar
        nn_halo=zeros(size(z_az));  % preallocate bin counts
        for ii=1:numel(z_az)
            dpsi=diffAngleSph(k_sph(:,1),k_sph(:,2),z_az(ii),z_el(ii));
            nn_halo(ii)=sum(dpsi<dpsi_max);     % number of atoms in this zone
        end
        
    case 'gauss'
        % define spherical momentum zones
        az=linspace(-pi,pi,nAz);
        az=az(1:end-1);             % unique angles only
        el=linspace(-pi/2,pi/2,nEl);
        [z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing

        nn_halo=wHaloDensity(halo_k,nAz,nEl,sig,lim);        % gaussian weighted counting
        
    case 'latlon'
        % lat-lon
        az=linspace(-pi,pi,nAz);
        el=linspace(-pi/2,pi/2,nEl+1);
        [z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing
        nn_halo=cell(1,2);
        for ii=1:2
            nn_halo{ii}=halo_zone_density(halo_k(:,ii),az,el,NaN,1);
            nn_halo{ii}=cat(3,nn_halo{ii}{:});  % nazim x nelev x nshot array
        end
        nn_halo=cellfun(@(x)mean(x,3)',nn_halo,'UniformOutput',false);
        
    otherwise
        error('unknown histtype. histtype should be either gauss or latlon.');
end

end
