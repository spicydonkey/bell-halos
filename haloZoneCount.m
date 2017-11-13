function [nn_halo,z_az,z_el]=haloZoneCount(halo_k,nAz,nEl,sig,lim,histtype)
% Counts number of vectors around sph-polar zones
%   a wrapper for testing different zone-types
%
% [NN_HALO, Z_AZ, Z_EL] = HALOZONECOUNT(HALO_K, NAZ, NEL, SIG, LIM, HISTTYPE)
%
%

% Count each zone
% TODO
% [] compare against the lat-lon counting
% [x] how to handle BAD regions: bad regions are post-processed by NaN-padding
% [] accept SHOT data (array)
%   [x] gauss
%       [] TEST
%   [] latlon
% [] document code


switch histtype
    case 'gauss'
        % define spherical momentum zones
        az=linspace(-pi,pi,nAz);
        az=az(1:end-1);             % unique angles only
        el=linspace(-pi/2,pi/2,nEl);
        [z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing
        
        % % simplify by collating all shots
        % halo_k_combined=cell(1,2);
        % for ii=1:2
        %     halo_k_combined{ii}=vertcat(halo_k{:,ii});
        % end
        % 
        % % Gaussian
        % nn_halo=cellfun(@(halo) wHaloDensity(halo,nAz,nEl,sig,lim),halo_k_combined,'UniformOutput',false);
        % 

        nn_halo=wHaloDensity(halo_k,nAz,nEl,sig,lim);

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
