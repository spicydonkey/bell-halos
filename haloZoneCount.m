function [nn_halo,z_az,z_el]=haloZoneCount(halo_k,nAz,nEl,sig,lim,histtype)
% [NN_HALO, Z_AZ, Z_EL] = HALOZONECOUNT(HALO_K, NAZ, NEL, SIG, LIM, HISTTYPE)
%
%
% NN_HALO: 1x2 cell of counts in zone
%

% Count each zone
% TODO
% [] compare against the lat-lon counting
% [] how to handle BAD regions

switch histtype
    case 'gauss'
        % define spherical momentum zones
        az=linspace(-pi,pi,nAz);
        az=az(1:end-1);             % unique angles only
        el=linspace(-pi/2,pi/2,nEl);
        [z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing
        
        % TODO
        % do each shot and get error statistics
        
        % simplify by collating all shots
        halo_k_combined=cell(1,2);
        for ii=1:2
            halo_k_combined{ii}=vertcat(halo_k{:,ii});
        end
        
        % Gaussian
        nn_halo=cellfun(@(halo) wHaloDensity(halo,nAz,nEl,sig,lim),halo_k_combined,'UniformOutput',false);
        
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
