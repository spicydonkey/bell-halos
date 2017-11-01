function [P_rabi,z_az,z_el]=rabiAnalyse(halo_k,nAz,nEl,sig,lim,histtype)
% [P_RABI, Z_AZ, Z_EL] = RABIANALYSE(HALO_K, NAZ, NEL, SIG, LIM, HISTTYPE)
%
%

% Count each zone
% TODO
% [] compare against the lat-lon counting
% [] ignore BAD regions

switch histtype
    case 'gauss'
        % define spherical momentum zones
        az=linspace(-pi,pi,nAz);
        az=az(1:end-1);             % unique angles only
        el=linspace(-pi/2,pi/2,nEl);
        [z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing
        
        % TODO
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

% TODO
% smooth nn_halo before processing P_rabi

% Rabi state population (population in mf=1)
P_rabi=nn_halo{2}./(nn_halo{1}+nn_halo{2});

end
