function [P_rabi,z_az,z_el]=rabiAnalyse(halo_k,nAz,nEl,sig,lim)
% [P_RABI, Z_AZ, Z_EL] = RABIANALYSE(HALO_K, NAZ, NEL, SIG, LIM)
%
%

% define spherical momentum zones
az=linspace(-pi,pi,nAz);
az=az(1:end-1);             % unique angles only
el=linspace(-pi/2,pi/2,nEl);
[z_az,z_el]=ndgrid(az,el);      % grid to preserve dimensional ordering in array indexing


% TODO
% here simplify by collating all shots
halo_k_combined=cell(1,2);
for ii=1:2
    halo_k_combined{ii}=vertcat(halo_k{:,ii});
end

nn_halo=cellfun(@(halo) wHaloDensity(halo,nAz,nEl,sig,lim),halo_k_combined,'UniformOutput',false);
% nn_halo=wHaloDensity(zxy,nAz,nEl,sig,lim);
% tnn=cell(1,2);
% for ii=1:2
%     tnn{ii}=cat(3,nn_halo{:,ii});
% end
% nn_halo=tnn;
% clearvars tnn;

% TODO
% smooth nn_halo before processing P_rabi


% Rabi state population (population in mf=1)
P_rabi=nn_halo{2}./(nn_halo{1}+nn_halo{2});

end
