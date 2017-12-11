function [nr,r_ct]=get_halo_sph_r_dist(k_zxy,th,phi,dpsi,r_ed)
% Get halo radial distribution around different spherical zones
%
% [nr,r_ct] = get_halo_sph_r_dist(k_zxy,th,phi,dpsi,r_ed)
%
% TODO
% [ ] write documentation

%% config
% radial density profiling
r_ct=r_ed(1:end-1)+0.5*diff(r_ed);
nBins_r=numel(r_ct);
nSmooth_r=round(nBins_r/10);    % simple moving average smoothing


%% Radial profiling around sph-zones
% bin-size of sph-section
vol_bin_r=4*pi*(r_ct.^2).*diff(r_ed);

nr=NaN([size(th),nBins_r]);     % preallocate radial density profile array

ss=zxy2sphpol(k_zxy);       % counts in sph-polar coord

nZ=numel(th);       % number of partitioned sph zones
sizth=size(th);
for ii=1:nZ
    [Ith,Iphi]=ind2sub(sizth,ii);
    
    % get this zone's azim-elev coords
    tth=th(ii);
    tphi=phi(ii);
    
    % get bool array to counts in this zone
    tb_in=inZone(ss(:,1),ss(:,2),tth,tphi,dpsi);
    
    % vec-norms of counts in zone
    tr_in=ss(tb_in,3);
    tNr=nhist(tr_in,{r_ed})';   % histogram
    % transposed to re-shape as row-vector
    
    % normalise to 1D-radial PDF
    tnr=tNr./vol_bin_r;
    tnr=tnr./(sum(tnr)*diff(r_ed));
    
    % simple moving average
    tnr=smooth(tnr,nSmooth_r);
    
    % normalise peak amplitude to unity
    tnr=tnr./max(tnr);
    
    nr(Ith,Iphi,:)=tnr;
end

end