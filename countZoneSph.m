function [nn_thphi,az,el] = countZoneSph(vs,nlatlon,sigPsi,limPsi)
%
% NN_THPHI = COUNTZONESPH(VS,NLATLON,SIG,LIM)
%
% NN_THPHI: normalised effective number in lat-lon zone
% AZ: azimuthal angles to zones
% EL: elevation angles to zones
%
% VS: N-by-3 array of spherical vector of counts
% NLATLON: 1x2 array of sphere subdivisions in long/lat-titudinal dimensions
% SIG: 1x2 array of sigmas for gaussian weights (angular, norm)
% LIM: 1x2 array to delimit outliers in weighted counting defined in
% standard score (i.e. number of sigmas)
%

azimUnique=linspace(-pi,pi,nlatlon(1)+1);
azimUnique=azimUnique(1:end-1);
elevUnique=linspace(-pi/2,pi/2,nlatlon(2));
[az,el]=ndgrid(azimUnique,elevUnique);
nZones=numel(az);

nCountsTot=size(vs,1);
weightedCounts=zeros(size(az));
for ii=1:nZones
    this_zone=[az(ii),el(ii),1];
    this_ww=weightedCountSph(vs,this_zone,[sigPsi,Inf],[limPsi,1]);
    weightedCounts(ii)=this_ww;
end

normFactor=nCountsTot*normSphCount(weightedCounts,az,el);
nn_thphi=weightedCounts*normFactor;

end