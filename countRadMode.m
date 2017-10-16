function ww_rad = countRadMode(vs,rr,mode_thphi,sig,lim)
%
% WW_RAD = COUNTRADMODE(VS,RR,MODE_THPHI,SIG,LIM)
%
% WW_RAD: Mx1 weighted sum of counts in the radially resolved momentum mode
%
% VS: N-by-3 array of spherical vector of counts
% RR: Mx1 array of mode norms (modes colinear)
% MODE_THPHI: 1x2 array (lon,lat) defining mode direction
% SIG: 1x2 array of sigmas for gaussian weights (angular, norm)
% LIM: 1x2 array to delimit outliers in weighted counting defined in
% standard score (i.e. number of sigmas)
%

ww_rad=zeros(size(rr));
for ii=1:numel(rr)
    ww_rad(ii)=weightedCountSph(vs,[mode_thphi,rr(ii)],sig,lim);
end

end