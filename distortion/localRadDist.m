function [r0,r_sig]=localRadDist(zxy,az,el,theta,verbose)
%Summary of local radial distribution.
%
%   zxy: cart-vect array
%   az: cone azim
%   el: cone elev
%   theta: cone half-angle
%
%   r0: mean radius
%   r_sig: rms radius
%
%
%   USAGE
%   [az,el]=sphgrid(100,50);
%   [r0,dr]=arrayfun(@(a,e)
%       localRadDist(double(vertcat(zxy{:})),a,e,0.3),az,el);
%

if ~exist('verbose','var')
    verbose=0;
end

% configs
min_n_cone=200;     % min vecs in cone reqd for dist analysis
n_per_bin=10;       % avg number of counts per bin
n_sm=5;             % moving avg filter span

% vectors in locally selected cone
[zxy_in]=inCone(zxy,az,el,theta);
% check for sufficient number of vectors
if size(zxy_in,1)<min_n_cone
    r0=NaN;
    r_sig=NaN;
    return
end

% get norms
r=vecnorm(zxy_in);

% radial profile
red=linspace(min(r),max(r),numel(r)/n_per_bin);
[nr,rcent]=radialprofile(r,red,verbose);
nr=smooth(nr,n_sm);      % smooth

% get approximate profile params
[amp0,Imax]=max(nr);
mu0=rcent(Imax);
sigma0=std(r);

% gaussian fit
peq={[],[],[],0};       % constrain to zero offset (c)
p0=[amp0,mu0,sigma0];   % initial param for optimiser
pfit=fit_gauss_1d(rcent,nr,p0,peq,verbose);

% get distribution statistics
r0=pfit(2,1);
r_sig=pfit(3,1);

end