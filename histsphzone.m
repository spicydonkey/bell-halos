function [nSter,f_cents] = histsphzone(F,Az,El,f_edges)
% Solid-anglular histogram of values defined on a sphere
%
% NSTER = HISTSPHZONE(F, AZ, EL, F_EDGES)
%
%
% NOTE
% - ensure the zones are defined "uniformly" on a sphere:
%   - Az is linspace(0,pi,Nazim+1) (exclude last element)
%   - El is linspace(-pi/2,pi/2,Nelev)
% - TODO: azimuthal angles are unused by the algorithm - error checking in future
%

dSter=2*(pi^2)/numel(El)*cos(El);   % solid angle at each zone
df=mean(diff(f_edges));
f_cents=f_edges(1:end-1)+0.5*diff(f_edges);

nbin=length(f_edges)-1;
nSter=NaN(nbin-1,1);    % preallocate histogram

for ii=1:nbin
    this_bin_F=(F>f_edges(ii))&(F<=f_edges(ii+1));
    nSter(ii)=sum(dSter(this_bin_F));   % sum all differential solid angles at the corresponding zones
end
nSter=nSter/df;     % normalise by bin width

end