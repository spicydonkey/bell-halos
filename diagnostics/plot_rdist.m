function [n_r,r_cent,b]=plot_rdist(kk,Nbins)
% Plot normalised radial histogram
%
% [n,r] = plot_rdist(kk,Nbins)
%
% kk:       nCountsx3 vector array of counts
% Nbins:    number of bins to divide min-max (default 50)
%
% n_r:        normalised radial count density
% r_cent:        radial bin centers
% b: handle to bar plot
%

% parse inputs
if ~exist('Nbins','var')
    Nbins=50;
end

%%% Radial PDF
knorm=vnorm(kk);      % get vector norms

% build bins for histogram 
r_ed=linspace(min(knorm),max(knorm),Nbins+1);

[n_r,r_cent,b] = radialprofile(knorm,r_ed,1);

end