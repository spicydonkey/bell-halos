function [n,r,b]=plot_rdist(kk,Nbins)
% Plot normalised radial histogram
%
% [n,r] = plot_rdist(kk,Nbins)
%
% kk:       nCountsx3 vector array of counts
% Nbins:    number of bins to divide min-max (default 50)
%
% n:        normalised radial count density
% r:        radial bin centers
%

% parse inputs
if ~exist('N','var')
    Nbins=50;
end

%%% Histogram
% knorm=vecnorm(kk,2,2);      % get vector norms
knorm=myvecnorm(kk);      % get vector norms

% build bins for histogram 
r_ed=linspace(min(knorm),max(knorm),Nbins+1)';
% r_ed=linspace(min(knorm),max(knorm),Nbins+1);
r=r_ed(1:end-1)+0.5*diff(r_ed);

% do the 1D histogram in vector norms
Nr=nhist(knorm,{r_ed});
% Nr=nhist(knorm,{r_ed})';

%%% normalise distribution
% normalise radially dependent bin volume ~ 4*pi*r^2*dr
vol_r_bin=4*pi*(r.^2).*diff(r_ed);
n=Nr./vol_r_bin;

% PDF: normalise by total and bin-size
n=n./(sum(n)*diff(r_ed));


%%% plot
b=bar(r,n,1,'FaceAlpha',0.5);

% annotate
box on;
xlabel('r');
ylabel('PDF');

end