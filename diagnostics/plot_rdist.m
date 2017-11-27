function [nn,r]=plot_rdist(kk,Nbins)
% Plot normalised radial histogram
%
% [n,r]=plot_rdist(kk)
%
% kk:   nCountsx3 vector array of counts
% N:    number of bins to divide min-max (default 50)
%
% nn:   normalised radial count density
% r:    radial bin centers
%

if ~exist('N','var')
    Nbins=50;
end

knorm=vecnorm(kk,2,2);      % get vector norms

r_ed=linspace(min(knorm),max(knorm),Nbins+1)';
r=r_ed(1:end-1)+0.5*diff(r_ed);

Nr=nhist(knorm,{r_ed});      % do the 1D histogram in vector norms

%%% normalise distribution
% normalise radially dependent bin volume ~ 4*pi*r^2*dr
vol_r_bin=4*pi*(r.^2).*diff(r_ed);
nn=Nr./vol_r_bin;

% PDF: normalise by total and bin-size
nn=nn./(sum(nn)*diff(r_ed));


%%% plot
bar(r,nn,1,'FaceAlpha',0.5);

% annotate
box on;

end