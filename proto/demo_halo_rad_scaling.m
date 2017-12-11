%%% load data
load('run1_K_20171127.mat');

%%% config
% zones
nAz=50;
nEl=25;
vaz=linspace(-pi,pi,nAz);
vaz=vaz(1:end-1);
vel=linspace(-pi/2,pi/2,nEl);
[az,el]=ndgrid(vaz,vel);

% sph-polar sampling
dpsi=0.33;
r_ed=linspace(0.8,1.2,100);

%% main
kk=vertcat(K{1:1000,1});        % sample of mf=0 halo

% get normed r-profile
[nr,rc]=get_halo_sph_r_dist(kk,az,el,dpsi,r_ed);
% halo caps to NaN
b_poles=abs(el)>asin(0.8);
% nr(repmat(b_poles,[1,1,size(nr,3)]))=NaN;

% get peak radius at each zone
% TODO 
% * [ ] fit 1d-gaussian to get peak location
r0=NaN(size(az));
for ii=1:numel(r0)
    [I,J]=ind2sub(size(az),ii);
    [~,idx_max_nr]=max(nr(I,J,:));
    r0(ii)=rc(idx_max_nr);
end
r0(b_poles)=NaN;        % halo caps to NaN

% summary
figure();
plotFlatMapWrappedRad(az,el,r0,'eckert4');
title('peak radius around halo');

% %% transform