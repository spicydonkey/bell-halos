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
mm=2;

KK=K(1:1000,mm);         % sample data for prototyping
kk=vertcat(KK{:});

%%% sph-zone radial profiling
% get normed r-profile
[nr,rc]=get_halo_sph_r_dist(kk,az,el,dpsi,r_ed);
% halo caps to NaN
b_poles=abs(el)>asin(0.8);      % 2d boolean array to poles
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

%%% transform
KK_new=KK;
for ii=1:numel(KK)
    ts=zxy2sphpol(KK{ii});
    ts_new=ts;
    
    % th-phi dependent radial scaling
    for jj=1:size(ts,1)
        % get sph-pol zone for this count
        %   simple nearest point
        tth=ts(jj,1);       % this count azim angle
        tphi=ts(jj,2);      % this count elev angle
        [~,idx_th]=min(abs(vaz-tth));
        [~,idx_phi]=min(abs(vel-tphi));
        
        % get normalising radius
        tr0=r0(idx_th,idx_phi);     
        
        % transform this vector
        ts_new(jj,3)=ts(jj,3)/tr0;
    end
    
    % update new KK in cartesian 
    KK_new{ii}=sphpol2zxy(ts_new);
end

% summary
figure();
plot_zxy(KK_new,1e4,10);
axis equal;
title('sph-zone radially transformed halo');
