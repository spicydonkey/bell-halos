%%% load data
% add MACHINE_DEPENDENT\bell\out\20171128 to path
load('run2_K_20171127.mat');

%%% config
% zones
nAz=100;
nEl=50;
vaz=linspace(-pi,pi,nAz);
vaz=vaz(1:end-1);
vel=linspace(-pi/2,pi/2,nEl);
[az,el]=ndgrid(vaz,vel);

% sph-polar sampling
dpsi=0.3;
r_ed=linspace(0.8,1.2,100);

%% main
h_r0_map=figure();

K_new=cell(size(K));
for mm=1:2
    
    tKK=K(:,mm);         % sample data for prototyping
    tkk=vertcat(tKK{:});
    
    %%% sph-zone radial profiling
    % get normed r-profile
    [tnr,trc]=get_halo_sph_r_dist(tkk,az,el,dpsi,r_ed);
    % halo caps to NaN
    b_poles=abs(el)>asin(0.8);      % 2d boolean array to poles
    % nr(repmat(b_poles,[1,1,size(nr,3)]))=NaN;
    
    % get peak radius at each zone
    % TODO
    % * [ ] fit 1d-gaussian to get peak location
    r0=NaN(size(az));
    for ii=1:numel(r0)
        [I,J]=ind2sub(size(az),ii);
        [~,idx_max_nr]=max(tnr(I,J,:));
        r0(ii)=trc(idx_max_nr);
    end
    r0(b_poles)=NaN;        % halo caps to NaN
    
    % summary
    figure(h_r0_map); hold on;
    subplot(1,2,mm);
    plotFlatMapWrappedRad(az,el,r0,'eckert4');
    hcb=colorbar('Southoutside');
    hcb.Label.String='K_{pk}';
    drawnow;
    
    %%% transform
    tKK_new=tKK;
    for ii=1:numel(tKK)
        ts=zxy2sphpol(tKK{ii});
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
            if isnan(tr0)
                tr0=1;      % hack - won't transform counts with undefined radius to zone
                % TODO - could fit r-profile to the whole halo
            end
            
            % transform this vector
            ts_new(jj,3)=ts(jj,3)/tr0;
        end
        
        % update new KK in cartesian
        tKK_new{ii}=sphpol2zxy(ts_new);
    end
    K_new(:,mm)=tKK_new;    % store
end
% summary
h_Knew=figure();
plot_zxy(K_new);
axis equal; 
box on;
title('radially transformed halos');
drawnow; 