%%% load data
load('run1_K_20171127.mat');

%%% config
% sph-polar sampling points
nAz=100;
nEl=50;
dpsi=0.2;   % angular diff to neighbour

% radial density profile
nrBins=100;
r_ed=linspace(0.8,1.2,nrBins+1);
r_c=r_ed(1:end-1)+0.5*diff(r_ed);
nsmooth=9;     % simple moving average smoothing


%% diagnostics
% radial dist
figure();
for ii=1:2
    hold on;
    plot_rdist(vertcat(K{:,ii}));
end
drawnow

% sph dist
figure();
for ii=1:2
    subplot(1,2,ii);
    plot_sphdist(vertcat(K{:,ii}));
end
drawnow

%% get radial distribution at each sph-polar zones
% build sph-polar grid
az_vec=linspace(-pi,pi,nAz+1);
az_vec=az_vec(1:end-1);     % unique points
el_vec=linspace(-pi/2,pi/2,nEl);
[az,el]=ndgrid(az_vec,el_vec);

% TODO 
% * [x] do for each mf -- kk
%

% bin volumes
vol_r_bin=4*pi*(r_c.^2).*diff(r_ed);

nr=cell(1,2);
for kk=1:2
    ttnr=NaN([size(az),nrBins]);  % preallocate - radial density profile at each grid point
    
    tk=vertcat(K{:,kk});     % collated k-vecs (cart-zxy)
    ts=zxy2sphpol(tk);      % collated k-vecs (spol)
    for ii=1:nAz
        for jj=1:nEl
            %%% get counts near each zone
            % this zone
            taz=az_vec(ii);
            tel=el_vec(jj);
            
            % bool-array to counts near this zone
            tb_in=inZone(ts(:,1),ts(:,2),taz,tel,dpsi);
            
            % vec-norms of counts in this zone
            tr_in=ts(tb_in,3);
            tNr=nhist(tr_in,{r_ed})';    % histogram - transpose to keep as row vecs
            
            % normalise to 1D-radial PDF
            tnr=tNr./vol_r_bin;      % radial volume
            tnr=tnr./(sum(tnr)*diff(r_ed));     % PDF: normalise by total and bin-size
            
            % TODO - gaussian convolution
            % smooth
            tnr=smooth(tnr,nsmooth);    % moving average
            
            % store
            ttnr(ii,jj,:)=tnr;
        end
    end
    nr{kk}=ttnr;
end

% plot radial distributions
figure();
for kk=1:2
    subplot(1,2,kk);
    p=[];   % initialise plots array
    hold on;
    for ii=1:nAz
        for jj=1:nEl
            if rand()<0.002     % decimate spol region to plot
                tp=plot(r_c,squeeze(nr{kk}(ii,jj,:)),...
                    'LineWidth',2,...
                    'DisplayName',num2str([ii,jj]));
                p=[p,tp];
            end
        end
    end
    box on;
    xlabel('r');
    ylabel('PDF');
    leg=legend(p);
    leg.Title.String='IDX: AZ EL';
end

%%% get distance at peak density
for kk=1:2
    
end