% Process and plot Raman rotation in 2D plane

% NOTE: 
%   Run run_raman_mixing or load output
%   here we let mF=0 angles equal to +1's since the former has significant
%   errors in signal


%% configs
% process flags
do_graphics=1;
    save_graphics=0;

% signal filtering and preprocessing
numElevEdgesNaN=11;
% 2D gaussian filter (see '>> help fspecial')
g2d_hsize=[5,5];
g2d_sigma=0.5;

% 1D zonal histogram
nbin_dth=300;
g1d_hsize=15;
g1d_sigma=3;

%% main
lon=rad2deg(azim{1});
lat=rad2deg(elev{2});
nlatlonzones=numel(lon);

hFiltGauss=fspecial('gaussian',g2d_hsize,g2d_sigma);

% Process rotation angle signal
% NOTE - all signals used to arrive at the LOOP rotation angle was RAW and
% unsmoothed - TODO smoothing numbers in zone would be nice
thetaFilt=Theta{2};
nLoopConfig=size(thetaFilt,3);

for ii=1:nLoopConfig
    tTheta=thetaFilt(:,:,ii);
    
    % remove data near BEC/poles
    tTheta([1:numElevEdgesNaN,end-numElevEdgesNaN:end],:)=NaN;
    
    % 2D smoothing gaussian filter
    tTheta=imfilter(tTheta,hFiltGauss);     
    
    thetaFilt(:,:,ii)=tTheta;
end

% relative rotation angle
dTheta=thetaFilt;
dThetaFilt=zeros(size(dTheta));
for ii=1:nLoopConfig
    tDTheta=dTheta(:,:,ii);
    tDTheta=tDTheta-flip_bb(tDTheta,1);
    tDThetaFilt=imfilter(tDTheta,hFiltGauss);
    
    dTheta(:,:,ii)=tDTheta;
    dThetaFilt(:,:,ii)=tDThetaFilt;
end

% loop distribution in halo
sterPerZone=(4*pi)/nlatlonzones;
ed_dth=linspace(-pi,pi,nbin_dth);
ct_dth=ed_dth(1:end-1)+0.5*diff(ed_dth);
ster_dth=zeros(nLoopConfig,length(ct_dth));
g1d_filt=gaussFilter(g1d_hsize,g1d_sigma);
for ii=1:nLoopConfig
    tdth=dThetaFilt(:,:,ii);
    tdth=tdth(:);
    nn=histcounts(tdth,ed_dth);
    tster=nn*sterPerZone;
    tster=conv(tster,g1d_filt,'same');
    ster_dth(ii,:)=tster;
end

%% Plot - rotation angles
if do_graphics>0
    %%% ECKERT4 projection
    for ii=1:nLoopConfig
        tTheta=thetaFilt(:,:,ii);
        
        figure;
        plotFlatMap(lat,lon,tTheta,'eckert4');
        
        % misc
        ax=gca;
        
        tstr=sprintf('$m_F=+1$, Raman Amp = %0.2g',ampRaman_mf{2}(ii));
        title(tstr);
        colormap('magma');
        hcb=colorbar('Southoutside');
        hcb.Label.String='$\psi(\theta,\phi)$';
        hcb.TickLabelInterpreter='latex';
        hcb.Label.Interpreter='latex';
        ax.FontSize=12;
        
        drawnow;
        
        % save
        if save_graphics
            fnamethis=sprintf('theta_eck_mf1_%d.png',ii);
            saveas(gcf,fnamethis);
        end
    end
    
    %%% rectangle
    for ii=1:nLoopConfig
        tTheta=thetaFilt(:,:,ii);
        
        figure;
        plotFlatMap(lat,lon,tTheta,'rect');
        
        % misc
        ax=gca;
        
        tstr=sprintf('$m_F=+1$, Raman Amp = %0.2g',ampRaman_mf{2}(ii));
        title(tstr);
        
        xlim([-180,180]);
        ylim([-90,90]);
        
        colormap('magma');
        hcb=colorbar('Southoutside');
        hcb.Label.String='$\psi(\theta,\phi)$';
        hcb.TickLabelInterpreter='latex';
        hcb.Label.Interpreter='latex';
        ax.FontSize=12;
        
        drawnow;
        
        % save
        if save_graphics
            fnamethis=sprintf('theta_rect_mf1_%d.png',ii);
            saveas(gcf,fnamethis);
        end
    end
end

%% Plot - relative rotation angles
if do_graphics>0
    %%% ECKERT4 projection
    for ii=1:nLoopConfig
        tDTh=dTheta(:,:,ii);
        
        figure;
        plotFlatMap(lat,lon,tDTh,'eckert4');
        
        % misc
        ax=gca;
        
        tstr=sprintf('Raman Amp = %0.2g',ampRaman_mf{2}(ii));
        title(tstr);
        colormap('magma');
        hcb=colorbar('Southoutside');
        hcb.Label.String='$\Delta\psi(\theta,\phi)$';
        hcb.TickLabelInterpreter='latex';
        hcb.Label.Interpreter='latex';
        ax.FontSize=12;
        
        drawnow;
        
        % save
        if save_graphics
            fnamethis=sprintf('dtheta_eck_%d.png',ii);
            saveas(gcf,fnamethis);
        end
    end
    
    %%% rectangle
    for ii=1:nLoopConfig
        tDTh=dTheta(:,:,ii);
        
        figure;
        plotFlatMap(lat,lon,tDTh,'rect');
        
        % misc
        ax=gca;
        
        tstr=sprintf('Raman Amp = %0.2g',ampRaman_mf{2}(ii));
        title(tstr);
        
        xlim([-180,180]);
        ylim([-90,90]);
        
        colormap('magma');
        hcb=colorbar('Southoutside');
        hcb.Label.String='$\Delta\psi(\theta,\phi)$';
        hcb.TickLabelInterpreter='latex';
        hcb.Label.Interpreter='latex';
        ax.FontSize=12;
        
        drawnow;
        
        % save
        if save_graphics
            fnamethis=sprintf('dtheta_rect_%d.png',ii);
            saveas(gcf,fnamethis);
        end
    end
end

%% Plot - rotation angle distribution in halo
if do_graphics>0
    cc=distinguishable_colors(nLoopConfig);
    
    h_ss=figure(); hold on;

    pp=[];
    for ii=1:nLoopConfig
        ss=ster_dth(ii,:);
        
        figure(h_ss);
        tstr=sprintf('%0.2g',ampRaman_mf{2}(ii));
        pp(ii)=plot(ct_dth,ss,...
            'Color',cc(ii,:),'LineWidth',1.5,...
            'DisplayName',tstr);
    end
    figure(h_ss); hold off;
    
    
    % test angles - lines
    testAngles=[-pi/4,pi/4,3/4*pi];
    
    hold on;
    yylim=ylim;
    ylim(yylim);        % fix ylim
    plines=[];
    gray_col=0.8*ones(1,3);         % gray data points
    for ii=1:numel(testAngles)
        plines(ii)=line(testAngles(ii)*[1,1],yylim,...
            'Color',gray_col,'LineStyle','--','LineWidth',2);
        uistack(plines(ii),'bottom');
    end

    
    % annotation
    figure(h_ss);
    ax=gca;
    xlabel('$\Delta\psi$');
    ylabel('Solid angle in halo [sr]');
    box on;
    xlim([-pi,pi]);
    leg=legend(pp);
    leg.Title.String='Raman amp.';
    ax.FontSize=12;
    leg.FontSize=10;
    
    % save
    if save_graphics
        saveas(gcf,'dth_dist.png');
    end
end
