% Anomalous Rabi-flopping
%   Investigate non-flopping population oscillations in some regions of the
%   halo under a single-beam Raman process.
%
% DK Shin
%


% Load data
load('theta_rabi_20180131.mat');


% get sph-grid
n_sphgrid=size(az);
n_az=n_sphgrid(1);
n_el=n_sphgrid(2);


% select parts of halo sph-zones
el_max=0.8;
el_in=find(abs(el(1,:))<el_max);        % polar coords away from poles

idx_az=1:6:n_az;
idx_el=el_in(1):3:el_in(end);


% Plot oscillations
% TODO - do for mf=0 data after background subtraction
this_mf_idx=2;      % 1 for mf=0; 2 for mf=1 

data_color=distinguishable_colors(numel(idx_az));          % thetha
data_linestyle={'-','--',':','-.'};      % phi

figure('Name','rabi_cycle');
for ii=1:numel(idx_az)
    for jj=1:numel(idx_el)
        hold on;
        plot(1e6*traman,squeeze(P_mf{this_mf_idx}(idx_az(ii),idx_el(jj),:)),...
            'Color',data_color(ii,:),...
            'LineStyle',data_linestyle{mod(jj-1,numel(data_linestyle))+1},...
            'LineWidth',1.5);
        
    end
end
xlabel('$\tau$ ($\mu$s)');
ylabel('$P$');
box on;


%% Fit Rabi-cycle
this_mf_idx=2;      % mf=1 works well for simple model

tP=P_mf{this_mf_idx};

% preallocate Rabi-cycle params
rabiAmp=NaN(size(az));
rabiOmega=NaN(size(az));

% simple model
f_rabi='y~amp*sin(omega*x1/2)^2';

figure('Name','fit_rabi');

for ii=1:size(az,1)
    for jj=1:size(az,2)
        tp=squeeze(tP(ii,jj,:));
        
        % filter for any zone with NaN
        if sum(isnan(tp))>0
            continue
        end
        
        [pmax,I]=max(tp);     % max amplitude
        t_pmax=traman(I);               % time at max amplitude
        
        param0=[pmax,pi/t_pmax];
        tfitobj=fitnlm(traman,tp,f_rabi,param0,'CoefficientNames',{'amp','omega'});
        
        % get Rabi-cycle params
        rabiAmp(ii,jj)=tfitobj.Coefficients.Estimate(1);               
        rabiOmega(ii,jj)=tfitobj.Coefficients.Estimate(2);
        
        % debug
        if rand()>0.99
            hold on;
            ax=plot(1e6*traman,tp,'o');
            
            tfit=linspace(0,10e-6);
            pfit=feval(tfitobj,tfit);
            
            plot(1e6*tfit,pfit,'Color',ax.Color,...
                'LineWidth',1.5);
        end
    end
end
xlabel('$\tau$ ($\mu$s)');
ylabel('$P$');
box on;
axis tight;


% Plot Rabi-params around halo
figure('Name','rabi_params');
subplot(1,2,1);
plotFlatMapWrappedRad(az,el,rabiAmp);
cb=colorbar('southoutside');
cb.Label.String='Rabi amplitude';

subplot(1,2,2);
plotFlatMapWrappedRad(az,el,1e-6/(2*pi)*rabiOmega);
cb=colorbar('southoutside');
cb.Label.String='Rabi frequency [MHz]';

colormap('parula');


%% Two-photon Raman transition simplified
% Process unphysical raman amplitude
%   raman amplitude needs to be in [0,1]
% NOTE apply a saturation at 1.
% TODO
%   [ ] need to constrain fit to [0,1]
rabiAmpConst=rabiAmp;
rabiAmpConst(rabiAmpConst>1)=1;     % all fitted amps greater than 1 is set to 1

% calculated simple 2-level model characteristics
rabiOmega0_calc=sqrt(rabiAmpConst).*rabiOmega;
rabiDelta_calc=sqrt(1-rabiAmpConst).*rabiOmega;

figure('Name','rabi_param0_calc');
subplot(1,2,1);
plotFlatMapWrappedRad(az,el,1e-6/(2*pi)*rabiOmega0_calc);
cb=colorbar('southoutside');
cb.Label.String='\Omega_0 Resonant Rabi frequency [MHz]';

subplot(1,2,2);
plotFlatMapWrappedRad(az,el,1e-6/(2*pi)*rabiDelta_calc);
cb=colorbar('southoutside');
cb.Label.String='\Delta Detuning [MHz]';

colormap('viridis');


% characteristic detuning
%   $\Delta / \Omega_0$ is a dim-less characteristic
%   parameter for Rabi solution
rabiZ=rabiDelta_calc./rabiOmega0_calc;

figure('Name','rabi_char_detuning');
plotFlatMapWrappedRad(az,el,rabiZ,'eckert4','texturemap');
geoshow(rad2deg(el),rad2deg(az),rabiZ,'DisplayType','contour','LineColor','white');
cb=colorbar('southoutside');
cb.Label.String='\Delta / \Omega';
colormap('viridis');

