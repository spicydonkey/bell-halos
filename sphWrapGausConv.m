function ff_filt = sphWrapGausConv(thphi,ff,sig,lim)
%
% FF_FILT = SPHWRAPGAUSSCONV(THPHI,FF,SIG,LIM)
%
% FF_FILT: Gaussian convolution filtered values on sphere
%
% THPHI: N-by-2 array of latitude-longitude coords on sphere
% FF: Nx1 array of function/map values on spherical coords to filter
% SIG: angular sigma for Gaussian profile
% LIM: limit (in standard score) from zero to evaluate Gaussian profile
% (i.e. extent of convolution)
%

th=thphi(:,1);
ph=thphi(:,2);

nZone=size(thphi,1);
ff_filt=zeros(nZone,1);
for ii=1:nZone
    z_psi=diffAngleSph(th,ph,th(ii),ph(ii))/sig;
    bool_lim=z_psi<lim;
    
    z_psi_lim=z_psi(bool_lim);  % diff angles in lim
    w_psi_lim=exp(-0.5*z_psi(bool_lim).^2);     % Gaussian weights of zones in lim
    ff_lim=ff(bool_lim);        % function values in lim
    
    % integrate weighted function values - convolution
    % NOTE: normalise for distribution of PSI
    
    tff_filt=sum(ff_lim.*w_psi_lim);
    % TODO - not normalised for zone distribution!
    
    ff_filt(ii)=tff_filt;
end