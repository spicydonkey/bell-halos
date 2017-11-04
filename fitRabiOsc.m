function [afreq,afreq_se,xfit,yfit]=fitRabiOsc(amp,pp)
% AFREQ = FITRABIOSC(AMP, PP)
%
% evaluates Rabi amplitude frequency fitted from oscillation
%
% AMP: AOM amplitude
% PP: normalised original state population
%

% define amp Rabi oscillation model
aRabiFun=@(aOmega,amp) cos(aOmega*amp/2).^2;  % simple, ideal amplitude
coefname={'aOmega'};
param0=1;

% aRabiFun=@(p,amp) cos(0.5*p(1)*(sqrt(amp)-p(2))).^2;  % ideal amplitude, scale
% coefname={'aOmega','A'};
% param0=[1,0.33];

% aRabiFun=@(p,amp) (amp>sqrt(p(2))).*cos(0.5*p(1)*(sqrt(amp)-p(2))).^2 + (amp<=sqrt(p(2)))*NaN;  % ideal amplitude, scale
% coefname={'aOmega','A'};
% param0=[1,0.33];

fopts=statset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e2);

fit_rabi=fitnlm(amp,pp,aRabiFun,param0,...
    'CoefficientNames',coefname,...
    'Options',fopts);

% get amplitude Rabi frequency from fit
afreq=fit_rabi.Coefficients.Estimate;
afreq_se=fit_rabi.Coefficients.SE;

% eval fitted oscillation
xfit=linspace(min(0,min(amp)),max(0.55,max(amp)),1000);
yfit=feval(fit_rabi,xfit);

end