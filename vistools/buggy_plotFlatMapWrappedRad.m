% debugging plotFlatMapWrappedRad
%   
%
% DKS 2019


%%% config data
% grid
az_lim = 0.9999*pi*[-1,1];
el_lim = 0.8*pi*0.5*[-1,1];

n_az = 10;
n_el = 5;

az = linspace_lim(az_lim,n_az);
el = linspace_lim(el_lim,n_el);

[AZ,EL] = ndgrid(az,el);

% image value
X = cos(AZ).*sin(EL);


%%% the bug
figure;

pp = plotFlatMapWrappedRad(AZ,EL,X,'rect','texturemap');
pp.DisplayName = 'plotFlatMapWrappedRad';

xlim(180*[-1,1]);
ylim(90*[-1,1]);

box on;


% markers for datapoints
hold on;

p_data = plot(rad2deg(AZ(:)),rad2deg(EL(:)),'g+','DisplayName','data');


% info
titlestr = 'Bug in plotFlatMapWrappedRad';
title(titlestr);

legend([p_data,pp])