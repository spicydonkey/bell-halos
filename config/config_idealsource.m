% Default configs for ideal source
%

clearvars configs;

%% HALO
configs.halo{1}.string='$m_F=0$';

configs.halo{2}.string='$m_F=1$';

%% Spherical zones
configs.zone.nazim=100;
configs.zone.nelev=50;

configs.zone.binmethod=1;
configs.zone.binwidth=0.05;


%% g2 corr
% % clear already predefined confgs
% configs.corr={};
% 
% % 1) X-halo Cart BB
% configs.corr{1}.type.comp=[1,2];           % components to analysis: cross halo 1,2
% configs.corr{1}.type.coord='cart';         % Cartesian (ZXY)
% configs.corr{1}.type.opt='BB';             % BB / CL
% configs.corr{1}.lim=0.2*repmat([-1,1],[3,1]);
% configs.corr{1}.nBin=15*[1,1,1];   % number of bins
% 
% % 2) Single-halo cart BB - m_J=0
% configs.corr{2}.type.comp=1;
% configs.corr{2}.type.coord='cart';
% configs.corr{2}.type.opt='BB';
% configs.corr{2}.lim=0.2*repmat([-1,1],[3,1]);
% configs.corr{2}.nBin=11*[1,1,1];   % number of bins
% 
% % 3) Single-halo cart BB - m_J=1
% configs.corr{3}.type.comp=2;
% configs.corr{3}.type.coord='cart';
% configs.corr{3}.type.opt='BB';
% configs.corr{3}.lim=0.2*repmat([-1,1],[3,1]);
% configs.corr{3}.nBin=11*[1,1,1];   % number of bins