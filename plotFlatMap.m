function mapimg = plotFlatMap(lat,lon,Z,axproj,dispType)
% H = PLOTFLATMAP(LAT,LON,Z,AXPROJ)
%
% plot spherical map on a 2D plane
%
% LAT, LON should be in degrees
%
% NaN is transparent
%
%   dispType
%   	'texturemap'    no smooth
%       'surface'       smooth
%

% input checking
if ~exist('axproj','var') || isequal(axproj,[])
    warning('axproj is undefined. Default to eckert4 map.');
    axproj='eckert4';
end
if ~exist('dispType','var') || isequal(dispType,[])
    warning('dispType is undefined. Default to surface.');
    dispType='surface';
end

% set up axes
switch axproj
    case 'eckert4'
        axesm eckert4;
        framem; gridm;
        axis off;

    case 'rect'                
        xlim([-180,180]);
        ylim([-90,90]);
        
    otherwise
        warning('axproj should be set to either "eckert4" or "rect".');
end
% draw data like map
mapimg=geoshow(lat,lon,Z,'DisplayType',dispType);

% default annotations
colormap('magma');
% hcb=colorbar('Southoutside');
% hcb.TickLabelInterpreter='latex';
% hcb.Label.Interpreter='latex';

%% NaN to transparent
if dispType=='texturemap'
    mapimg.AlphaDataMapping='none';     % interpet alpha values as transparency values
    mapimg.FaceAlpha='texturemap';      % Indicate that the transparency can be different each pixel
    mapimg.AlphaData=double(~isnan(Z)); % set transparency matrix: NaN-->1, else-->0.
end

end
