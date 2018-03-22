function h = plotFlatMap(lat,lon,Z,axproj,dispType)
% H = PLOTFLATMAP(LAT,LON,Z,AXPROJ)
%
% plot spherical map on a 2D plane
%
% LAT, LON should be in degrees
%
%   dispType
%   	'texturemap'    NaN as minimum value; smooth
%       'surface'       NaN transparent; no smooth
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
h=geoshow(lat,lon,Z,'DisplayType',dispType);

% default annotations
colormap('magma');
% hcb=colorbar('Southoutside');
% hcb.TickLabelInterpreter='latex';
% hcb.Label.Interpreter='latex';

end
