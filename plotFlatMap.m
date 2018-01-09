function h = plotFlatMap(lat,lon,Z,axesOpts)
% H = PLOTFLATMAP(LAT,LON,Z,AXESOPTS)
%
% plot spherical map on a 2D plane
%
% LAT, LON should be in degrees
%

if ~exist('axesOpts','var')
    warning('axesOpts is undefined. Default to rectangle map.');
    axesOpts='rect';
end

% set up axes
switch axesOpts
    case 'eckert4'
        axesm eckert4;
        framem; gridm;
        axis off;

    case 'rect'                
        xlim([-180,180]);
        ylim([-90,90]);
        
    otherwise
        warning('axesOpts should be set to either "eckert4" or "rect".');
end
% draw data like map
% h=geoshow(lat,lon,Z,'DisplayType','texturemap');  % NaN as minimum value; smooth
h=geoshow(lat,lon,Z,'DisplayType','surface');   % NaN transparent; no smooth


% default annotations
colormap('magma');
% hcb=colorbar('Southoutside');
% hcb.TickLabelInterpreter='latex';
% hcb.Label.Interpreter='latex';

end
