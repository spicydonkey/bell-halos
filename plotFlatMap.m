function h = plotFlatMap(lat,lon,Z,axesOpts)
% H = PLOTFLATMAP(LAT,LON,Z,AXESOPTS)
%
% plot spherical map on a 2D plane
%
% LAT, LON should be in degrees
%

% TODO
% [] tidy switch cases below
% [] wrap azimuthal angle (currently chopped)
% [] create a wrapper for current zone convention (azim,elev) in rads
%

% axes type
switch axesOpts
    case 'eckert4'
        axesm eckert4;
        framem; gridm;
        axis off;
        
        h=geoshow(lat,lon,Z,'DisplayType','texturemap');
    case 'rect'
        xlim([-180,180]);
        ylim([-90,90]);
        
        h=geoshow(lat,lon,Z,'DisplayType','texturemap');
    otherwise
        warning('axesOpts should be set to either "eckert4" or "rect".');
end

% default annotations
colormap('magma');
hcb=colorbar('Southoutside');
hcb.TickLabelInterpreter='latex';
hcb.Label.Interpreter='latex';

end
