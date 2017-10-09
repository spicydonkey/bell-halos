% plot spherical map on a 2D plane
function h = plotFlatMap(lat,lon,Z,axesOpts)
% axes type
if isequal(axesOpts,'eckert4')
    axesm eckert4;
    framem; gridm;
    axis off;
end

h=geoshow(lat,lon,Z,'DisplayType','texturemap');

end
