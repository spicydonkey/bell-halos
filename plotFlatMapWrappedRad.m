function h = plotFlatMapWrappedRad(Az,El,Z,axproj)
% H = PLOTFLATMATWRAPPEDRAD(AZ,EL,Z,AXPROJ)
%
% AZ: 2d-grid of azimuthal angles (rad); ordered and evenly distributed
% in [-pi,pi); must be unique modulo 2*pi;
% EL: 2d-grid of elevation angles (meas from XY plane) (rad); ordered 
% and in range [-pi/2,pi/2]; must be unique;
% Z: data defined on grid
% AXPROJ: eckert4 or rect
%

% input checking
if ~exist('axproj','var')
    warning('axproj is undefined. Default to rectangle map.');
    axproj='rect';
end

% wrap azimuthal angle to close surface
Az(end+1,:)=Az(1,:)+2*pi;
El(end+1,:)=El(1,:);
Z(end+1,:)=Z(1,:);

h=plotFlatMap(rad2deg(El),rad2deg(Az),Z,axproj);

end