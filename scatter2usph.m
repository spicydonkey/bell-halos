function [kk,ecent,erad,evecs]=scatter2usph(k)
% simple unit sphere mapping by ellipsoid fit to distorted scatter
%
% K: Nx1 cell-array (all shots) of ZXY counts
%
% 
% Fits to the whole dataset and "undistorts" all data by mappint go usphere based on the fitted
% ellipsoid
%

[ecent,erad,evecs]=ellipsoid_fit(circshift(vertcat(k{:}),-1,2),'');
kk=cellfun(@(x)ellip2usph(x,ecent,erad,evecs),k,'UniformOutput',false);

end