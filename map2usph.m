function q=map2usph(k)
% simple transformation to ~unit sphere by the ellipsoid_fit routine
% q=map2usph(k)
%
%

% fit ellipsoid using the complete dataset
K=vertcat(k{:});
[ec,er,ev]=ellipsoid_fit(circshift(K,-1,2),'');     % in XYZ frame

% transform to the unit-sphere!
q=cellfun(@(x)ellip2usph(x,ec,er,ev),k,'UniformOutput',false);

end