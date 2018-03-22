function z_padded=pad_sphgrid_poles(el,z,el_pol,pval)
%Pads polar regions of a sph-grid map with special value.
%
%   el: elev vals of sph-grid
%   z: sph-grid map
%   el_pol: elev polar limits [rad]
%   pval: value to pad regions in pole (default: NaN)
%
%   z_padded: sph-grid map with polar vals padded
%

if ~exist('pval','var')
    pval=NaN;
end

b_pole=(abs(el)>el_pol);      % bool to region in poles to pad

z_padded=z;
z_padded(b_pole)=pval;

end