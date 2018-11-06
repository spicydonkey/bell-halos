function [E,E0] = g_to_E(g)
%G_TO_E evaluates E correlator from g2 vector (g_++,g_--,g_+-)
%
%   g: vector of g2 correlations ordered: ++, --, +-,
%   
%   E: corr coefficient
%   E0: corrected for finite mode occupancy
%

[E,E0] = g2toE(0.5*(g(:,1)+g(:,2)),g(:,3));
end