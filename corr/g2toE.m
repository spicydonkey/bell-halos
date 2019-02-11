function [E,E0] = g2toE(g2_corr,g2_anti)
% G2TOE evaluates spin-correlation from g2
%   calculated between correlated and anti-corr pairs
%
%   g2_corr: ++/-- 
%   g2_anti: +-/-+
%
%   E   : raw correlation coefficient as determined from g2
%   E0  : mode occupancy corrected correlation (infinite g2 limit)
%

P_corr = g2_corr./(g2_corr+g2_anti);
P_anti = 1-P_corr;

E = P_corr - P_anti;        % measured correlation coefficient

g2=0.5*(g2_corr+g2_anti);       % spin-integrated two-particle corr amp
E_ideal_src=1-1./g2;             % correlator amplitude for pair source

E0=E./E_ideal_src;               % mode occupancy corrected correlation

end