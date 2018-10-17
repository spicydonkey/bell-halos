function [E,E0] = g2toE(g2_corr,g2_anti)
%   spin-correlation evaluated from two-particle g2 correlation amplitudes
%   calculated between correlated and anti-corr pairs
%
%   E   : raw correlation coefficient as determined from g2
%   E0  : mode occupancy corrected correlation (infinite g2 limit)
%

P_corr = g2_corr./(g2_corr+g2_anti);
P_anti = 1-P_corr;

E = P_corr - P_anti;        % measured correlation coefficient

g2=g2_corr+g2_anti-1;           % effective two-particle corr amplitude
E_ideal_src=(g2-1)./(g2+1);      % ideal source 2-particle correlation
E0=E./E_ideal_src;               % mode occupancy corrected correlation

end