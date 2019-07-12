function yhat = Bcorr_mdl(beta,x)
%BCORR_MDL is a a Bell triplet correlator model.
%   YHAT = BCORR_MDL(BETA,X) gives the predicted values of the diagonal
%   correlator, YHAT, as a function of Bell triplet's purity and detuning
%   of Rabi oscillation (a vector of parameters, BETA), and the rotation
%   angle, X.
%   BETA is a 2-length vector: [p, alpha]
%       purity: 0 <= p <= 1
%       polar angle of rotation axis: 0 <= alpha (a) <= pi
%
%   The model form is:
%   y = -p * (cos(a)^2 + sin(a)^2*cos(x)).^2 ...
%       + p * ( (sin(a)*cos(a)*(1-cos(x))).^2 + (sin(a)*sin(x)).^2)
%
%
%   2019-07-12
%   DKS

p = beta(1);        % purity
a = beta(2);        % alpha (polar angle of rotation axis)

yhat = -p * (cos(a)^2 + sin(a)^2*cos(x)).^2 ...
      + p * ( (sin(a)*cos(a)*(1-cos(x))).^2 + (sin(a)*sin(x)).^2);