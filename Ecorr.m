function E = Ecorr(N00,N11,N01,N10)
% ECORR correlator formula
% DKS
% 20181019

E = (N00 + N11 - N01 - N10)./(N00 + N11 + N01 + N10);

end