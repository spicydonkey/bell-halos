function E = Emixed(th1,th2)
%   correlation from rotating mixed state of anti-correlated pair
%   

    E = -P_anti(th1,th2) + P_corr(th1,th2);
end


function p = P_anti(th1,th2)
%   probability of detecting anti-correlation
%

    p = th2p(th1).*th2p(th2) + (1-th2p(th1)).*(1-th2p(th2));
end


function p = P_corr(th1,th2)
%   probability of detecting correlation
%

    p = (1-th2p(th1)).*th2p(th2) + th2p(th1).*(1-th2p(th2));
end


function p = th2p(th)
%   probability of finding original polarization upon rotation
%

    p = cos(th/2).^2;
end