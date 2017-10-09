% Quantum mechanical prediction for spin-1/2 entangled pair
function corr = corr_qm(thetaAlice,thetaBob)
dtheta=thetaBob-thetaAlice;    
corr = -cos(dtheta);
end