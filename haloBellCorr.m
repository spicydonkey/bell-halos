function [E]=haloBellCorr(N_halo_down,N_halo_up)
% Calculate Bell correlation E
% 
% [E] = haloBellCorr(N_halo)
%

% coincidence around halo
Jz_A=N_halo_up-N_halo_down;
N_A=N_halo_up+N_halo_down;

% entangled counterpart
% NOTE: assumes az-elev grid can be flipped to BB condition
Jz_B=flip_bb(Jz_A,1);
N_B=flip_bb(N_A,1);

% evaluate correlation
JJ=Jz_A.*Jz_B;
NN=N_A.*N_B;
E=mean(JJ,3)./mean(NN,3);

end