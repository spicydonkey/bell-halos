function [E]=haloBellCorr(N_halo_down,N_halo_up)
% Calculate Bell correlation E
% 
% [E] = haloBellCorr(N_halo)
%

s=size(N_halo_down);

% coincidence around halo
Jz_A=N_halo_up-N_halo_down;
N_A=N_halo_up+N_halo_down;

% entangled counterpart
% NOTE: assumes az-elev grid can be flipped to BB condition
% below (B) values are BB counterpart to above (A)
Jz_B=flip_bb(Jz_A,1);
N_B=flip_bb(N_A,1);

%%% evaluate Bell correlation
b_nan=zeros(s);   % initialise ignorable events
JJ=Jz_A.*Jz_B;
NN=N_A.*N_B;

% pair-detection filter
% get 4-mode detection events for which there is exactly 2 atoms with BB momenta
b_pair_event=((N_A==1)&(N_B==1));
b_nan=b_nan|~b_pair_event;   % this line eliminates all obviously non-ideal events from contributing to the corr analysis

% ensure any event of nan is shared across the two (pseudo)coincidences
b_nan=b_nan|(isnan(JJ)|isnan(NN));      % update nan events
JJ(b_nan)=NaN;      % update (pseudo)coincidences
NN(b_nan)=NaN;

E=mean(JJ,3,'omitnan')./mean(NN,3,'omitnan');   % bell correlation

end