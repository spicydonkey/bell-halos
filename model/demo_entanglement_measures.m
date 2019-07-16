% 2-qubit density matrix and entanglement measures
% DKS
% 2019-07-16


% config
sig_y = [0,-1i;1i,0];

v_psim = 1/sqrt(2)*[0;1;-1;0];
v_psip = 1/sqrt(2)*[0;1;1;0];
v_phip = 1/sqrt(2)*[1;0;0;1];
v_phim = 1/sqrt(2)*[1;0;0;-1];

rho_psim = v_psim'.*v_psim;
rho_psip = v_psip'.*v_psip;
rho_phip = v_phip'.*v_phip;
rho_phim = v_phim'.*v_phim;

rho_incoh = diag(ones(1,4))/4;


% psi+ with imperfections
p = 0.89        % fidelity

% rho = p* rho_psim + (1-p)*rho_incoh
rho = p* rho_psip + (1-p)*rho_incoh


% 1. entanglement of formation
spinflip = @(rho) kron(sig_y,sig_y)*conj(rho)*kron(sig_y,sig_y);
C = @(L) max([0, 2*L(1) - sum(L)])          % concurrence
h = @(x) -x*log2(x) - (1-x)*log2(1-x)
E_eof = @(C) h(0.5*(1 + sqrt(1 - C^2)))     % formula for entanglement of formation of arbitrary 2-qubit states

X = rho*spinflip(rho);
lambdas = sort(sqrt(eig(X)),'descend')

% equivalent way to evaluate eigenvalues:
% R = sqrt( sqrt(rho) .* rho_bar .* sqrt(rho));
% lambdas = sort(eig(R),'descend')


E = E_eof(C(lambdas))             % entanglement of formation


% plot
H = figure;
p_dm = bar3(rho);

ax=gca;

dm_basis_label = {'$\uparrow\uparrow$';'$\uparrow\downarrow$';'$\downarrow\uparrow$';'$\downarrow\downarrow$'};
set(ax,'XTickLabel',dm_basis_label);
set(ax,'YTickLabel',dm_basis_label);
zlabel('$P$');
titlestr = sprintf('Fidelity $= %0.2g$',p);
title(titlestr);