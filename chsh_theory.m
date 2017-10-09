% CHSH value 'S' for quantum correlation
function s_chsh = chsh_theory(thetas)
% local rotation angles: a1,a2,b1,b2
a1=thetas(1);
a2=thetas(2);
b1=thetas(3);
b2=thetas(4);

s_chsh = abs(corr_qm(a1,b1)-corr_qm(a1,b2)+corr_qm(a2,b1)+corr_qm(a2,b2));
end