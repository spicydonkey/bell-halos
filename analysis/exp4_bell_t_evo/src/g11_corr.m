function g = g11_corr(k,dK)

% dK=dK*1e-4;      % to unit-scale

% CONFIG
dth_lim=[pi-0.2,pi];
n_dth=30;
% construct diff-angle bins for g2 histogram
dth_ed=linspace(dth_lim(1),dth_lim(2),n_dth);
dth=edge2cent(dth_ed);      % angle between k1 and k2 (0 for CL)
thbb=pi-dth;


az_k=0;
el_k=0;
alpha=pi/8;

K=boost_zxy(k,dK);
v_k=cellfun(@(x) inDoubleCone(x,az_k,el_k,alpha),K,'UniformOutput',false);
g2=g2_ang(v_k,dth_ed);

% fit amplitude
% Constraints: MU = 0, OFFSET = 1
fgauss.mdl='y~amp*exp(-1*x^2/(2*sigma^2))+1';
fgauss.cname={'amp','sigma'};
fgauss.fopt=statset('TolFun',10^-10,'TolX',10^-10,'MaxIter',10^6,'UseParallel',0);

sigma0=0.02;

% FIT
% g2_fit=cell(nfiles,3);
g2_fitmdl=fitnlm(thbb,g2,fgauss.mdl,[g2(end) ,sigma0],...
    'CoefficientNames',fgauss.cname,'Options',fgauss.fopt);

% get fit params
fparam=g2_fitmdl.Coefficients.Estimate;
famp=fparam(1);
fsig=fparam(2);

% g=g2(end);
g=famp;
end