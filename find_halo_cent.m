function [cent_best,dr_best] = find_halo_cent(vs)
% locate center of halo (mf=0) by minimising rms width
%
%   vs: Nx3 ZXY array of counts 
%       shots are collated
%   
%   cent_best: 1x3 best centre to minimise rms width
%   dr_best: minimum rms of vec norms
%

fun_normrms = @(x) std(zxy2rdist(vs,x));      % rms function
cent0 = [0,0,0];       % initial param at zero
fopts = optimset('Display','iter','PlotFcns',@optimplotfval,...
    'TolX',1e-6,'TolFun',1e-6);

[cent_best,dr_best] = fminsearch(fun_normrms,cent0,fopts);

end