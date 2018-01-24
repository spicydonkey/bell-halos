function [g2_bb_00,g2_bb_11,g2_bb_01] = summary_disthalo_g2dist(N,Az,El,verbose)
%Summarise BB-pair g2 across all regions
%
%   N: 1x2 cell; count distribution around halo
%   Az: azim sph-grid
%   El: elev sph-grid
%   verbose: flag for summary plots
%
%   Usage:
%   
%
%   TODO
%   [ ] document code
%   


%%% input check
if nargin<4
    verbose=0;  % quiet by default
end


%%% calculate g2
g2_bb_00=g2_ncorr(N{1},flip_bb(N{1},1));
g2_bb_11=g2_ncorr(N{2},flip_bb(N{2},1));
g2_bb_01=g2_ncorr(N{1},flip_bb(N{2},1));

N_p=N{1}+N{2};	% only momenta resolved
g2_bb_comb=g2_ncorr(N_p,flip_bb(N_p,1));


%%% g2-distribution
% g2 saturated at 3.5 SD
if verbose>0
    % 1. BB; (0,0)
    figure;
    y=g2_bb_00;
    plotFlatMapWrappedRad(Az,El,y,'eckert4');
    cb=colorbar('southoutside');
    cb.Label.String='g^{(2)}';
    title('BB; (0,0)');
    yavg=mean(y(:),'omitnan');
    ystd=std(y(:),'omitnan');
    caxis([0,yavg+3.5*ystd]);
    
    % 2. BB; (1,1)
    figure;
    y=g2_bb_11;
    plotFlatMapWrappedRad(Az,El,y,'eckert4');
    cb=colorbar('southoutside');
    cb.Label.String='g^{(2)}';
    title('BB; (1,1)');
    yavg=mean(y(:),'omitnan');
    ystd=std(y(:),'omitnan');
    caxis([0,yavg+3.5*ystd]);
    
    % 3. BB; (0,1)
    figure;
    y=g2_bb_01;
    plotFlatMapWrappedRad(Az,El,y,'eckert4');
    cb=colorbar('southoutside');
    cb.Label.String='g^{(2)}';
    title(' BB; (0,1)');
    yavg=mean(y(:),'omitnan');
    ystd=std(y(:),'omitnan');
    caxis([0,yavg+3.5*ystd]);
    
    % 4. BB; m_F combined
    figure;
    y=g2_bb_comb;
    plotFlatMapWrappedRad(Az,El,y,'eckert4');
    cb=colorbar('southoutside');
    cb.Label.String='g^{(2)}';
    title('BB; $m_F$ combined');
    yavg=mean(y(:),'omitnan');
    ystd=std(y(:),'omitnan');
    caxis([0,yavg+3.5*ystd]);
    
    % 5. BB; (0,0) + (0,1)
    figure;
    y=g2_bb_00+g2_bb_01;
    plotFlatMapWrappedRad(Az,El,y,'eckert4');
    cb=colorbar('southoutside');
    cb.Label.String='g^{(2)}';
    title('BB; (0,0) + (0,1)');
    yavg=mean(y(:),'omitnan');
    ystd=std(y(:),'omitnan');
    caxis([0,yavg+3.5*ystd]);
    
    % 6. BB; (1,1) + (0,1)
    figure;
    y=g2_bb_11+g2_bb_01;
    plotFlatMapWrappedRad(Az,El,y,'eckert4');
    cb=colorbar('southoutside');
    cb.Label.String='g^{(2)}';
    title('BB; (1,1) + (0,1)');
    yavg=mean(y(:),'omitnan');
    ystd=std(y(:),'omitnan');
    caxis([0,yavg+3.5*ystd]);
end

end