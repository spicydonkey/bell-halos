function [g2,dk,g2mdl,h]=summary_disthalo_g2(K,dk_ed,bfastrun,doplot,dofit,verbose)
% Summarise distinguishable (2-species) halos with various g2 corrleation
% functions
% Useful diagnostic to check how well halos were processed
%
%   Simplest usage:
%       [g2,dk,h]=summary_disthalo_g2(K)
%
% TODO
% [ ] documentation
%

tstart=tic;


%% parse inputs
if ~exist('dk_ed','var')
    dk_ed={};           % default rel-vector bin edges
end

if ~exist('bfastrun','var')
    bfastrun=false;     % default fully evaluates uncorrelated G2 for normalisation
end

if ~exist('doplot','var')
    doplot=true;        % default is to plot
end

if ~exist('dofit','var')
    dofit=true;         % default is to fit with gaussian
end

if ~exist('verbose','var')
    verbose=1;          % default is verbose
end


%% configs
%%% set-up g2 histogramming bins
if isempty(dk_ed)
    % set up default domain (symmetric in Cart axes)
    nbin_g2=29;             % num. bins
    dk_lim=[-0.2,0.2];      % dim limits
    dk_ed_vec=linspace(dk_lim(1),dk_lim(2),nbin_g2+1);  % bin edges per dim
    
    % set up default 3D bins
    dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
end

%%% parse bin edges
dk_cent=cellfun(@(ed) ed(1:end-1)+0.5*diff(ed),dk_ed,'UniformOutput',false);    % bin centers
[~,idx_dk0]=cellfun(@(c) min(abs(c)),dk_cent);      % idx bin nearest zero
[dk{1},dk{2},dk{3}]=ndgrid(dk_cent{:});             % grid of dk centers


%%% G2 smoothing
filter_toggle=true;       % turn filter on/off
filter_type='gaussian';
filter_size=5;
filter_sd=1.5;


%% g2 analysis
% get scattering type
n_halo=size(K,2);
n_g2_type=2*n_halo-1;   % in fact all permutations gives n_halo^2 (but symmetry!)

% preallocate
g2=cell(1,n_g2_type);
G2s=cell(1,n_g2_type);
G2n=cell(1,n_g2_type);

counter=1;

%%% Single halo: (0,0)/(1,1)
for mm=1:n_halo
    if bfastrun
        [g2{counter},G2s{counter},G2n{counter}]=g2_bb_fast(K(:,mm),dk_ed,0.05);
    else
        [g2{counter},G2s{counter},G2n{counter}]=g2_bb(K(:,mm),dk_ed);
    end
    if filter_toggle
        % filter g2
        G2s{counter}=smooth3(G2s{counter},filter_type,filter_size,filter_sd);
        G2n{counter}=smooth3(G2n{counter},filter_type,filter_size,filter_sd);
        g2{counter}=G2s{counter}./G2n{counter};
    end
    counter=counter+1;
end


%%% X-species: (0,1)
if n_halo==2
    if bfastrun
        [g2{counter},G2s{counter},G2n{counter}]=g2x_bb_fast(K,dk_ed,0.05);
    else
        [g2{counter},G2s{counter},G2n{counter}]=g2x_bb(K,dk_ed);
    end
    if filter_toggle
        % filter g2
        G2s{counter}=smooth3(G2s{counter},filter_type,filter_size,filter_sd);
        G2n{counter}=smooth3(G2n{counter},filter_type,filter_size,filter_sd);
        g2{counter}=G2s{counter}./G2n{counter};
    end
end


%% Fit g2 - gaussian model
g2_fit=cell(1,n_g2_type);
if dofit
    %   3D gaussian model
    g2_gauss3_mdl='y~1+amp*exp(-0.5*(((x1-mu1)/sig1)^2+((x2-mu2)/sig2)^2+((x3-mu3)/sig3)^2))';
    g2_gauss3_mdl_param={'amp','mu1','mu2','mu3','sig1','sig2','sig3'};
    
    g2mdl=cell(1,n_g2_type);
    
    % rel-vec centers as 1D list of vectors
    X=cellfun(@(d) d(:),dk,'UniformOutput',false);
    X=cat(2,X{:});
    
    for ii=1:n_g2_type
        % get evaluated g2
        Y=g2{ii}(:);
        
        % estimate params
        tamp0=g2{ii}(idx_dk0(1),idx_dk0(2),idx_dk0(3));
        tmu0=[0,0,0];
        tsig0=[0.03,0.03,0.03];
        tparam0=[tamp0,tmu0,tsig0];
        
        %%% fit model
        tfopts=statset('Display','off');
        tg2mdl=fitnlm(X,Y,g2_gauss3_mdl,tparam0,'CoefficientNames',g2_gauss3_mdl_param,...
            'Options',tfopts);
        tparam_fit=tg2mdl.Coefficients.Estimate;
        
        g2mdl{ii}=tg2mdl;
    end
    
    % evaluate fitted model
    n_bins_fit=100;
    kk_1d=cellfun(@(c) linspace(min(c),max(c),n_bins_fit),dk_cent,'UniformOutput',false);
    [~,idx_dk0_fit]=cellfun(@(c) min(abs(c)),kk_1d);      % idx bin nearest zero
    [kk{1},kk{2},kk{3}]=ndgrid(kk_1d{:});
    
    ksqueezed=cellfun(@(k) k(:),kk,'UniformOutput',false);
    ksqueezed=cat(2,ksqueezed{:});
    % how to go back from 1D flattened to 3D-grid?
    
    for ii=1:n_g2_type
        tg2_squeezed=feval(g2mdl{ii},ksqueezed);
        g2_fit{ii}=reshape(tg2_squeezed,size(kk{1}));     % form as 3D grid;
    end
end


%% Plot
[cc,clight,cdark]=palette(3);

h=[];       % initialise fig handle
if doplot
    
    titlestr={'(0,0)',...
        '(1,1)',...
        '(0,1)',...
        };

    dispname={'$x$','$y$','$z$'};     % axis
    mark_type={'o','^','d'};    % markers for each axis
    
    
    h=figure();
    for ii=1:n_g2_type
        subplot(1,n_g2_type,ii);
        
        perm_ord=[2,3,1];
        temp_ord=[1,2,3];
        temp_g2_perm=g2{ii};    % temporary var to hold dimension permuted g2
        temp_gfit_perm=g2_fit{ii};  % fitted g2 model
        
        p=NaN(3,1);     % fig objects
        pfit=NaN(3,1);
        for jj=1:3
            hold on;
            % permute dims
            %   NOTE: first permutation turns ZXY --> XYZ
            temp_ord=temp_ord(perm_ord);
            temp_g2_perm=permute(temp_g2_perm,perm_ord);
            
            p(jj)=plot(dk_cent{temp_ord(1)},temp_g2_perm(:,idx_dk0(temp_ord(2)),idx_dk0(temp_ord(3))),...
                'MarkerEdgeColor',cc(jj,:),'MarkerFaceColor',clight(jj,:),...
                'LineStyle','none','Marker',mark_type{jj},'DisplayName',dispname{jj});
            
            % plot fitted
            if dofit
                temp_gfit_perm=permute(temp_gfit_perm,perm_ord);
                pfit(jj)=plot(kk_1d{ii},temp_gfit_perm(:,idx_dk0_fit(temp_ord(2)),idx_dk0_fit(temp_ord(3))),...
                    'Color',cc(jj,:),'LineStyle','-');
                uistack(pfit(jj),'bottom');
            end
        end
        hold off;
        
        % annotate
        legend(p);
        title(titlestr{ii});
        xlabel('$\Delta k$');
        ylabel('$g^{(2)}_{BB}$');
        box on;
        grid on;
    end
end


%% DONE
if verbose>0
    toc(tstart);    % report elapsed time
end

end