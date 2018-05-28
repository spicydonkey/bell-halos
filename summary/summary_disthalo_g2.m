function [g2,dk,h]=summary_disthalo_g2(K,bfastrun,doplot,verbose)
% Summarise distinguishable (2-species) halos with various g2 corrleation
% functions
% Useful diagnostic to check how well halos were processed
%
% [g2,dk,h]=summary_disthalo_g2(K,bfastrun)
%
% TODO
% [ ] documentation

tstart=tic;

%% parse inputs
if ~exist('bfastrun','var')
    bfastrun=false;     % default fully evaluates uncorrelated G2 for normalisation
end

if ~exist('doplot','var')
    doplot=true;        % default is to plot
end

if ~exist('verbose','var')
    verbose=1;          % default is verbose
end


%% configs
%%% g2
nbin_g2=30;
dk_lim=[-0.2,0.2];
dk_ed_vec=linspace(dk_lim(1),dk_lim(2),nbin_g2);
dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
[~,idx_dk0]=min(abs(dk_cent_vec));
dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};      % symmetric in cart axis
dk_cent={dk_cent_vec,dk_cent_vec,dk_cent_vec};
% dk=ndgrid(dk_cent{:});      % grid of dk centers
[dk{1},dk{2},dk{3}]=ndgrid(dk_cent{:});      % grid of dk centers

%%% G2 smoothing
filter_toggle=true;       % turn filter on/off
filter_type='gaussian';
filter_size=5;
filter_sd=1.5;

%% g2 analysis
% preallocate 
g2=cell(1,3);
G2s=cell(1,3);
G2n=cell(1,3);
counter=1;

%%% Single halo: (0,0)/(1,1)
for mm=1:2
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


%% plot
h=[];       % initialise fig handle
if doplot
    
    titlestr={'(0,0)',...
        '(1,1)',...
        '(0,1)',...
        };
    linestyle={'o-','^-','*-'};
    dispname={'X','Y','Z'};
    
    % h=NaN(3,1);
    h=figure();
    for ii=1:3
        %     h(ii)=figure();
        subplot(1,3,ii);
        
        temp_g2_perm=g2{ii};    % temporary var to hold dimension permuted g2
        
        p=NaN(3,1);     % fig objects
        for jj=1:3
            hold on;
            temp_g2_perm=permute(temp_g2_perm,[2,3,1]);
            
            p(jj)=plot(dk_cent_vec,temp_g2_perm(:,idx_dk0,idx_dk0),...
                linestyle{jj},'DisplayName',dispname{jj});
            % NOTE: dk_cent_vec needs to be symmetric in cart-dims otherwise
            % need to permute this too to match against g2 dim
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