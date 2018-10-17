%% Halo centering sensitivity analysis
%   
% DKS
% 20181015

% nparam=1;     % DEBUGGING

%% CONFIG
flag_analysis_3d=true;       % true for 3D; false for 1-bin g2
flag_g2_fastrun=false;      % fast-run gives noisy g2 so fit doesn't work

% some constants
dk_bb=0.04;     % back-to-back pair correlation length (rms width) [norm unit]


%% load data
if ~exist('flag_load_override','var') || ~flag_load_override
    path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';
    fname_mat='exp5_1_20180813_1.mat';
    path_data=fullfile(path_dir,fname_mat);
    
    load(path_data);
end

%% define g2 analysis
lim_Dk=1*dk_bb*[-1,1];          % limits for shift [normed]
n_Dk=19;                        % n-points to scan (MUST BE ODD!)
Dk_0=linspace(lim_Dk(1),lim_Dk(2),n_Dk);     % shift in halo center in each dim (symm)
idx_0=find(Dk_0==0);        % index to Dk=0

if ~flag_analysis_3d
    %%% 1-bin
    dk_ed={dk_bb*[-1,1],dk_bb*[-1,1],dk_bb*[-1,1]};   % define BB pair single-pixel    
else
    %%% full 3D analysis
    lim_dk=0.2;     % original: 0.2
    n_bins=29;      % original: 29
    dk_ed_vec=linspace(-lim_dk,lim_dk,n_bins+1);
    dk_cent_vec=dk_ed_vec(1:end-1)+0.5*diff(dk_ed_vec);
    [~,idx_dk0]=min(abs(dk_cent_vec));    % g2 bin-idx nearest 0
    dk_ed={dk_ed_vec,dk_ed_vec,dk_ed_vec};
    dk_grid_size=cellfun(@(k) length(k)-1,dk_ed);
end

%% evaluate g2 while shifting center
g2_shift=cell(nparam,1);        % scalar "amplitude" of g2 (3D-->fitted; 1-bin-->BB value)
g2_3d=cell(nparam,1);           % g2 function for all params
g2mdl=cell(nparam,1);           % fitted g2 model
% NOTE: Cell structure indexing:
%   {TAU} {ZXY} {DK, SPINSPIN}
%   TODO
%       [ ] {TAU} {ZXY, DK, SPINSPIN}

for ii=1:nparam
    tk=k_par{ii};
    
    % allocate temporary vars for storage
    tg2_3d=cell(3,n_Dk,3);          % DIMS: ZXY, DK, SPINSPIN
    tg2mdl=cell(3,n_Dk,3);
    tg2_shift=NaN(3,n_Dk,3);
    
    for jj=1:3          % ZXY dim
        for kk=1:n_Dk   % DK center shift
            tDk=circshift([Dk_0(kk),0,0],jj-1);     % get new origin shift
            
            % evaluate g2 at near BB
            if ~flag_analysis_3d
                %%% 1D g2 at BB (single-bin)
                tg2_shift(jj,kk,1)=g2_bb(boost_zxy(tk(:,1),tDk),dk_ed);       % (up,up)
                tg2_shift(jj,kk,2)=g2_bb(boost_zxy(tk(:,2),tDk),dk_ed);       % (down,down)
                
                ttk=tk;
                ttk(:,1)=boost_zxy(tk(:,1),tDk);
                tg2_shift(jj,kk,3)=g2x_bb(ttk,dk_ed);           % (up,down): shift only UP
            else
                %%% g2 in 3D. get fitted amplitude
                tk_shift=tk;
                tk_shift(:,1)=boost_zxy(tk(:,1),tDk);       % shift mJ=1
                % upup and updown
                [tg2,~,ttg2mdl]=summary_disthalo_g2(tk_shift,dk_ed,flag_g2_fastrun,0,1,0);
                tg2_shift(jj,kk,1)=ttg2mdl{1}.Coefficients.Estimate(1);       % UP-UP
                tg2_shift(jj,kk,3)=ttg2mdl{3}.Coefficients.Estimate(1);       % UP-DOWN
                tg2_3d{jj,kk,1}=tg2{1};
                tg2_3d{jj,kk,3}=tg2{3};
                tg2mdl{jj,kk,1}=ttg2mdl{1};
                tg2mdl{jj,kk,3}=ttg2mdl{3};
                
                % downdown
                tk_shift=tk;
                tk_shift(:,2)=boost_zxy(tk(:,2),tDk);       % shift mJ=0
                [tg2,~,ttg2mdl]=summary_disthalo_g2(tk_shift,dk_ed,flag_g2_fastrun,0,1,0);     
                tg2_shift(jj,kk,2)=ttg2mdl{2}.Coefficients.Estimate(1);       % DOWN-DOWN
                tg2_3d{jj,kk,2}=tg2{2};                
                tg2mdl{jj,kk,2}=ttg2mdl{2};
                
                %%% check acceptibility of fit
                %   peak location (mu) from 0 <~ corr length
                %   spread should not be >> than orig corr length
                for ll=1:3
                    % get fitted params
                    tparam=tg2mdl{jj,kk,ll}.Coefficients.Estimate;
                    tmu=tparam(2:4);
                    tsig=tparam(5:7);
                    % check if ok
                    if vnorm(tmu')>2*dk_bb || geomean(abs(tsig))>2*dk_bb
                        % NOT OK
                        tg2_shift(jj,kk,ll)=NaN;
                    end
                end
            end
        end
    end
    % store analysis results
    g2_shift{ii}=tg2_shift;
    if flag_analysis_3d
        g2_3d{ii}=tg2_3d;
        g2mdl{ii}=tg2mdl;
    end
end

% evaluate g2 normed wrt Dk=0
g2_shift_norm=cellfun(@(x) x./x(:,idx_0,:),g2_shift,'UniformOutput',false);

%% data vis
% configs
[c,cl,cd]=palette(3);
line_sty={'-','--',':'};
mark_typ={'o','s','^'};
% str_ss={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
% str_dim={'$z$','$x$','$y$'};

hfig=cell(nparam,1);
for ii=1:nparam
    %%
    tT=1e6*par_T(ii);
    figname=sprintf('g2shifted_T%0.2g',tT);
    
    hfig{ii}=figure('Name',figname,'Units','normalized',...
        'Position',[0.2,0.2,0.2,0.3],'Renderer','painters');    
    hold on;
    for jj=1:3      % Spin-Spin corr type
        for kk=1:3  % offset Cartesian dim
            % g2
            %             plot(Dk_0,g2_shift{ii}{jj}(:,kk),...
            %                 'Color',c(kk,:),'LineStyle',line_sty{jj},...
            %                 'LineWidth',1.5);
            
            % normalised to Dk=0
            plot(Dk_0,g2_shift_norm{ii}(kk,:,jj),...
                'Color',c(kk,:),'LineStyle',line_sty{jj},...
                'Marker',mark_typ{jj},...
                'LineWidth',1.5);
        end
    end
    ax=gca;
    xlabel('$\mathcal{K}_i$');
    if ~flag_analysis_3d
        ylabel('Normalised correlation $\bar{g}^{(2)}_{\textrm{BB}}(0)$');
    else
        ylabel('Correlation amplitude $\mathcal{G}$');      % fitted
    end
    box on;
    ax.FontSize=12;
    ax.LineWidth=1.2;
    
    str_tau=sprintf('%s%0.3g %s','$\tau=$',tT,'$\mu$s');
%     text('Units','normalized','Position',[0.95 0.9],'String',str_tau,...
%         'FontSize',14,'VerticalAlignment','middle','HorizontalAlignment','right');
    text('Units','normalized','Position',[0.5 1.03],'String',str_tau,...
        'FontSize',14,'VerticalAlignment','baseline','HorizontalAlignment','center');
end

%% save data
if exist('flag_save_data','var')   % otherwise must be called manually
    if flag_save_data
        %% Manually run this section
        path_save_dir='C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\figS\k_calibration';
        
%         % workspace       % VERY LARGE
%         fname=sprintf('data_%d_%s',II,getdatetimestr);
%         save(fullfile(path_save_dir,fname));
        
        % figs
        for ii=1:nparam
            tT=1e6*par_T(ii);
            fname=sprintf('fig_%0.2g_%s',tT,getdatetimestr);
            fpath=fullfile(path_save_dir,fname);
            
            saveas(hfig{ii},strcat(fpath,'.fig'),'fig');
            print(hfig{ii},strcat(fpath,'.svg'),'-dsvg');
        end
    end
end