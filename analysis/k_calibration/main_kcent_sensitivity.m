%% Halo centering sensitivity analysis
%   
% DKS
% 20181015

%% load data
if ~exist('flag_load_override','var') || ~flag_load_override
    path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot';
    fname_mat='exp5_1_20180813_1.mat';
    path_data=fullfile(path_dir,fname_mat);
    
    load(path_data);
end

%% define g2 analysis
dk_bb=0.04;     % back-to-back pair correlation length (rms width) [norm unit]
dk0={dk_bb*[-1,1],dk_bb*[-1,1],dk_bb*[-1,1]};   % define BB pair single-pixel

lim_Dk=[-dk_bb,dk_bb];          % limits for shift [normed]
n_Dk=15;                        % n-points to scan (MUST BE ODD!)

Dk_0=linspace(lim_Dk(1),lim_Dk(2),n_Dk);     % shift in halo center in each dim (symm)
idx_0=find(Dk_0==0);        % index to Dk=0

%% evaluate g2 while shifting center
g2_shift=cell(nparam,1);
for ii=1:nparam
    tk=k_par{ii};
    g2_shift{ii}=cell(1,3);
    for tt=1:3
        g2_shift{ii}{tt}=NaN(n_Dk,3);
    end
    for jj=1:3
        for kk=1:n_Dk
            tDk=circshift([Dk_0(kk),0,0],jj-1);     % get new origin shift
            
            % evaluate g2 at near BB
            tg2_11=g2_bb(boost_zxy(tk(:,1),tDk),dk0);       % (up,up)
            tg2_00=g2_bb(boost_zxy(tk(:,2),tDk),dk0);       % (down,down)
            
            ttk=tk;
            ttk(:,1)=boost_zxy(tk(:,1),tDk);
            tg2_01=g2x_bb(ttk,dk0);           % (up,down): shift only UP
            
            % store
            g2_shift{ii}{1}(kk,jj)=tg2_11;
            g2_shift{ii}{2}(kk,jj)=tg2_00;
            g2_shift{ii}{3}(kk,jj)=tg2_01;
        end
    end
end

% evaluate g2 normed wrt Dk=0
g2_shift_norm=cell(size(g2_shift));
for ii=1:nparam
    g2_shift_norm{ii}=cellfun(@(x) x./x(idx_0,:),g2_shift{ii},'UniformOutput',false);
end

%% data vis
[c,cl,cd]=palette(3);
line_sty={'-','--',':'};
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
            plot(Dk_0,g2_shift_norm{ii}{jj}(:,kk),...
                'Color',c(kk,:),'LineStyle',line_sty{jj},...
                'LineWidth',1.5);
        end
    end
    ax=gca;
    xlabel('$\mathcal{K}_i$');
    ylabel('$g^{(2)}_{\textrm{BB}}(0)$');
    box on;
    ax.FontSize=12;
    ax.LineWidth=1.2;
end

%% save data
if exist('flag_save_data','var')   % otherwise must be called manually
    if flag_save_data
        %% Manually run this section
        path_save_dir='C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\figS\k_calibration';
        
        for ii=1:nparam
            tT=1e6*par_T(ii);
            fname=sprintf('fig_%0.2g_%s',tT,getdatetimestr);
            fpath=fullfile(path_save_dir,fname);
            
            saveas(hfig{ii},strcat(fpath,'.fig'),'fig');
            print(hfig{ii},strcat(fpath,'.svg'),'-dsvg');
        end
    end
end