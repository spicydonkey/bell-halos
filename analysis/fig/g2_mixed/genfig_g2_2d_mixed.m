%% Display various g2 functions under unitary rotation (mixing spins)
% DKS
% 20181008

% config
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp5_bell_yrot\bell_signal';
data_fname='bell_20181009_bs_fix.mat';
path_data=fullfile(path_dir,data_fname);

% load data
load(path_data);    % load into workspace

%%% config
cmap='viridis';


%% g2 plotting
n_param_data=sum(cellfun(@(x) numel(x),{S.par_T}));
hfig=cell(n_param_data,1);
counter=0;
for ii=1:numel(S)
    tS=S(ii);
    tn_config=numel(S(ii).par_T);
    for jj=1:tn_config
        counter=counter+1;
        %% g2 2D profile
        tdk=unique(tS.dk{1}{1});    % dk grid (symmetric) in each dim
        
        tg2=tS.g2{jj};
        tg2_se=tS.g2_sdev{jj};
        I0=tS.idx_dk0;        
        
        tg2_max=max(cellfun(@(x) max(max(x(:,:,I0))),tg2));
        tclims=[0,tg2_max];
        
        hfig{counter}=figure('Units', 'normalized', 'Position', [0.2,0.2,0.45,0.26]);
        hfig{counter}.Renderer='painters';
        
        bgAxes=axes('Position',[0,0,1,1],'XColor','none','YColor','none',...
            'XLim',[0,1],'YLim',[0,1]);     % bg axes
      
        % positions
        x=linspace(0.1,0.9,4);
        w=0.9*diff(x(1:2));
        y=linspace(0.17,1.05,2);
        h=0.8*diff(y(1:2));
        
        %%% labels
        % param
        str_param=sprintf('%s%0.3g %s','$\tau=$',1e6*tS.par_T(jj),'$\mu$s');
        text(0.015,y(1)+0.5*h,str_param,'FontSize',18,'Rotation',90,...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontWeight','bold');
        
        % corr type
        str_corr_type={'$\uparrow\uparrow$','$\downarrow\downarrow$','$\uparrow\downarrow$'};
        for kk=1:3
            text(x(kk)+0.5*diff(x(1:2)),y(2)-0.1,str_corr_type{kk},...
                'FontSize',18,'FontWeight','bold','HorizontalAlignment','center',...
                'VerticalAlignment','middle');
        end
        
        % colorbar
        str_cbar='$g^{(2)}_{\textrm{BB}}$';
        tt=text(0.975,y(1)+0.5*h,str_cbar,'FontSize',16,'Rotation',90,...
            'HorizontalAlignment','center','VerticalAlignment','middle');
        
        %%% data
        for kk=1:3
%             subplot(1,3,kk);
            axes('Position',[x(kk),y(1),w,h]);      % set up axes

            timg=imagesc(tdk,tdk,tg2{kk}(:,:,I0),tclims);
            
            ax=gca;
            ax.YDir='normal';       % conventional 2D axes direction
            ax.XTick=-0.2:0.1:0.2; ax.YTick=-0.2:0.1:0.2;   % toggle axis ticks off
            ax.XTickLabel=[]; ax.YTickLabel=[];
            ax.TickLength=0.033*ones(1,2);
            ax.LineWidth=1.5;
            ax.XColor='k';  ax.YColor='k';
            ax.FontSize=12;
            colormap(cmap);
            
            % axes label
            if kk==1
                ax.YTickLabelMode='auto';
                ylabel('$\Delta \mathrm{k}_x$');
            end
            ax.XTickLabelMode='auto';
            if kk==2
                ax.XTickLabelMode='auto';
                xlabel('$\Delta \mathrm{k}_z$');
            end
        end
        
        % separate colorbar
        axes('Position',[x(kk)+0.07,y(1),w,h]);         % same as rightmost + dX
        timg=imagesc(tdk,tdk,tg2{kk}(:,:,I0),tclims);   
        timg.Visible='off';     % hide it
        ax=gca;
        ax.Visible='off';
        cbar=colorbar('Position',[0.89,y(1),0.02,h],'FontSize',14);
    end
end

%% save data
if false        % must be called manually
    %% Manually run this section
    path_save_dir='C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\figS\g2_2d_mix';
    
    for ii=1:n_param_data
        fname=sprintf('fig_%d_%s',ii,getdatetimestr);
        fpath=fullfile(path_save_dir,fname);
        
%         saveas(hfig{ii},fpath,'fig');     % HUGE!
        print(hfig{ii},fpath,'-dsvg');
    end
end