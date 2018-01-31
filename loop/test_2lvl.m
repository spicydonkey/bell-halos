% Two-level system approximation of mF=0,1 atoms in the Raman process
%
% DK Shin
%


%% load data
% load k-space point-cloud data
load('k_par_20180130.mat');

% load experimental params
load('param_20180131.mat');


%% Distribution of scattered atoms in He* triplet

% number of *detected* and *captured* scattered atoms
N=cellfun(@(k) shotSize(k),k_par,'UniformOutput',false);
N_nshot=cellfun(@(k) size(k,1),k_par,'UniformOutput',false);  % number of shots repd per param

N_mean=cellfun(@(n) mean(n,1),N,'UniformOutput',false);
N_std=cellfun(@(n) std(n,0,1),N,'UniformOutput',false);
N_se=cellfun(@(s,n) s/sqrt(n),N_std,N_nshot,'UniformOutput',false);
% tidy data format
temp_N_mean=cell(1,2);
temp_N_std=cell(1,2);
temp_N_se=cell(1,2);
for ii=1:2
    temp_N_mean{ii}=vertcat(N_mean{ii,:});
    temp_N_std{ii}=vertcat(N_std{ii,:});
    temp_N_se{ii}=vertcat(N_se{ii,:});
end
N_mean=temp_N_mean;
N_std=temp_N_std;
N_se=temp_N_se;


% plot
data_name={'1','0','-1'};
data_linestyle={'-','--','-.'};
data_color=distinguishable_colors(3);

for ii=1:2
    figname=sprintf('src_%d',mf(ii));
    figure('Name',figname);
    for jj=1:3
        hold on;
        plot(1e6*traman,N_mean{ii}(:,jj),...
            'Color',data_color(jj,:),...
            'LineStyle',data_linestyle{jj},...
            'LineWidth',1.5,...
            'DisplayName',data_name{jj});
%         h=ploterr(traman,N_mean{ii}(:,jj),[],N_se{ii}(:,jj),'o','hhxy',0);
    end
    xlabel('$\tau$ ($\mu$s)');
    ylabel('$N$');
    lgd=legend('show');
    title(lgd,'$m_F$')
    box on;
end


%% loss to mf=-1

% loss into -1 state zero-referenced to non-ideal source
N_loss=cellfun(@(n) n(:,3)-n(1,3),N_mean,'UniformOutput',false);

N_tot=cellfun(@(n) sum(n,2),N_mean,'UniformOutput',false);

r_loss=cellfun(@(n1,n2) n1./n2,N_loss,N_tot,'UniformOutput',false);


% plot
data_name={'0','1'};
data_linestyle={'-','--'};
data_color=distinguishable_colors(2);

figure;
for ii=1:2
    hold on;
    plot(1e6*traman,r_loss{ii},...
        'Color',data_color(ii,:),...
        'LineStyle',data_linestyle{ii},...
        'LineWidth',1.5,...
        'DisplayName',data_name{ii});
        
    xlabel('$\tau$ ($\mu$s)');
    ylabel('$\eta_{loss}$');
    lgd=legend('show');
    title(lgd,'source $m_F$')
    box on;
end

