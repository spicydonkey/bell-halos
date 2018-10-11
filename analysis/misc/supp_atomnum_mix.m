%% Atom number mixing from spin,k processed halo
% DKS
% 20181010

%% load data
% get spin-momentum processed data
path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\bell_epr_2018\proc\exp3_yrot_char';
fname_data='exp3_s1_20180823_1.mat';
path_data=fullfile(path_dir,fname_data);

load(path_data);        % load data

%% process data
K=k_par;
n_shot_par=cellfun(@(x) size(x,1),K);

N=cell(nparam,1);
for ii=1:nparam
    N{ii}=cellfun(@(x) size(x,1),K{ii});
end
Navg=cellfun(@(x) mean(x,1),N,'UniformOutput',false);
Navg=cat(1,Navg{:});
Nstd=cellfun(@(x) std(x,0,1),N,'UniformOutput',false);
Nstd=cat(1,Nstd{:});
Nse=Nstd./sqrt(n_shot_par)';

% NUMBER LOSS
Ntot=cellfun(@(x) sum(x,2),N,'UniformOutput',false);
Ntot_avg=cellfun(@(x) mean(x), Ntot);
Ntot_std=cellfun(@(x) std(x), Ntot);
Ntot_se=Ntot_std./sqrt(n_shot_par)';
Nloss=Ntot_avg(1)-Ntot_avg;     % decrease in counts from NO PULSE

%%% collate loss term to mJ
Navg(:,4)=Nloss;
Nstd(:,4)=NaN(size(Nstd,1),1);
Nse(:,4)=NaN(size(Nse,1),1);


%%% smooth fit to data
polyfit_N=cell(3,1);
ninterp=1e4;
tau_sm_fit=linspace(min(par_T),max(par_T),ninterp);
Nfit=NaN(3,ninterp);
for ii=1:4
    polyfit_N{ii}=polyfit(par_T,Navg(:,ii),6);
    Nfit(ii,:)=polyval(polyfit_N{ii},tau_sm_fit);
end

%% vis
% grahics configs
[c0,clight,cdark]=palette(4);
mark_ss={'o','s','^','d'};
mark_siz=8;
line_wid=2;
str_ss={'$m_J=1$','$m_J=0$','$m_J=-1$','loss'};

hfig=figure('Units', 'normalized', 'Position', [0.2,0.2,0.3,0.3]);
hfig.Renderer='painters';
hold on;
pleg=NaN(1,4);      % data for legend
for ii=1:4
    tp=myploterr(1e6*par_T,Navg(:,ii),[],Nse(:,ii),mark_ss{ii},c0(ii,:));  
    set(tp(1),'MarkerSize',mark_siz,'LineWidth',line_wid,'DisplayName',str_ss{ii});
    pleg(ii)=tp(1);     % line data to show in legend
    set(tp(2),'LineWidth',line_wid,'DisplayName','');
    
    % fit
    tpfit=plot(1e6*tau_sm_fit,Nfit(ii,:),'LineWidth',2.2,'Color',0.5*(c0(ii,:)+clight(ii,:)));
    uistack(tpfit,'bottom');
end

ax=gca;
axis tight;
lgd=legend(pleg);
lgd.FontSize=11;
lgd.Location='East';
% lgd.Title.String='$m_J$';
box on;
xlabel('Pulse duration $\tau$ [$\mu$s]');
ylabel('Atom counts');
ax.FontSize=11;
ax.LineWidth=1.2;